"""
generate_data.py  Data generation for simglucose experiments
================================================================================
Simulates N_DAYS=4 days under policy c=C_REF=150.
Decision rule: at each eating event t, A_t = 1{Z_t >= C_REF},
where Z_t is the CGM at the exact meal-onset sub-step (pre-prandial).

Reward: negative Kovatchev risk level at obs.CGM instead of step.reward
Y_{t+1} is normalized by substep count N_t (average risk per substep)

Output: rep{NNNN}.parquet with columns
  i          : unit index within rep (0 .. 10*n_seeds_per_rep - 1)
  t          : eating event index, 0-based (not clock time)
  Z          : pre-prandial CGM at event t  (= Z_t in oracle notation)
  A          : action at event t  (= A_t = 1{Z >= C_REF})
  Y          : normalized reward Y_{t+1} = mean substep risk from event t
               to event t+1 (attributable to A_t; = Ys[t+1] in oracle notation)
  outer_step : outer 30-min step index when event t fired
               (used for day-cutoff filtering: outer_step < d * 48 for d days)

Seed design (no collision with oracle seeds 0 .. 9,999,999):
  seed = DATA_BASE + rep_global * STRIDE + unit_i
  DATA_BASE = 10,000,000,  STRIDE = 10,000
"""

import argparse, os, time, warnings, concurrent.futures
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from simglucose.simulation.env import T1DSimEnv
from simglucose.controller.base import Action
from simglucose.controller.basal_bolus_ctrller import CONTROL_QUEST, PATIENT_PARA_FILE
from simglucose.sensor.cgm import CGMSensor
from simglucose.actuator.pump import InsulinPump
from simglucose.patient.t1dpatient import T1DPatient
from simglucose.simulation.scenario_gen import RandomScenario

warnings.filterwarnings('ignore')

# ── Fixed constants ────────────────────────────────────────────────────────────
OBS_MIN      = 30
N_INT        = OBS_MIN // 3         # 10 sub-steps of 3 min each
_TARGET      = 140.
N_DAYS       = 4
STEPS_PER_DAY = 24 * 60 // OBS_MIN  # = 48
STEPS_MAX    = N_DAYS * STEPS_PER_DAY  # = 192 (4 days) 384 (8 days)
C_REF        = 150.  # data-generating policy threshold
ADULTS       = [f'adult#{i:03d}' for i in range(1, 11)]
_START_TIME  = datetime(2023, 1, 1, 5, 0, 0)

# ── Seed isolation from oracle ────────────────────────────────────────────────
STRIDE     = 10_000
DATA_BASE  = 10_000_000

_PNAMES = ['BW','u2ss','kmax','kmin','kabs','f','b','d','kp1','kp2','kp3',
           'Fsnc','ke1','ke2','k1','k2','Vm0','Vmx','Km0','m1','m2','m4',
           'm30','ka1','ka2','kd','Vi','Ib','p2u','ki','ksc','Vg']


class _FP:
    def __init__(self, s):
        object.__setattr__(self, '_s', s)
        for k in _PNAMES:
            try: object.__setattr__(self, k, float(s[k]))
            except: pass
    def __getattr__(self, n): return getattr(object.__getattribute__(self, '_s'), n)
    def __getitem__(self, k):
        try: return object.__getattribute__(self, k)
        except AttributeError: return object.__getattribute__(self, '_s')[k]
    def items(self): return object.__getattribute__(self, '_s').items()


def _build_pcache():
    q = pd.read_csv(CONTROL_QUEST)
    p = pd.read_csv(PATIENT_PARA_FILE)
    cache = {}
    for pname in ADULTS:
        qr = q[q.Name.str.match(pname)].iloc[0]
        pr = p[p.Name.str.match(pname)].iloc[0]
        cache[pname] = dict(basal=float(pr.u2ss * pr.BW / 6000),
                            CF=float(qr.CF), CR=float(qr.CR))
    return cache


def _risk_reward(BG_last, **kwargs):
    """
    Negative Kovatchev instantaneous risk level at BG_last (mg/dL).

    Formula: f(BG) = 1.509 * (ln(BG)^1.084 - 5.381)
             r(BG) = 10 * f(BG)^2
    Reward  = -r(BG)

    Source: Kovatchev et al. (1997), Diabetes Care 20(11):1655-1658.
    Constants alpha=1.084, beta=5.381, gamma=1.509 are for BG in mg/dL.
    """
    fBG = 1.509 * (np.log(max(BG_last, 1.)) ** 1.084 - 5.381)
    return -10. * fBG ** 2


def simulate(pname, seed, c, pcache):
    """
    Simulate one N_DAYS-day trajectory.

    Returns
    -------
    Zs  : float32 (STEPS_MAX,)  pre-prandial CGM at step t if meal fired, else -inf.
    As  : int8    (STEPS_MAX,)  A_t = 1{Z_t >= c} at meal steps, 0 otherwise.
    Ys  : float32 (STEPS_MAX,)  Y_t = _risk_reward(CGM_{t+1}), the instantaneous
                                 reward based on the last CGM of outer step t.
    """
    cache = pcache[pname]
    basal = cache['basal'];  CF = cache['CF'];  CR = cache['CR']

    p    = T1DPatient.withName(pname, random_init_bg=True, seed=seed)
    p._params = _FP(p._params)
    s    = CGMSensor.withName('Dexcom', seed=seed)
    pump = InsulinPump.withName('Insulet')
    env  = T1DSimEnv(p, s, pump,
                     RandomScenario(start_time=_START_TIME, seed=seed))

    step0 = env.reset()
    obs   = step0.observation
    done  = step0.done
    meal  = step0.info['meal']

    Zs = np.full(STEPS_MAX, -np.inf, dtype=np.float32)
    As = np.zeros(STEPS_MAX, dtype=np.int8)
    Ys = np.zeros(STEPS_MAX, dtype=np.float32)

    pending_corr = 0.

    for step_idx in range(STEPS_MAX):
        meal_fired = False
        Z_onset    = None
        last_cgm   = float(obs.CGM)

        for k in range(N_INT):
            bb    = meal / (CR * 3.) if meal > 0 else 0.
            bolus = (pending_corr / 3. + bb) if k == 0 else bb
            if k == 0:
                pending_corr = 0.

            step     = env.step(Action(basal=basal, bolus=bolus))
            obs      = step.observation
            done     = step.done
            new_meal = step.info['meal']
            last_cgm = float(obs.CGM)

            if meal == 0 and new_meal > 0 and not meal_fired:
                meal_fired = True
                Z_onset    = np.float32(obs.CGM).item()

            meal = new_meal
            if done:
                break

        Ys[step_idx] = _risk_reward(last_cgm)

        if meal_fired:
            Zs[step_idx] = Z_onset
            A = int(Z_onset >= c)
            As[step_idx] = A
            pending_corr = max(0., (Z_onset - _TARGET) / CF) if A else 0.

        if done:
            break

    return Zs, As, Ys

def _data_worker(args):
    i_unit, pname, seed, c, pcache = args
    Zs, As, Ys = simulate(pname, seed, c, pcache)
    return i_unit, Zs, As, Ys


def run(job_id, n_reps_per_job, n_seeds_per_rep, out_dir, max_workers):
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    pcache  = _build_pcache()
    n_units = len(ADULTS) * n_seeds_per_rep

    print(f"Job {job_id}: {n_reps_per_job} reps x {n_units} units  -> {out_dir}")
    print(f"  N_DAYS={N_DAYS}  STEPS_MAX={STEPS_MAX}  C_REF={C_REF}  "
          f"E[events/unit]={N_DAYS*3.75:.1f}  workers={max_workers}")

    t_wall = time.time()

    for rep_local in range(n_reps_per_job):
        rep_global = job_id * n_reps_per_job + rep_local
        base_seed  = DATA_BASE + rep_global * STRIDE
        t0 = time.time()

        dargs = [(i, ADULTS[i % len(ADULTS)], base_seed + i, C_REF, pcache)
                 for i in range(n_units)]

        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as ex:
            chunks = list(ex.map(_data_worker, dargs,
                                 chunksize=max(1, n_units // max_workers)))

        rows = []
        for i_unit, Zs, As, Ys in chunks:
            rows.append(pd.DataFrame({
                'i': np.full(STEPS_MAX, i_unit, np.int16),
                't': np.arange(STEPS_MAX,       dtype=np.int16),  # outer-step index
                'Z': Zs,   # -inf at non-meal steps
                'A': As,
                'Y': Ys,
            }))

        df = pd.concat(rows, ignore_index=True)

        # Sanity check: A_t = 1{Z_t >= C_REF} at only meal times
        meal_mask = np.isfinite(df['Z'].values)
        assert (df['A'].values[meal_mask] ==
                (df['Z'].values[meal_mask] >= C_REF).astype(np.int8)).all(), \
            "Action/threshold mismatch at meal steps — check C_REF and simulate()"
        assert (df['A'].values[~meal_mask] == 0).all(), \
            "Non-zero action at non-meal step — check simulate()"

        path = out_dir / f'rep{rep_global:04d}.parquet'
        df.to_parquet(path, compression='snappy', index=False)

        elapsed = time.time() - t_wall
        eta     = elapsed / (rep_local + 1) * (n_reps_per_job - rep_local - 1)
        print(f"  rep {rep_global:>4}: {len(df):>7,} rows  "
              f"E[T_i]={len(df)/n_units:.1f}  "
              f"A_rate={df['A'].mean():.3f}  "
              f"{time.time()-t0:.1f}s  ETA {eta:.0f}s")

    total = time.time() - t_wall
    print(f"Done: job {job_id}  {n_reps_per_job} reps  "
          f"{total:.1f}s ({total/60:.1f}min)")


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--job_id',          type=int,
                    default=int(os.environ.get('SLURM_ARRAY_TASK_ID', 0)))
    ap.add_argument('--n_reps_per_job',  type=int, default=100)
    ap.add_argument('--n_seeds_per_rep', type=int, default=50)
    ap.add_argument('--out_dir',         type=str, default='runs/data')
    ap.add_argument('--max_workers',     type=int,
                    default=int(os.environ.get('SLURM_CPUS_PER_TASK', 32)))
    a = ap.parse_args()
    run(a.job_id, a.n_reps_per_job, a.n_seeds_per_rep, a.out_dir, a.max_workers)
