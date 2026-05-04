"""
generate_oracle.py generates trajectories for each threshold and computes value functions as an oracle
=============================================================================
Each individual i is observed for N_DAYS=4 days. T_i = number of eating
events they experience (random, expected ~15). The oracle computes discounted
sums truncated at each of DAY_CUTOFFS = [1, 2, 4] days:

  VY(gamma, d) = mean_i  sum_{t=0}^{T_i(d)-1} gamma^t * Y_{i,t+1}
  VA(gamma, d) = mean_i  sum_{t=0}^{T_i(d)-1} gamma^t * A_{i,t}

where T_i(d) = number of events for individual i that fired within the
first d days (= first d * STEPS_PER_DAY outer steps).

Output: oracle_c{int(c)}_{n_seeds}seeds.csv
  columns: gamma, days, VY, VA  (9 rows = 3 gammas x 3 day cutoffs)

*** CHANGE ONLY THIS BLOCK TO ALTER THE DESIGN ***
  N_DAYS      = 8             # simulation window per patient
  DAY_CUTOFFS = [1, 2, 4, 8]  # truncation horizons (must all be <= N_DAYS)
  GAMMAS      = [0, 1-1/24, 1-1/48, 1-1/96, 1-1/192, 1]  # effective horizons 0.5/1/2/4/Inf days
  THRESHOLDS  = [120.0, 125.0, ..., 150.0, 155.0, ..., 185.0]
"""

import argparse, os, time, warnings, concurrent.futures
from collections import defaultdict
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

# ── Fixed simulation constants ────────────────────────────────────────────────
OBS_MIN      = 30
N_INT        = OBS_MIN // 3          # 10 sub-steps of 3 min each
_TARGET      = 140.
ADULTS       = [f'adult#{i:03d}' for i in range(1, 11)]
_START_TIME  = datetime(2023, 1, 1, 5, 0, 0)

# ── CHANGE ONLY THIS BLOCK ───────────────────────────────────────────────────
N_DAYS      = 8             # observation window (in number of days)
DAY_CUTOFFS = [1, 2, 4, 8]  # milestone horizons, must all be <= N_DAYS
GAMMAS      = [0, *(1 - 1/(24 * x) for x in (1, 2, 4, 8)), 1]  # effective horizons 0, 0.5, 1, 2, 4, Inf days
THRESHOLDS  = list(range(120, 190, 5)) # 120, 130, ..., 185 (14 numbers)
# ─────────────────────────────────────────────────────────────────────────────

STEPS_PER_DAY = 24 * 60 // OBS_MIN   # = 48 outer steps per day
STEPS_MAX     = N_DAYS * STEPS_PER_DAY  # = 192

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
    Zs  : float32 (STEPS_MAX,)  pre-prandial CGM at step t if a meal fired,
                                 else -inf.
    As  : int8    (STEPS_MAX,)  A_t = 1{Z_t >= c} at meal steps, 0 otherwise.
    Ys  : float32 (STEPS_MAX,)  Y_t = _risk_reward(CGM_{t+1}), i.e. the
                                 instantaneous reward based on the last CGM
                                 observation of outer step t (= "next-period"
                                 CGM in 30-min units).

    Discounting is gamma^t with t the outer-step index (not meal-event index).
    For day cutoff d:
      T_d   = d * STEPS_PER_DAY   (no searchsorted needed)
      disc  = gamma ** arange(T_d)
      VY    = disc @ Ys[:T_d]
      VA    = disc @ As[:T_d].astype(float)
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
        last_cgm   = float(obs.CGM)   # fallback if done fires before first sub-step

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

        # Y_t = reward based on next-period CGM (last sub-step of this outer step)
        Ys[step_idx] = _risk_reward(last_cgm)

        if meal_fired:
            Zs[step_idx] = Z_onset
            A = int(Z_onset >= c)
            As[step_idx] = A
            pending_corr = max(0., (Z_onset - _TARGET) / CF) if A else 0.

        if done:
            break

    return Zs, As, Ys

def _oracle_worker(args):
    pname, seed, c, pcache = args
    Zs, As, Ys = simulate(pname, seed, c, pcache)
    meal_mask = np.isfinite(Zs)   # True at meal steps, False elsewhere
    out = {}
    for g in GAMMAS:
        for d in DAY_CUTOFFS:
            T_d = d * STEPS_PER_DAY
            if abs(g) < 1e-12:
                meal_steps = np.where(meal_mask[:T_d])[0]
                if len(meal_steps) == 0:
                    out[(g, d)] = (0., 0.)    
                else:
                    first = meal_steps[0]
                    out[(g, d)] = (float(Ys[first]), float(As[first]))
            else:
                disc = g ** np.arange(T_d)
                VY   = float(disc @ Ys[:T_d])
                VA   = float(disc @ As[:T_d].astype(float))
                out[(g, d)] = (VY, VA)
    return out


def _save(acc_Y, acc_A, n_units_done, c, out_dir):
    rows = []
    for g in GAMMAS:
        for d in DAY_CUTOFFS:
            rows.append({
                'gamma': round(g, 10),
                'days':  d,
                'VY':    acc_Y[(g, d)] / n_units_done,
                'VA':    acc_A[(g, d)] / n_units_done,
            })
    df = pd.DataFrame(rows)
    n_seeds = n_units_done // len(ADULTS)
    fname = out_dir / f'oracle_c{int(c)}_{n_seeds}seeds.csv'
    df.to_csv(fname, index=False)
    # Diagnostic: middle gamma, all day cutoffs
    mid_g = GAMMAS[len(GAMMAS) // 2]
    sub = df[df['gamma'].sub(mid_g).abs() < 1e-9]
    for _, r in sub.iterrows():
        print(f"  {fname.name}  "
              f"[g={'0' if mid_g == 0 else f'1-1/{round(1/(1-mid_g))}'}, d={int(r.days)}d: "
              f"VY={r.VY:.4f}  VA={r.VA:.4f}]")


def run(job_id, max_seeds, save_at, out_dir, max_workers):
    c       = THRESHOLDS[job_id]
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    pcache  = _build_pcache()

    n_units_total = len(ADULTS) * max_seeds
    save_at_units = [len(ADULTS) * s for s in save_at]

    print(f"Oracle v6  job {job_id}: c={c}  {max_seeds} seeds ({n_units_total} units)")
    print(f"  N_DAYS={N_DAYS}  STEPS_MAX={STEPS_MAX}  STEPS_PER_DAY={STEPS_PER_DAY}")
    print(f"  DAY_CUTOFFS={DAY_CUTOFFS}")
    print(f"  GAMMAS={[round(g,6) for g in GAMMAS]}  "
          f"eff horizons ~{[round(1/((1-g)*STEPS_PER_DAY),1) if g < 1 else 'inf' for g in GAMMAS]} days")
    print(f"  Save at seeds: {save_at}  workers={max_workers}")

    acc_Y = defaultdict(float)   # keyed by (gamma, days)
    acc_A = defaultdict(float)

    args = [(ADULTS[i % len(ADULTS)], i // len(ADULTS), c, pcache)
            for i in range(n_units_total)]

    t0 = time.time(); n_done = 0; next_save_idx = 0

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as ex:
        cs = max(1, n_units_total // (max_workers * 8))
        for result in ex.map(_oracle_worker, args, chunksize=cs):
            for (g, d), (vy, va) in result.items():
                acc_Y[(g, d)] += vy
                acc_A[(g, d)] += va
            n_done += 1

            if (next_save_idx < len(save_at_units) and
                    n_done == save_at_units[next_save_idx]):
                elapsed = time.time() - t0
                print(f"  checkpoint {n_done}/{n_units_total}  {elapsed:.0f}s")
                _save(acc_Y, acc_A, n_done, c, out_dir)
                next_save_idx += 1

            if n_done % max(1, n_units_total // 20) == 0:
                elapsed = time.time() - t0
                eta = elapsed / n_done * (n_units_total - n_done)
                print(f"  {n_done:>7}/{n_units_total}  {elapsed:.0f}s  ETA {eta:.0f}s")

    if n_done not in save_at_units:
        _save(acc_Y, acc_A, n_done, c, out_dir)

    total = time.time() - t0
    print(f"Done: job {job_id} (c={c})  {total:.1f}s ({total/60:.1f}min)")


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--job_id',      type=int,
                    default=int(os.environ.get('SLURM_ARRAY_TASK_ID', 1)))
    ap.add_argument('--max_seeds',   type=int, default=10000)
    ap.add_argument('--save_at',     type=str, default='1000,5000,10000')
    ap.add_argument('--out_dir',     type=str, default='runs/oracle')
    ap.add_argument('--max_workers', type=int,
                    default=int(os.environ.get('SLURM_CPUS_PER_TASK', 32)))
    a = ap.parse_args()
    save_at = [int(x) for x in a.save_at.split(',')]
    run(a.job_id, a.max_seeds, save_at, a.out_dir, a.max_workers)
