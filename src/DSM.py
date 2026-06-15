r'''
Model of Disease state from Sagar
---

x(t): Hidden disease state
    x' = a*A(t)(1-x) - (b + cs*Ts(t) + cl*Tl(t))*x

g(t): Observable GI symptom burden
    g' = k*x*(1-g) - (r + qs*Ts(t))*g

d(t): Long term chronic degredation
    d' = \eps*[m*(g - g_c)_+(1-d) - nl*Tl(t)*d]

A(t): Adverse load from observed data (fit from data)
    A(t) = w_1 * BP_d(t) + w_2 * HR_d(t) + w_3 * BM_d(t) +
         + w_4 * W_d(t) + w_5 * C

Treatment is split into a short-term and a long-term component:

Ts(t): Short-term / acute treatment, acts on active disease (cs*Ts term) and
       on symptom recovery (qs*Ts term)
Tl(t): Long-term / maintenance treatment + adherence, acts on active disease
       (cl*Tl term) and on chronic degradation (nl*Tl term)
    Ts(t), Tl(t) \in [0, 1]
    where 0 -> poor adherence
          1 -> effective treatment / strong adherence

Model Paramters:
a   -> effect of adverse load on disease activity
b   -> natural recovery rate
cs  -> short-term treatment effect on active disease
cl  -> long-term treatment effect on active disease
qs  -> short-term treatment effect on symptoms
k   -> symptom production rate from disease burden
r   -> symptom recovery rate
m   -> flare effect on long-term degredation
nl   -> long-term treatment effect on long term degredation
eps -> slow progression rate
gc  -> flare threshold


'''

import numpy as np

def ode(t, y, a, b, cs, cl, qs, k, r, eps, m, gc, nl, A, Ts, Tl):
    x = y[0]
    g = y[1]
    d = y[2]

    # A, Ts and Tl may be passed as constants or as callables of time.
    A_t  = A(t)  if callable(A)  else A
    Ts_t = Ts(t) if callable(Ts) else Ts
    Tl_t = Tl(t) if callable(Tl) else Tl

    dx = (a*A_t*(1-x) - (b + cs*Ts_t + cl*Tl_t)*x)
    dg = (k*x*(1-g) - (r + qs*Ts_t)*g)
    dd = (eps*(m*max(0, g-gc)*(1-d) - nl*Tl_t*d))

    z = np.zeros(3,)
    z[0] = dx
    z[1] = dg
    z[2] = dd

    return z

'''
Helper Functions

'''

# Healthy baselines and characteristic scales for each observed signal.
# Scales are used to non-dimensionalize each signal in the adverse load.
BASELINES = {"BP": 120.0, "HR": 70.0, "BM": 1.5, "W": 70.0}
SCALES    = {"BP": 15.0,  "HR": 10.0, "BM": 1.5, "W": 5.0}

# Default weights for the adverse load.  All strictly positive and summing to 1
# (equal weighting) for now; a direction-aware measure for each signal (e.g.
# adverse weight loss) will be handled separately later.  Comorbidity term
# (w5 * C) is omitted for now.
WEIGHTS = {"BP": 0.25, "HR": 0.25, "BM": 0.25, "W": 0.25}


def make_observations(t, baselines=None, flare_times=(), flare_width=0.75,
                      flare_mag=None, noise=0.0, seed=None):
    '''Synthetic time series for the four observed signals on the grid `t`.

    Each signal is built as
        baseline + (flare bumps) + (gaussian measurement noise).

    During a GI flare, heart rate and bowel movements rise, blood pressure
    rises mildly, and weight drops.  `flare_times` is a list of bump centers.

    Returns a dict of arrays keyed "BP", "HR", "BM", "W", aligned to `t`.
    '''
    t = np.asarray(t, float)
    baselines = {**BASELINES, **(baselines or {})}
    rng = np.random.default_rng(seed)

    default_mag = {"BP": 8.0, "HR": 15.0, "BM": 3.0, "W": -3.0}
    flare_mag = {**default_mag, **(flare_mag or {})}

    obs = {}
    for key in ("BP", "HR", "BM", "W"):
        sig = np.full_like(t, baselines[key])
        for tf in flare_times:
            sig = sig + flare_mag[key]*np.exp(-0.5*((t - tf)/flare_width)**2)
        if noise:
            sig = sig + rng.normal(0.0, noise*SCALES[key], size=t.shape)
        obs[key] = sig
    return obs


def adverse_load(obs, weights=None, baselines=None, scales=None):
    '''Combine the observed signals into a dimensionless adverse load A.

    Each signal is taken as a deviation from baseline, non-dimensionalized by
    its scale, then weighted and summed:

        A = sum_i w_i * (s_i - baseline_i) / scale_i

    Weight magnitudes are normalized to sum to 1 so A stays an order-1,
    interpretable weighted average regardless of the raw weights passed in.
    The comorbidity term (w5 * C) is omitted for now.
    '''
    weights   = {**WEIGHTS,   **(weights or {})}
    baselines = {**BASELINES, **(baselines or {})}
    scales    = {**SCALES,    **(scales or {})}

    total = sum(abs(weights[key]) for key in ("BP", "HR", "BM", "W"))
    if total:
        weights = {key: weights[key]/total for key in ("BP", "HR", "BM", "W")}

    A = np.zeros_like(np.asarray(next(iter(obs.values())), float))
    for key in ("BP", "HR", "BM", "W"):
        A = A + weights[key]*(np.asarray(obs[key], float) - baselines[key])/scales[key]
    return A


def as_function(t, values):
    '''Wrap sampled `values` on grid `t` into a callable f(tau) (linear interp).

    Lets a time series (e.g. A) be passed to `ode`/`solve_ivp`, which evaluate
    at arbitrary times.  Values are held constant outside the grid endpoints.
    '''
    t = np.asarray(t, float)
    values = np.asarray(values, float)

    def f(tau):
        return np.interp(tau, t, values)
    return f


def treatment_short(t, courses=(), level=1.0):
    '''Short-term / acute treatment Ts(t) in [0, 1].

    Modeled as discrete acute courses, each given as a (start, duration) pair;
    Ts = `level` while a course is active and 0 otherwise.  Acts on active
    disease (cs*Ts) and on symptom recovery (qs*Ts).

    Accepts a scalar or array `t`.
    '''
    t = np.asarray(t, float)
    Ts = np.zeros_like(t)
    for start, dur in courses:
        Ts = np.where((t >= start) & (t < start + dur), level, Ts)
    return Ts


def treatment_long(t, level=0.8, amp=0.0, period=7.0, phase=0.0):
    '''Long-term / maintenance treatment + adherence Tl(t) in [0, 1].

    Baseline adherence `level` with an optional periodic fluctuation
    (amplitude `amp`, period `period`) to mimic imperfect long-term adherence.
    Acts on active disease (cl*Tl) and on chronic degradation (nl*Tl).

    Accepts a scalar or array `t`.
    '''
    t = np.asarray(t, float)
    Tl = level + amp*np.sin(2*np.pi*(t - phase)/period)
    return np.clip(Tl, 0.0, 1.0)


'''
Output Functions

'''
