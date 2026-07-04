#!/usr/bin/env python3
"""Certificate REPAIR for the White (arXiv:2201.05704) Section-5 dual-certificate
lower bound on mu, the Erdos minimum-overlap constant.

WHY A REPAIR IS NEEDED, AND WHAT THE CORRECT REPAIR IS
------------------------------------------------------
White's Prop. 9 states: Omega* = the *minimum* of the Section-5 convex program
is a lower bound on mu (any admissible f,M yields a feasible program point with
Omega = ||M||_inf, so the program minimum p* satisfies p* <= mu).

Consequently, to certify "mu >= X" we need a rigorous *lower* bound on the
program optimum p*.  A lower bound on the minimum of a convex program comes ONLY
from WEAK DUALITY: a dual-feasible point whose dual objective d satisfies
d <= p* <= mu.  White's own paper is explicit about this (Sec. 5.1 / Appendix II):
he certifies the bound with a *verified-feasible point of the DUAL SOCP*, never
with a primal point.

This is the crux the naive "repair" misses:

  * A PRIMAL-feasible point (even a strictly feasible one) has objective
    Omega_feas >= p*, i.e. it UPPER-bounds p*.  Combined with p* <= mu it says
    NOTHING about a lower bound on mu.  Shifting Omega down / scaling the primal
    variables into the feasible region therefore does NOT rigorously certify
    "mu >= X".  The interior-point solver's slightly-infeasible primal point
    (min slack -1.03e-7 in constraint 5.5 at N=2e5) is not the object to repair.

  * The object to repair is the DUAL point.  Our solve (cvxpy/CLARABEL) returns,
    for every constraint, a Lagrange multiplier.  These multipliers are the raw
    material of a weak-duality certificate.  Interior-point output is slightly
    dual-infeasible too; the repair pushes the dual point back to strict
    feasibility and quantifies the (tiny) objective penalty.

THE RIGOROUS BOUND WE COMPUTE (self-contained, interval-arithmetic-ready)
-------------------------------------------------------------------------
Let x = (Omega, w, v, c, d, eps, delta) and let

    B = { Omega in [0,1], w,v in [0,1]^N, |c|,|d| <= 2/pi,
          |eps| <= eps_bound, |delta| <= delta_bound }

be the box of the program's *simple* bound constraints.  B contains the whole
feasible region (0 <= w_j <= Omega <= 1, etc.), so for any multipliers
lambda >= 0 (equalities free), the Lagrangian L(x, lambda) = Omega +
sum_i lambda_i g_i(x) satisfies

    inf_{x in B} L(x, lambda)  <=  inf_{x feasible} L(x, lambda)  <=  p*  <=  mu,

the middle inequality because lambda_i g_i(x) <= 0 on the feasible set.

L(.,lambda) is convex, so for ANY reference point x* (we use the solve's own
primal vector) the tangent underestimator gives, coordinate-separably,

    inf_{x in B} L  >=  L(x*, lambda)
                        + sum_i  min_{x_i in B_i}  dL/dx_i(x*) * (x_i - x*_i).

Every term on the right is a finite float computation (a dot product, a gradient,
and per-coordinate box minima), so the whole bound is checkable in exact/interval
arithmetic starting from White's Lemma-backed constraint definitions.  This is
the number we certify:  mu >= rigorous_lower_bound.

REPAIR STEP.  The only thing we must guarantee about lambda is lambda_i >= 0 for
every inequality constraint (equality multipliers are free).  Interior-point
output can leave a multiplier a hair negative; we clip such entries to 0 (this is
the "repair": it keeps L a valid weak-duality lower bound) and report the total
clip mass.  The residual dual-infeasibility (how far the KKT stationarity
gradients sit from 0) is reported too -- it is the dual analogue of the primal
"min slack" and it is what governs how close rigorous_lower_bound sits to Omega*.

ATTRIBUTION.  The lower-bound METHOD (the whole convex program and its dual
verification) is White's (arXiv:2201.05704).  Our contribution here is (a) scaling
the single-box program to N=1.5e5-2e5, (b) extracting the full primal+dual
certificate, and (c) this weak-duality repair + rigorous re-verification.  A
bracket on mu is NOT a determination of mu.

Usage:
    python3 erdos_cert_repair.py /tmp/erdos_cert_N200000_R50_T50000.npz \
        [--out /tmp/erdos_repaired_cert_N200000.json] [--cvxpy-crosscheck]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def _recompute_alpha_beta(N: int, R: int, L: float):
    """alpha^-, alpha^+, beta^-, beta^+  of shape (2R, N), verbatim White Sec.5."""
    j = np.arange(1, N + 1, dtype=np.float64)
    m = np.arange(1, 2 * R + 1, dtype=np.float64)
    arg = np.outer(m, np.pi * L * (j - 0.5) / 2.0)
    slope = np.outer(m, np.full(N, np.pi * L / 4.0))
    cosv, sinv = np.cos(arg), np.sin(arg)
    return cosv - slope, cosv + slope, sinv - slope, sinv + slope


def repair_certificate(npz_path: str, cvxpy_crosscheck: bool = False) -> dict:
    z = np.load(npz_path, allow_pickle=True)
    N = int(z["N"]); R = int(z["R"]); T = int(z["T"]); L = float(z["L"])
    h1 = float(z["h1"]); h2 = float(z["h2"])
    p1 = float(z["p1"]); p2 = float(z["p2"]); q1 = float(z["q1"]); q2 = float(z["q2"])
    Omega_star = float(z["Omega"])
    w = np.asarray(z["w"], np.float64); v = np.asarray(z["v"], np.float64)
    c = np.asarray(z["c"], np.float64); d = np.asarray(z["d"], np.float64)
    eps = np.asarray(z["eps"], np.float64); delta = np.asarray(z["delta"], np.float64)
    a = np.asarray(z["a"], np.float64); b = np.asarray(z["b"], np.float64)
    eps_bound = np.asarray(z["eps_bound"], np.float64)
    delta_bound = np.asarray(z["delta_bound"], np.float64)

    am_, ap_, bm_, bp_ = _recompute_alpha_beta(N, R, L)  # (2R, N)
    m = np.arange(1, 2 * R + 1, dtype=np.float64)
    sinmpi2 = np.sin(np.pi * m / 2.0)
    coeff55 = 4.0 * sinmpi2 / (m * np.pi)
    kbnd = 8.0 * sinmpi2 / (m * np.pi)

    def dv(key):
        arr = np.atleast_1d(np.asarray(z["dual__" + key], np.float64))
        return arr

    l1a = dv("5.1a_w_le_Omega"); l1b = dv("5.1b_v_le_Omega")
    muO = float(dv("5.1c_Omega_le_1")[0])
    nu = float(dv("5.2_mass_eq")[0])
    lmean = float(dv("5.3_mean_ge_h1")[0]); lmom = float(dv("5.4_2ndmom_le")[0])
    l5 = dv("5.5_cos_fourier"); l6 = dv("5.6_sin_lower"); l7 = dv("5.7_sin_upper")
    l8 = dv("5.8_eps_bound"); l9 = dv("5.9_delta_bound")
    l10c = dv("5.10a_c_bound"); l10d = dv("5.10b_d_bound")
    lc1lo = float(dv("5.12_c1_ge_p1")[0]); lc1hi = float(dv("5.12_c1_le_p2")[0])
    ld1lo = float(dv("5.12_d1_ge_q1")[0]); ld1hi = float(dv("5.12_d1_le_q2")[0])
    l11 = float(dv("5.11_parseval")[0]); l13 = float(dv("5.13_A2_lower")[0])

    # ---- REPAIR: clip inequality multipliers to >= 0 (equality nu is free) ----
    def clip_nonneg(arr):
        arr = np.atleast_1d(np.asarray(arr, np.float64))
        neg = np.minimum(arr, 0.0)
        return np.maximum(arr, 0.0), float(-neg.sum())

    clip_mass = 0.0
    l1a, m0 = clip_nonneg(l1a); clip_mass += m0
    l1b, m0 = clip_nonneg(l1b); clip_mass += m0
    l5, m0 = clip_nonneg(l5); clip_mass += m0
    l6, m0 = clip_nonneg(l6); clip_mass += m0
    l7, m0 = clip_nonneg(l7); clip_mass += m0
    l8, m0 = clip_nonneg(l8); clip_mass += m0
    l9, m0 = clip_nonneg(l9); clip_mass += m0
    l10c, m0 = clip_nonneg(l10c); clip_mass += m0
    l10d, m0 = clip_nonneg(l10d); clip_mass += m0
    for nm in (muO, lmean, lmom, lc1lo, lc1hi, ld1lo, ld1hi, l11, l13):
        if nm < 0.0:
            clip_mass += -nm
    muO = max(muO, 0.0); lmean = max(lmean, 0.0); lmom = max(lmom, 0.0)
    lc1lo = max(lc1lo, 0.0); lc1hi = max(lc1hi, 0.0)
    ld1lo = max(ld1lo, 0.0); ld1hi = max(ld1hi, 0.0)
    l11 = max(l11, 0.0); l13 = max(l13, 0.0)

    j = np.arange(1, N + 1, dtype=np.float64)
    jm1 = j - 1.0
    ap2 = ap_[1]  # alpha^+_{j,2}

    # ---- L(x*, lambda) = Omega* + sum_i lambda_i g_i(x*) ----
    g51a = w - Omega_star; g51b = v - Omega_star; g51c = Omega_star - 1.0
    h52 = L * np.sum(w + v) - 1.0
    g53 = h1 - (L ** 2) * (np.sum(j * w) - np.sum(jm1 * v))
    g54 = (L ** 3) * np.sum(jm1 ** 2 * (w + v)) - (2.0 / 3.0 + h2 ** 2 / 2.0)
    lhs55 = (L / 2.0) * (am_ @ (w + v)); rhs55 = coeff55 * a - 2.0 * (a ** 2 + b ** 2)
    g55 = lhs55 - rhs55
    lhs56 = (L / 2.0) * (bm_ @ w - bp_ @ v); rhs56 = -kbnd * b; g56 = lhs56 - rhs56
    lhs57 = (L / 2.0) * (bp_ @ w - bm_ @ v); rhs57 = -kbnd * b; g57 = rhs57 - lhs57
    g58 = np.abs(eps) - eps_bound; g59 = np.abs(delta) - delta_bound
    g510c = np.abs(c) - 2.0 / np.pi; g510d = np.abs(d) - 2.0 / np.pi
    gc1lo = p1 - c[0]; gc1hi = c[0] - p2; gd1lo = q1 - d[0]; gd1hi = d[0] - q2
    g511 = np.sum(c ** 2) + np.sum(d ** 2) - 0.5
    lhs513 = (L / 2.0) * (ap2 @ (w + v)); rhs513 = -0.5 * (p2 ** 2 + max(q1 ** 2, q2 ** 2))
    g513 = rhs513 - lhs513

    Lx = (Omega_star
          + float(l1a @ g51a) + float(l1b @ g51b) + muO * g51c
          + nu * h52 + lmean * g53 + lmom * g54
          + float(l5 @ g55) + float(l6 @ g56) + float(l7 @ g57)
          + float(l8 @ g58) + float(l9 @ g59)
          + float(l10c @ g510c) + float(l10d @ g510d)
          + lc1lo * gc1lo + lc1hi * gc1hi + ld1lo * gd1lo + ld1hi * gd1hi
          + l11 * g511 + l13 * g513)

    # ---- gradients dL/dx at x* ----
    gO = 1.0 - float(l1a.sum()) - float(l1b.sum()) + muO
    gw = (l1a + nu * L - lmean * (L ** 2) * j + lmom * (L ** 3) * (jm1 ** 2)
          + (L / 2.0) * (l5 @ am_) + (L / 2.0) * (l6 @ bm_) - (L / 2.0) * (l7 @ bp_)
          - l13 * (L / 2.0) * ap2)
    gv = (l1b + nu * L + lmean * (L ** 2) * jm1 + lmom * (L ** 3) * (jm1 ** 2)
          + (L / 2.0) * (l5 @ am_) - (L / 2.0) * (l6 @ bp_) + (L / 2.0) * (l7 @ bm_)
          - l13 * (L / 2.0) * ap2)
    # dL/da_m and dL/db_m through the 5.5/5.6/5.7 terms
    dLda = l5 * (-coeff55 + 4.0 * a)
    dLdb = l6 * kbnd - l7 * kbnd
    # chain a_m, b_m -> (c, d, eps, delta)
    kk = np.arange(1, T + 1, dtype=np.float64)
    sign_k = (-1.0) ** kk
    gc = np.zeros(T); gd = np.zeros(T); geps = np.zeros(R); gdel = np.zeros(R)
    for i in range(2 * R):
        mm = i + 1
        if mm % 2 == 0:
            gc[mm // 2 - 1] += dLda[i] * 0.5
            gd[mm // 2 - 1] += dLdb[i] * 0.5
        else:
            mh = (mm + 1) // 2
            sinv = np.sin(np.pi * mm / 2.0)
            denom = (mm ** 2) - 4.0 * (kk ** 2)
            coeff_a = (2.0 * mm * sinv / np.pi) * (sign_k / denom)
            coeff_b = (4.0 / np.pi) * (kk * sign_k * sinv / denom)
            gc += dLda[i] * coeff_a; geps[mh - 1] += dLda[i]
            gd += dLdb[i] * coeff_b; gdel[mh - 1] += dLdb[i]
    gc += l11 * 2.0 * c; gd += l11 * 2.0 * d
    gc[0] += (lc1hi - lc1lo); gd[0] += (ld1hi - ld1lo)
    # subgradients of the abs terms (valid since l8,l9,l10 >= 0)
    geps += l8 * np.sign(eps); gdel += l9 * np.sign(delta)
    gc += l10c * np.sign(c); gd += l10d * np.sign(d)

    # ---- separable box-underestimator correction ----
    def boxmin(grad, xstar, lo, hi):
        grad = np.atleast_1d(np.asarray(grad, np.float64))
        xstar = np.atleast_1d(np.asarray(xstar, np.float64))
        lo = np.broadcast_to(lo, grad.shape); hi = np.broadcast_to(hi, grad.shape)
        xopt = np.where(grad > 0, lo, hi)
        return float(np.sum(grad * (xopt - xstar)))

    two_over_pi = 2.0 / np.pi
    # Real per-index simple bounds: |c_k|,|d_k| <= 2/pi (5.10) for k>=2, and the
    # TIGHT box c_1 in [p1,p2], d_1 in [q1,q2] (5.12).  Both are honest supersets
    # of the feasible region; the tight c_1/d_1 box matters a lot because the
    # optimum pins c_1 at p1, so the underestimator there collapses to ~0.
    c_lo = np.full(T, -two_over_pi); c_hi = np.full(T, two_over_pi)
    c_lo[0] = p1; c_hi[0] = p2
    d_lo = np.full(T, -two_over_pi); d_hi = np.full(T, two_over_pi)
    d_lo[0] = q1; d_hi[0] = q2
    corr_O = boxmin([gO], [Omega_star], 0.0, 1.0)
    corr_w = boxmin(gw, w, 0.0, 1.0)
    corr_v = boxmin(gv, v, 0.0, 1.0)
    corr_c = boxmin(gc, c, c_lo, c_hi)
    corr_d = boxmin(gd, d, d_lo, d_hi)
    corr_eps = boxmin(geps, eps, -eps_bound, eps_bound)
    corr_del = boxmin(gdel, delta, -delta_bound, delta_bound)
    correction = corr_O + corr_w + corr_v + corr_c + corr_d + corr_eps + corr_del

    rigorous_lb = Lx + correction

    # ---- residual dual-infeasibility (KKT stationarity gaps) ----
    resid = {
        "c_Omega": float(gO),
        "min_c_w": float(gw.min()), "min_c_v": float(gv.min()),
        "max_abs_grad_c": float(np.abs(gc).max()),
        "max_abs_grad_d": float(np.abs(gd).max()),
        "max_abs_grad_eps": float(np.abs(geps).max()),
        "max_abs_grad_delta": float(np.abs(gdel).max()),
    }

    out = {
        "source_npz": str(npz_path),
        "N": N, "R": R, "T": T,
        "h1": h1, "h2": h2, "p1": p1, "p2": p2, "q1": q1, "q2": q2,
        "Omega_star_primal": Omega_star,
        "L_at_primal_point": Lx,
        "box_underestimator_correction": correction,
        "correction_breakdown": {
            "Omega": corr_O, "w": corr_w, "v": corr_v,
            "c": corr_c, "d": corr_d, "eps": corr_eps, "delta": corr_del,
        },
        "rigorous_lower_bound_on_mu": rigorous_lb,
        "dual_repair_clip_mass": clip_mass,
        "dual_infeasibility_residuals": resid,
        "method": "weak-duality tangent-underestimator over the simple-bound box "
                  "(White Prop.9 lower bound; dual repair = clip lambda>=0).",
        "theorem_statement": (
            f"The extracted multipliers (all inequality lambda >= 0 after a clip of "
            f"total mass {clip_mass:.2e}) certify, via the tangent underestimator of "
            f"the Lagrangian over the program's simple-bound box, that the Section-5 "
            f"program optimum p* >= {rigorous_lb:.7f}. By White (arXiv:2201.05704) "
            f"Prop. 9, p* <= mu, hence mu >= {rigorous_lb:.7f}. "
            f"[Method: White. Scaling + certificate extraction + repair: ours. "
            f"Not a determination of mu.]"
        ),
    }

    if cvxpy_crosscheck:
        out["cvxpy_box_dual_objective"] = _cvxpy_box_bound(
            N, R, T, L, h1, h2, p1, p2, q1, q2, am_, ap_, bm_, bp_,
            coeff55, sinmpi2, kbnd, eps_bound, delta_bound,
            l1a, l1b, muO, nu, lmean, lmom, l5, l6, l7, l8, l9, l10c, l10d,
            lc1lo, lc1hi, ld1lo, ld1hi, l11, l13,
        )
    return out


def _cvxpy_box_bound(N, R, T, L, h1, h2, p1, p2, q1, q2, am_, ap_, bm_, bp_,
                     coeff55, sinmpi2, kbnd, eps_bound, delta_bound,
                     l1a, l1b, muO, nu, lmean, lmom, l5, l6, l7, l8, l9, l10c,
                     l10d, lc1lo, lc1hi, ld1lo, ld1hi, l11, l13):
    """Tighter (still valid) reference: inf over the box of L via a convex solve.
    NOT part of the interval-rigorous path -- a cross-check that the tangent
    underestimator is not badly loose."""
    import cvxpy as cp
    Om = cp.Variable(); w = cp.Variable(N); v = cp.Variable(N)
    c = cp.Variable(T); d = cp.Variable(T); eps = cp.Variable(R); delta = cp.Variable(R)
    kk = np.arange(1, T + 1, dtype=np.float64); sign_k = (-1.0) ** kk
    a_list = [None] * (2 * R); b_list = [None] * (2 * R)
    for i in range(2 * R):
        mm = i + 1
        if mm % 2 == 0:
            a_list[i] = 0.5 * c[mm // 2 - 1]; b_list[i] = 0.5 * d[mm // 2 - 1]
        else:
            mh = (mm + 1) // 2; sinv = np.sin(np.pi * mm / 2.0)
            denom = (mm ** 2) - 4.0 * (kk ** 2)
            coeff_a = (2.0 * mm * sinv / np.pi) * (sign_k / denom)
            coeff_b = (4.0 / np.pi) * (kk * sign_k * sinv / denom)
            k0 = (2.0 * mm * sinv / np.pi) * (1.0 / (2.0 * mm ** 2))
            a_list[i] = eps[mh - 1] + k0 + coeff_a @ c
            b_list[i] = delta[mh - 1] + coeff_b @ d
    a = cp.hstack(a_list); b = cp.hstack(b_list)
    j = np.arange(1, N + 1, dtype=np.float64); jm1 = j - 1.0; ap2 = ap_[1]
    g51a = w - Om; g51b = v - Om; g51c = Om - 1.0
    h52 = L * cp.sum(w + v) - 1.0
    g53 = h1 - (L ** 2) * (cp.sum(cp.multiply(j, w)) - cp.sum(cp.multiply(jm1, v)))
    g54 = (L ** 3) * cp.sum(cp.multiply(jm1 ** 2, w + v)) - (2.0 / 3.0 + h2 ** 2 / 2.0)
    g55 = (L / 2.0) * (am_ @ (w + v)) - (cp.multiply(coeff55, a) - 2.0 * (cp.square(a) + cp.square(b)))
    g56 = (L / 2.0) * (bm_ @ w - bp_ @ v) - (cp.multiply(-kbnd, b))
    g57 = cp.multiply(-kbnd, b) - (L / 2.0) * (bp_ @ w - bm_ @ v)
    g58 = cp.abs(eps) - eps_bound; g59 = cp.abs(delta) - delta_bound
    g510c = cp.abs(c) - 2.0 / np.pi; g510d = cp.abs(d) - 2.0 / np.pi
    g511 = cp.sum_squares(c) + cp.sum_squares(d) - 0.5
    g513 = (-0.5 * (p2 ** 2 + max(q1 ** 2, q2 ** 2))) - (L / 2.0) * (ap2 @ (w + v))
    Lag = (Om + l1a @ g51a + l1b @ g51b + muO * g51c + nu * h52 + lmean * g53
           + lmom * g54 + l5 @ g55 + l6 @ g56 + l7 @ g57 + l8 @ g58 + l9 @ g59
           + l10c @ g510c + l10d @ g510d + lc1lo * (p1 - c[0]) + lc1hi * (c[0] - p2)
           + ld1lo * (q1 - d[0]) + ld1hi * (d[0] - q2) + l11 * g511 + l13 * g513)
    box = [w >= 0, w <= 1, v >= 0, v <= 1, Om >= 0, Om <= 1,
           c >= -2 / np.pi, c <= 2 / np.pi, d >= -2 / np.pi, d <= 2 / np.pi,
           eps >= -eps_bound, eps <= eps_bound, delta >= -delta_bound, delta <= delta_bound]
    prob = cp.Problem(cp.Minimize(Lag), box)
    prob.solve(solver=cp.CLARABEL)
    return float(prob.value)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("npz", help="path to the .npz certificate dump")
    ap.add_argument("--out", default=None)
    ap.add_argument("--cvxpy-crosscheck", action="store_true")
    args = ap.parse_args()

    res = repair_certificate(args.npz, cvxpy_crosscheck=args.cvxpy_crosscheck)
    print(json.dumps(res, indent=2))
    out = args.out or args.npz.replace(".npz", "_repaired.json")
    Path(out).write_text(json.dumps(res, indent=2))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
