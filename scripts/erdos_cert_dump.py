#!/usr/bin/env python3
"""Solve White's Section-5 convex program and DUMP the full primal variable
vector as a certificate JSON (schema consumed by erdos_cert_verify.py).

The stock solver (erdos_white_dual_certificate.py) only saves a summary; this
produces the full {w,v,c,d,eps,delta,Omega} so the independent interval
harness has something to check.  Construction mirrors the solver (that is
fine -- the verifier is the independent party).
"""
from __future__ import annotations
import argparse, json, time
import numpy as np
import cvxpy as cp


def solve_dump(N, R, T, h1, h2, p1, p2, q1, q2, out, solver="CLARABEL", margin=0.0):
    # margin>0 tightens the (5.5)/(5.6)/(5.7) inequalities inward by `margin`
    # so the returned point is STRICTLY interior to the true feasible region
    # (a simple, explicitly-quantified "repair to strict feasibility").
    L = 2.0 / N
    j = np.arange(1, N + 1, dtype=np.float64)
    m = np.arange(1, 2 * R + 1, dtype=np.float64)
    arg = np.outer(m, np.pi * L * (j - 0.5) / 2.0)
    slope = np.outer(m, np.full(N, np.pi * L / 4.0))
    cosv, sinv = np.cos(arg), np.sin(arg)
    alpha_minus, alpha_plus = cosv - slope, cosv + slope
    beta_minus, beta_plus = sinv - slope, sinv + slope

    Omega = cp.Variable(name="Omega")
    w = cp.Variable(N, nonneg=True); v = cp.Variable(N, nonneg=True)
    c = cp.Variable(T); d = cp.Variable(T)
    eps = cp.Variable(R); delta = cp.Variable(R)
    cons = [w <= Omega, v <= Omega, Omega <= 1.0, L * cp.sum(w + v) == 1.0]

    jidx = np.arange(1, N + 1, dtype=np.float64); jm1 = jidx - 1.0
    cons.append((L ** 2) * (cp.sum(cp.multiply(jidx, w)) - cp.sum(cp.multiply(jm1, v))) >= h1)
    jm1sq = jm1 ** 2
    cons.append((L ** 3) * cp.sum(cp.multiply(jm1sq, w + v)) <= 2.0 / 3.0 + h2 ** 2 / 2.0)

    m_all = np.arange(1, 2 * R + 1)
    kk = np.arange(1, T + 1, dtype=np.float64); sign_k = (-1.0) ** kk
    a_list = [None] * (2 * R); b_list = [None] * (2 * R)
    for idx, mm in enumerate(m_all):
        if mm % 2 == 0:
            a_list[idx] = 0.5 * c[mm // 2 - 1]; b_list[idx] = 0.5 * d[mm // 2 - 1]
        else:
            mh = (mm + 1) // 2
            s = np.sin(np.pi * mm / 2.0)
            denom = (mm ** 2) - 4.0 * (kk ** 2)
            coeff_a = (2.0 * mm * s / np.pi) * (sign_k / denom)
            coeff_b = (4.0 / np.pi) * (kk * sign_k * s / denom)
            k0 = (2.0 * mm * s / np.pi) * (1.0 / (2.0 * mm ** 2))
            a_list[idx] = eps[mh - 1] + k0 + coeff_a @ c
            b_list[idx] = delta[mh - 1] + coeff_b @ d
    a_expr = cp.hstack(a_list); b_expr = cp.hstack(b_list)
    m_f = m_all.astype(np.float64); sinmpi2 = np.sin(np.pi * m_f / 2.0)

    cons.append((L / 2.0) * (alpha_minus @ (w + v))
                <= cp.multiply(4.0 * sinmpi2 / (m_f * np.pi), a_expr)
                - 2.0 * (cp.square(a_expr) + cp.square(b_expr)) - margin)
    cons.append((L / 2.0) * (beta_minus @ w - beta_plus @ v)
                <= cp.multiply(-8.0 * sinmpi2 / (m_f * np.pi), b_expr) - margin)
    cons.append((L / 2.0) * (beta_plus @ w - beta_minus @ v)
                >= cp.multiply(-8.0 * sinmpi2 / (m_f * np.pi), b_expr) + margin)

    m_r = np.arange(1, R + 1, dtype=np.float64)
    eps_bound = (1.0 / (4.0 - (m_r ** 2) / (T ** 2))) * (2.0 * m_r / (np.pi * np.sqrt(6.0 * T ** 3)))
    delta_bound = (1.0 / (4.0 - (m_r ** 2) / (T ** 2))) * (4.0 / (np.pi * np.sqrt(2.0 * T)))
    cons += [cp.abs(eps) <= eps_bound, cp.abs(delta) <= delta_bound,
             cp.abs(c) <= 2.0 / np.pi, cp.abs(d) <= 2.0 / np.pi,
             cp.sum_squares(c) + cp.sum_squares(d) <= 0.5,
             c[0] >= p1, c[0] <= p2, d[0] >= q1, d[0] <= q2]
    alpha_plus_2 = alpha_plus[1]
    cons.append((L / 2.0) * (alpha_plus_2 @ (w + v)) >= -0.5 * (p2 ** 2 + max(q1 ** 2, q2 ** 2)))

    prob = cp.Problem(cp.Minimize(Omega), cons)
    t0 = time.time()
    prob.solve(solver=getattr(cp, solver), verbose=False)
    dt = time.time() - t0
    print(f"status={prob.status} Omega={Omega.value} solve={dt:.1f}s")

    cert = {
        "N": N, "R": R, "T": T, "h1": h1, "h2": h2, "p1": p1, "p2": p2, "q1": q1, "q2": q2,
        "Omega": float(Omega.value), "claimed_bound": float(Omega.value),
        "status": prob.status, "solve_seconds": dt, "solver": solver, "repair_margin": margin,
        "w": [float(x) for x in np.asarray(w.value)],
        "v": [float(x) for x in np.asarray(v.value)],
        "c": [float(x) for x in np.asarray(c.value)],
        "d": [float(x) for x in np.asarray(d.value)],
        "eps": [float(x) for x in np.asarray(eps.value)],
        "delta": [float(x) for x in np.asarray(delta.value)],
    }
    json.dump(cert, open(out, "w"))
    print(f"wrote {out}  ({len(cert['w'])} w, {len(cert['c'])} c)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    for k, dv in [("N", 20000), ("R", 10), ("T", 5000)]:
        ap.add_argument(f"--{k}", type=int, default=dv)
    for k, dv in [("h1", 0.0), ("h2", 0.06), ("p1", 0.33), ("p2", 0.35), ("q1", -0.02), ("q2", 0.02)]:
        ap.add_argument(f"--{k}", type=float, default=dv)
    ap.add_argument("--solver", default="CLARABEL")
    ap.add_argument("--margin", type=float, default=0.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    solve_dump(a.N, a.R, a.T, a.h1, a.h2, a.p1, a.p2, a.q1, a.q2, a.out, a.solver, a.margin)
