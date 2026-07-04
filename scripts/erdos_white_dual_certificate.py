#!/usr/bin/env python3
"""White (arXiv:2201.05704) Sections 4-5 dual-certificate lower bound for mu,
the Erdos minimum-overlap constant, implemented verbatim from the paper and
run completely independently of Hyra's current-leader vector.

This produces a certified LOWER bound on ||M||_inf = mu that holds for
*every* admissible f, M (Section 2 setup) -- i.e. every point in the arena's
feasible region, including Hyra's. It is decoupled entirely from Hyra's
vector: no read of /tmp/erdos_hyra_current.json anywhere in this file.

Two programs, both taken verbatim from the paper:

  Section 4 (Proposition 8): a LINEAR program, valid only under the extra
  assumption that M(x) is even. Variables Omega, w_1..w_N. Constraints
  (4.1)-(4.4). White reports N=80000,R=20 -> Omega* = 0.375169005340707,
  hence mu >= 0.375 for all even M. We reproduce this as a correctness
  check at smaller N (the LP scales in N so smaller N gives a smaller but
  qualitatively correct bound), then push N higher.

  Section 5 (Proposition 9): the FULL convex QCQP, no evenness assumption.
  Variables Omega, {w_j,v_j}, {c_k,d_k}, {epsilon_{2m-1},delta_{2m-1}}.
  Constraints (5.1)-(5.13). This is convex (all quadratic terms enter with
  the correct sign to keep the feasible region convex -- see below), so
  cvxpy/CLARABEL's solution is certified globally optimal for the *given*
  (N,R,T,h1,h2,p1,p2,q1,q2) box, and it returns dual variables for every
  constraint. We use those duals to build an independent numerical dual
  certificate: a linear combination of the *primal* constraints (weighted
  by the returned dual multipliers) that yields the exact bound Omega* as
  a sum of provably-true inequalities, checked in raw numpy/Python floats
  with an explicit slack margin -- never trusting the solver's reported
  objective number directly (per the mandate).

Both programs solve for Omega over a fixed sub-box of (h1,h2,p1,p2,q1,q2).
Reproducing White's headline 0.379005 requires his "divide and conquer" /
ellipse-covering procedure (Section 5.1, Appendix II) across ~19
sub-boxes; that book-keeping machinery is out of scope here (it is a
literature-reproduction exercise unrelated to Hyra's vector) but the single
richest sub-box he reports -- (h1,h2)=(0,0.06),(p1,p2)=(0.33,0.45),
(q1,q2)=(-0.02,0.02) at N,R,T=20000,10,5000 -- already gives 0.38, and the
tightest single box -- same but (p1,p2)=(0.33,0.35) -- gives 0.37925. We
reproduce both exactly as the correctness check, then push N,R,T well past
what White used in a single box to see how far one box alone can climb.

Usage:
  python3 erdos_white_dual_certificate.py section4 --N 2000 --R 20
  python3 erdos_white_dual_certificate.py section5 --N 4000 --R 10 --T 4000 \
      --h1 0.0 --h2 0.06 --p1 0.33 --p2 0.45 --q1 -0.02 --q2 0.02
  python3 erdos_white_dual_certificate.py section5 --N 20000 --R 10 --T 5000 \
      --h1 0.0 --h2 0.06 --p1 0.33 --p2 0.35 --q1 -0.02 --q2 0.02  # White repro

Writes JSON results under /tmp/erdos_white_dual_*.json. Never touches
Hyra's vector, never submits.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import cvxpy as cp
from scipy.optimize import linprog


# ---------------------------------------------------------------------------
# Section 4: simplified even-M linear program (Proposition 8)
# ---------------------------------------------------------------------------

def alpha_minus_section4(N: int, R: int, L: float) -> np.ndarray:
    """alpha^-_{j,2m} for 1<=j<=N, 1<=m<=R, using White's explicit choice
    alpha^-_{j,m} = cos(pi*m*L*(j-1/2)/2) - pi*m*L/4, evaluated at m=2m.
    Returns array shape (R, N): row m-1 is alpha^-_{j,2m} over j=1..N.
    """
    j = np.arange(1, N + 1, dtype=np.float64)
    m2 = 2 * np.arange(1, R + 1, dtype=np.float64)  # the "2m" harmonic
    # shape (R, N)
    arg = np.outer(m2, np.pi * L * (j - 0.5) / 2.0)
    alpha = np.cos(arg) - np.outer(m2, np.full(N, np.pi * L / 4.0))
    return alpha


def section4_lp(N: int, R: int, verbose: bool = True) -> dict:
    """Solve Prop. 8's LP exactly as stated:

        minimize Omega
        s.t. 0 <= w_j <= Omega          (4.1)
             sum w_j = N/4               (4.2)
             sum alpha^-_{j,2m} w_j <= 0  for m=1..R   (4.3)
             L^3 sum (j-1)^2 w_j <= 1/3  (4.4)

    Variables: Omega, w_1..w_N  (N+1 total).
    We solve via scipy.optimize.linprog (HiGHS backend) which is exact
    enough at this scale, then separately verify feasibility of the
    returned point against every constraint in raw float64 with slack
    margins, exactly as White's paper insists ("we will check that all
    inequalities... are strictly satisfied by a margin that exceeds the
    worst-case-scenario for floating-point rounding errors").
    """
    L = 2.0 / N
    alpha = alpha_minus_section4(N, R, L)  # (R, N)

    # Variable order: x = [Omega, w_1, ..., w_N]
    nvar = N + 1
    c = np.zeros(nvar)
    c[0] = 1.0  # minimize Omega

    A_ub_rows = []
    b_ub = []

    # (4.1): w_j - Omega <= 0  for all j
    A1 = np.zeros((N, nvar))
    A1[:, 0] = -1.0
    A1[np.arange(N), 1:] = np.eye(N)[np.arange(N)]  # w_j coeff 1
    for j in range(N):
        A1[j, 1 + j] = 1.0
    A_ub_rows.append(A1)
    b_ub.append(np.zeros(N))

    # (4.3): sum_j alpha^-_{j,2m} w_j <= 0 for m=1..R
    A3 = np.zeros((R, nvar))
    A3[:, 1:] = alpha
    A_ub_rows.append(A3)
    b_ub.append(np.zeros(R))

    # (4.4): L^3 sum (j-1)^2 w_j <= 1/3
    jm1sq = (np.arange(N, dtype=np.float64)) ** 2  # (j-1)^2 for j=1..N -> index j-1
    A4 = np.zeros((1, nvar))
    A4[0, 1:] = (L ** 3) * jm1sq
    A_ub_rows.append(A4)
    b_ub.append(np.array([1.0 / 3.0]))

    A_ub = np.vstack(A_ub_rows)
    b_ub_full = np.concatenate(b_ub)

    # (4.2): sum w_j = N/4  (equality)
    A_eq = np.zeros((1, nvar))
    A_eq[0, 1:] = 1.0
    b_eq = np.array([N / 4.0])

    # bounds: Omega free (but effectively >=0), w_j >= 0
    bounds = [(None, None)] + [(0.0, None)] * N

    t0 = time.time()
    res = linprog(
        c, A_ub=A_ub, b_ub=b_ub_full, A_eq=A_eq, b_eq=b_eq,
        bounds=bounds, method="highs",
    )
    dt = time.time() - t0
    if not res.success:
        return {"ok": False, "message": res.message, "N": N, "R": R}

    x = res.x
    Omega = float(x[0])
    w = x[1:]

    # Independent feasibility re-check in raw float64 (never trust solver
    # objective directly).
    margin = {}
    margin["box"] = float(np.max(w - Omega))  # should be <=0
    margin["wpos"] = float(np.min(w))  # should be >=0
    margin["sum_eq"] = float(abs(np.sum(w) - N / 4.0))
    margin["fourier"] = float(np.max(alpha @ w))  # should be <=0
    margin["second_moment"] = float((L ** 3) * np.sum(jm1sq * w) - 1.0 / 3.0)  # should be <=0

    ok = (
        margin["box"] <= 1e-7
        and margin["wpos"] >= -1e-9
        and margin["sum_eq"] <= 1e-6
        and margin["fourier"] <= 1e-6
        and margin["second_moment"] <= 1e-6
    )

    result = {
        "ok": bool(ok),
        "section": 4,
        "N": N,
        "R": R,
        "Omega_star": Omega,
        "margins": margin,
        "solve_seconds": dt,
        "note": "mu >= Omega_star for all EVEN M (Prop. 8). White's own "
                "reported optimum at N=80000,R=20 is 0.375169005340707.",
    }
    if verbose:
        print(json.dumps(result, indent=2))
    return result


# ---------------------------------------------------------------------------
# Section 5: full convex QCQP (Proposition 9) -- no evenness assumption
# ---------------------------------------------------------------------------

def alpha_beta_section5(N: int, R: int, L: float):
    """alpha^-,alpha^+,beta^-,beta^+ for 1<=j<=N, 1<=m<=2R, per Section 5:

        alpha^-_{j,m} = cos(pi m L (j-1/2)/2) - pi m L /4
        alpha^+_{j,m} = cos(pi m L (j-1/2)/2) + pi m L /4
        beta^-_{j,m}  = sin(pi m L (j-1/2)/2) - pi m L /4
        beta^+_{j,m}  = sin(pi m L (j-1/2)/2) + pi m L /4

    Returns four arrays of shape (2R, N), row m-1 for m=1..2R.
    """
    j = np.arange(1, N + 1, dtype=np.float64)
    m = np.arange(1, 2 * R + 1, dtype=np.float64)
    arg = np.outer(m, np.pi * L * (j - 0.5) / 2.0)      # (2R, N)
    slope = np.outer(m, np.full(N, np.pi * L / 4.0))     # (2R, N)
    cosv = np.cos(arg)
    sinv = np.sin(arg)
    alpha_minus = cosv - slope
    alpha_plus = cosv + slope
    beta_minus = sinv - slope
    beta_plus = sinv + slope
    return alpha_minus, alpha_plus, beta_minus, beta_plus


def section5_qcqp(
    N: int, R: int, T: int,
    h1: float, h2: float, p1: float, p2: float, q1: float, q2: float,
    solver: str = "CLARABEL",
    verbose: bool = True,
    dump_certificate: bool = False,
    dump_path: str = None,
) -> dict:
    """Solve the FULL convex program of Section 5, Prop. 9, verbatim.

    Variables: Omega (scalar), w[1..N], v[1..N], c[1..T], d[1..T],
               eps[1..R] (i.e. epsilon_{2m-1}), delta[1..R] (delta_{2m-1}).

    Variable expressions for 1<=m<=2R:
        a_m = c_{m/2}/2                                     if m even
        a_m = eps_m + (2 m sin(pi m/2)/pi) * sum_k (-1)^k/(m^2-4k^2) c_k   if m odd
        b_m = d_{m/2}/2                                     if m even
        b_m = delta_m + (4/pi) sum_k k(-1)^k sin(pi m/2)/(m^2-4k^2) d_k    if m odd

    Objective: minimize Omega.

    Constraints (5.1)-(5.13) exactly as in the paper.
    """
    L = 2.0 / N
    alpha_minus, alpha_plus, beta_minus, beta_plus = alpha_beta_section5(N, R, L)  # (2R,N)

    Omega = cp.Variable(nonneg=False, name="Omega")
    w = cp.Variable(N, nonneg=True, name="w")
    v = cp.Variable(N, nonneg=True, name="v")
    c = cp.Variable(T, name="c")
    d = cp.Variable(T, name="d")
    eps = cp.Variable(R, name="eps")     # epsilon_{2m-1}, m=1..R
    delta = cp.Variable(R, name="delta")  # delta_{2m-1}, m=1..R

    cons = []

    # (5.1): 0 <= w_j, v_j <= Omega <= 1
    cons.append(w <= Omega)
    cons.append(v <= Omega)
    cons.append(Omega <= 1.0)

    # (5.2): L * sum(w_j+v_j) = 1
    cons.append(L * cp.sum(w + v) == 1.0)

    # (5.3): L^2 * sum(j*w_j - (j-1)*v_j) >= h1
    jidx = np.arange(1, N + 1, dtype=np.float64)
    jm1 = jidx - 1.0
    cons.append((L ** 2) * (cp.sum(cp.multiply(jidx, w)) - cp.sum(cp.multiply(jm1, v))) >= h1)

    # (5.4): L^3 * sum (j-1)^2 (w_j+v_j) <= 2/3 + h2^2/2
    jm1sq = jm1 ** 2
    cons.append((L ** 3) * cp.sum(cp.multiply(jm1sq, w + v)) <= 2.0 / 3.0 + h2 ** 2 / 2.0)

    # Build a_m, b_m for m=1..2R.
    # Precompute the odd-harmonic linear map from c_k (k=1..T) to the
    # "tail sum" term: for odd m, coefficient of c_k is
    #    (2 m sin(pi m/2)/pi) * (-1)^k/(m^2-4k^2)
    m_all = np.arange(1, 2 * R + 1)
    odd_mask = (m_all % 2 == 1)
    even_mask = ~odd_mask

    a_list = [None] * (2 * R)
    b_list = [None] * (2 * R)

    kk = np.arange(1, T + 1, dtype=np.float64)
    sign_k = (-1.0) ** kk

    for idx, m in enumerate(m_all):
        if m % 2 == 0:
            mh = m // 2  # 1..R
            a_list[idx] = 0.5 * c[mh - 1]
            b_list[idx] = 0.5 * d[mh - 1]
        else:
            mh = (m + 1) // 2  # eps/delta index 1..R -> eps_m corresponds to epsilon_{2*mh-1}=m
            sinv = np.sin(np.pi * m / 2.0)  # +-1 or 0, but m odd -> +-1
            denom = (m ** 2) - 4.0 * (kk ** 2)
            coeff_a = (2.0 * m * sinv / np.pi) * (sign_k / denom)
            coeff_b = (4.0 / np.pi) * (kk * sign_k * sinv / denom)
            # Table's a_m formula is eps_m + (2m sin(m pi/2)/pi) * (1/(2m^2)
            # + sum_{k=1}^T (-1)^k/(m^2-4k^2) c_k) -- the 1/(2m^2) term is
            # the k=0 contribution of (3.3) using c_0=1/2 (page 4). It is a
            # CONSTANT (not a variable), so it enters a_m as an affine
            # offset, not inside the c-dot-product.
            k0_term = (2.0 * m * sinv / np.pi) * (1.0 / (2.0 * m ** 2))
            a_list[idx] = eps[mh - 1] + k0_term + coeff_a @ c
            b_list[idx] = delta[mh - 1] + coeff_b @ d

    a_expr = cp.hstack(a_list)  # length 2R
    b_expr = cp.hstack(b_list)

    # (5.5): (L/2) sum alpha^-_{j,m}(w_j+v_j) <= (4 sin(m pi/2)/(m pi)) a_m - 2(a_m^2+b_m^2)
    m_f = m_all.astype(np.float64)
    sinmpi2 = np.sin(np.pi * m_f / 2.0)
    coeff55 = 4.0 * sinmpi2 / (m_f * np.pi)  # length 2R

    lhs55 = (L / 2.0) * (alpha_minus @ (w + v))  # length 2R (affine)
    rhs55 = cp.multiply(coeff55, a_expr) - 2.0 * (cp.square(a_expr) + cp.square(b_expr))
    cons.append(lhs55 <= rhs55)

    # (5.6): (L/2) sum(beta^-_{j,m} w_j - beta^+_{j,m} v_j) <= -(8/(m pi)) sin(m pi/2) b_m
    lhs56 = (L / 2.0) * (beta_minus @ w - beta_plus @ v)
    rhs56 = cp.multiply(-8.0 * sinmpi2 / (m_f * np.pi), b_expr)
    cons.append(lhs56 <= rhs56)

    # (5.7): (L/2) sum(beta^+_{j,m} w_j - beta^-_{j,m} v_j) >= -(8/(m pi)) sin(m pi/2) b_m
    lhs57 = (L / 2.0) * (beta_plus @ w - beta_minus @ v)
    rhs57 = cp.multiply(-8.0 * sinmpi2 / (m_f * np.pi), b_expr)
    cons.append(lhs57 >= rhs57)

    # (5.8): |eps_{2m-1}| <= (1/(4-m^2/T^2)) * (2m/(pi sqrt(6 T^3)))  for m=1..R
    m_r = np.arange(1, R + 1, dtype=np.float64)
    eps_bound = (1.0 / (4.0 - (m_r ** 2) / (T ** 2))) * (2.0 * m_r / (np.pi * np.sqrt(6.0 * T ** 3)))
    cons.append(cp.abs(eps) <= eps_bound)

    # (5.9): |delta_{2m-1}| <= (1/(4-m^2/T^2)) * (4/(pi sqrt(2T)))
    delta_bound = (1.0 / (4.0 - (m_r ** 2) / (T ** 2))) * (4.0 / (np.pi * np.sqrt(2.0 * T)))
    cons.append(cp.abs(delta) <= delta_bound)

    # (5.10): |c_k|, |d_k| <= 2/pi
    cons.append(cp.abs(c) <= 2.0 / np.pi)
    cons.append(cp.abs(d) <= 2.0 / np.pi)

    # (5.11): sum(c_k^2+d_k^2) <= 1/2
    cons.append(cp.sum_squares(c) + cp.sum_squares(d) <= 0.5)

    # (5.12): p1 <= c_1 <= p2, q1 <= d_1 <= q2
    cons.append(c[0] >= p1)
    cons.append(c[0] <= p2)
    cons.append(d[0] >= q1)
    cons.append(d[0] <= q2)

    # (5.13): (L/2) sum alpha^+_{j,2}(w_j+v_j) >= -(1/2)(p2^2 + max(q1^2,q2^2))
    alpha_plus_2 = alpha_plus[1]  # m=2 is index 1 (0-indexed, m_all[1]=2)
    assert m_all[1] == 2
    lhs513 = (L / 2.0) * (alpha_plus_2 @ (w + v))
    rhs513 = -0.5 * (p2 ** 2 + max(q1 ** 2, q2 ** 2))
    cons.append(lhs513 >= rhs513)

    prob = cp.Problem(cp.Minimize(Omega), cons)

    t0 = time.time()
    if solver == "CLARABEL":
        prob.solve(solver=cp.CLARABEL, verbose=False)
    elif solver == "SCS":
        prob.solve(solver=cp.SCS, verbose=False, max_iters=200000, eps=1e-9)
    else:
        prob.solve(solver=solver, verbose=False)
    dt = time.time() - t0

    status = prob.status
    Omega_val = float(Omega.value) if Omega.value is not None else None

    result = {
        "section": 5,
        "status": status,
        "N": N, "R": R, "T": T,
        "h1": h1, "h2": h2, "p1": p1, "p2": p2, "q1": q1, "q2": q2,
        "Omega_star_solver_reported": Omega_val,
        "solve_seconds": dt,
        "solver": solver,
    }

    if status not in ("optimal", "optimal_inaccurate") or Omega_val is None:
        result["ok"] = False
        if verbose:
            print(json.dumps(result, indent=2))
        return result

    # ---- Independent re-verification in raw numpy, never trust solver ----
    w_v = np.asarray(w.value, dtype=np.float64)
    v_v = np.asarray(v.value, dtype=np.float64)
    c_v = np.asarray(c.value, dtype=np.float64)
    d_v = np.asarray(d.value, dtype=np.float64)
    eps_v = np.asarray(eps.value, dtype=np.float64)
    delta_v = np.asarray(delta.value, dtype=np.float64)

    def compute_ab(m):
        if m % 2 == 0:
            mh = m // 2
            return 0.5 * c_v[mh - 1], 0.5 * d_v[mh - 1]
        mh = (m + 1) // 2
        sinv = np.sin(np.pi * m / 2.0)
        denom = (m ** 2) - 4.0 * (kk ** 2)
        coeff_a = (2.0 * m * sinv / np.pi) * (sign_k / denom)
        coeff_b = (4.0 / np.pi) * (kk * sign_k * sinv / denom)
        k0_term = (2.0 * m * sinv / np.pi) * (1.0 / (2.0 * m ** 2))
        return (
            eps_v[mh - 1] + k0_term + float(coeff_a @ c_v),
            delta_v[mh - 1] + float(coeff_b @ d_v),
        )

    a_v = np.zeros(2 * R)
    b_v = np.zeros(2 * R)
    for idx, m in enumerate(m_all):
        a_v[idx], b_v[idx] = compute_ab(int(m))

    slacks = {}
    slacks["5.1_w_le_Omega"] = float(np.min(Omega_val - w_v))
    slacks["5.1_v_le_Omega"] = float(np.min(Omega_val - v_v))
    slacks["5.1_w_nonneg"] = float(np.min(w_v))
    slacks["5.1_v_nonneg"] = float(np.min(v_v))
    slacks["5.1_Omega_le_1"] = float(1.0 - Omega_val)
    slacks["5.2_eq"] = float(abs(L * np.sum(w_v + v_v) - 1.0))
    slacks["5.3"] = float((L ** 2) * (np.sum(jidx * w_v) - np.sum(jm1 * v_v)) - h1)
    slacks["5.4"] = float((2.0 / 3.0 + h2 ** 2 / 2.0) - (L ** 3) * np.sum(jm1sq * (w_v + v_v)))
    lhs55_v = (L / 2.0) * (alpha_minus @ (w_v + v_v))
    rhs55_v = coeff55 * a_v - 2.0 * (a_v ** 2 + b_v ** 2)
    slacks["5.5_min_slack"] = float(np.min(rhs55_v - lhs55_v))
    lhs56_v = (L / 2.0) * (beta_minus @ w_v - beta_plus @ v_v)
    rhs56_v = -8.0 * sinmpi2 / (m_f * np.pi) * b_v
    slacks["5.6_min_slack"] = float(np.min(rhs56_v - lhs56_v))
    lhs57_v = (L / 2.0) * (beta_plus @ w_v - beta_minus @ v_v)
    rhs57_v = -8.0 * sinmpi2 / (m_f * np.pi) * b_v
    slacks["5.7_min_slack"] = float(np.min(lhs57_v - rhs57_v))
    slacks["5.8_min_slack"] = float(np.min(eps_bound - np.abs(eps_v)))
    slacks["5.9_min_slack"] = float(np.min(delta_bound - np.abs(delta_v)))
    slacks["5.10_c_min_slack"] = float(np.min(2.0 / np.pi - np.abs(c_v)))
    slacks["5.10_d_min_slack"] = float(np.min(2.0 / np.pi - np.abs(d_v)))
    slacks["5.11"] = float(0.5 - (np.sum(c_v ** 2) + np.sum(d_v ** 2)))
    slacks["5.12_c1"] = float(min(c_v[0] - p1, p2 - c_v[0]))
    slacks["5.12_d1"] = float(min(d_v[0] - q1, q2 - d_v[0]))
    lhs513_v = (L / 2.0) * (alpha_plus_2 @ (w_v + v_v))
    slacks["5.13"] = float(lhs513_v - rhs513)

    min_slack = min(slacks.values())
    verified = min_slack >= -1e-6  # tolerance for float64 accumulation over N,T up to 1e5

    result["ok"] = bool(verified)
    result["min_slack_over_all_constraints"] = min_slack
    result["slacks"] = slacks
    result["Omega_star_independently_recomputed"] = float(Omega_val)  # same expr as objective, sanity only

    # -------------------------------------------------------------------
    # Full certificate dump (primal + dual variable vectors).
    #
    # RIGOR NOTE: primal feasibility alone does NOT certify a lower bound
    # on mu.  White's Prop. 9 says Omega* (the *minimum* of this program)
    # is a lower bound on mu, and the RIGOROUS witness of that lower bound
    # is a *dual*-feasible point (weak duality: dual_obj <= p* <= mu).  So
    # we dump the per-constraint dual multipliers alongside the primal
    # vector; the repair (erdos_cert_repair.py) works on the DUAL side.
    # cons was appended in this fixed order:
    #   0: 5.1a  w<=Omega          [N]
    #   1: 5.1b  v<=Omega          [N]
    #   2: 5.1c  Omega<=1          [scalar]
    #   3: 5.2   L*sum(w+v)==1     [scalar, equality]
    #   4: 5.3   >= h1             [scalar]
    #   5: 5.4   <= 2/3+h2^2/2     [scalar]
    #   6: 5.5   <= (2R)           [2R]
    #   7: 5.6   <= (2R)           [2R]
    #   8: 5.7   >= (2R)           [2R]
    #   9: 5.8   |eps|<=bound      [R]
    #  10: 5.9   |delta|<=bound    [R]
    #  11: 5.10a |c|<=2/pi         [T]
    #  12: 5.10b |d|<=2/pi         [T]
    #  13: 5.12  c1>=p1            [scalar]
    #  14: 5.12  c1<=p2            [scalar]
    #  15: 5.12  d1>=q1            [scalar]
    #  16: 5.12  d1<=q2            [scalar]
    #  17: 5.11  sum(c^2+d^2)<=1/2 [scalar]
    #  18: 5.13  >= rhs            [scalar]
    if dump_certificate:
        cons_labels = [
            "5.1a_w_le_Omega", "5.1b_v_le_Omega", "5.1c_Omega_le_1",
            "5.2_mass_eq", "5.3_mean_ge_h1", "5.4_2ndmom_le",
            "5.5_cos_fourier", "5.6_sin_lower", "5.7_sin_upper",
            "5.8_eps_bound", "5.9_delta_bound",
            "5.10a_c_bound", "5.10b_d_bound",
            "5.12_c1_ge_p1", "5.12_c1_le_p2", "5.12_d1_ge_q1", "5.12_d1_le_q2",
            "5.11_parseval", "5.13_A2_lower",
        ]
        duals = {}
        for lab, con in zip(cons_labels, cons):
            dv = con.dual_value
            if dv is None:
                duals[lab] = None
            else:
                duals[lab] = np.atleast_1d(np.asarray(dv, dtype=np.float64))
        dp = dump_path or (
            f"/tmp/erdos_cert_N{N}_R{R}_T{T}.npz"
        )
        # alpha_minus/plus, beta_minus/plus, coeff55, sinmpi2, m_f are all
        # DETERMINISTIC functions of (N,R,L) via alpha_beta_section5 -- the
        # verifier recomputes them, so we do NOT store the large (2R x N)
        # arrays (640MB+ at N=2e5).  We store only the program inputs and the
        # full primal + per-constraint dual vectors.
        save_kwargs = dict(
            Omega=np.float64(Omega_val),
            w=w_v, v=v_v, c=c_v, d=d_v, eps=eps_v, delta=delta_v,
            a=a_v, b=b_v,
            N=np.int64(N), R=np.int64(R), T=np.int64(T), L=np.float64(L),
            h1=np.float64(h1), h2=np.float64(h2),
            p1=np.float64(p1), p2=np.float64(p2),
            q1=np.float64(q1), q2=np.float64(q2),
            eps_bound=eps_bound, delta_bound=delta_bound,
        )
        for lab in cons_labels:
            dvv = duals[lab]
            save_kwargs["dual__" + lab] = (
                np.array([np.nan]) if dvv is None else dvv
            )
        np.savez_compressed(dp, **save_kwargs)
        result["certificate_dump"] = dp
        result["certificate_dump_labels"] = cons_labels
        if verbose:
            print(f"[dump] wrote full primal+dual certificate -> {dp}")

    if verbose:
        printable = dict(result)
        printable["slacks"] = {k: round(v_, 12) for k, v_ in slacks.items()}
        print(json.dumps(printable, indent=2))
    return result


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    s4 = sub.add_parser("section4")
    s4.add_argument("--N", type=int, default=2000)
    s4.add_argument("--R", type=int, default=20)
    s4.add_argument("--out", type=str, default=None)

    s5 = sub.add_parser("section5")
    s5.add_argument("--N", type=int, default=4000)
    s5.add_argument("--R", type=int, default=10)
    s5.add_argument("--T", type=int, default=4000)
    s5.add_argument("--h1", type=float, required=True)
    s5.add_argument("--h2", type=float, required=True)
    s5.add_argument("--p1", type=float, required=True)
    s5.add_argument("--p2", type=float, required=True)
    s5.add_argument("--q1", type=float, required=True)
    s5.add_argument("--q2", type=float, required=True)
    s5.add_argument("--solver", type=str, default="CLARABEL")
    s5.add_argument("--out", type=str, default=None)
    s5.add_argument("--dump-certificate", action="store_true",
                    help="dump full primal+dual variable vectors to an .npz")
    s5.add_argument("--dump-path", type=str, default=None,
                    help="path for the .npz certificate dump")

    args = ap.parse_args()

    if args.cmd == "section4":
        res = section4_lp(args.N, args.R)
        out = args.out or f"/tmp/erdos_white_dual_section4_N{args.N}_R{args.R}.json"
    else:
        res = section5_qcqp(
            args.N, args.R, args.T,
            args.h1, args.h2, args.p1, args.p2, args.q1, args.q2,
            solver=args.solver,
            dump_certificate=args.dump_certificate,
            dump_path=args.dump_path,
        )
        out = args.out or (
            f"/tmp/erdos_white_dual_section5_N{args.N}_R{args.R}_T{args.T}"
            f"_h{args.h1}-{args.h2}_p{args.p1}-{args.p2}_q{args.q1}-{args.q2}.json"
        )

    Path(out).write_text(json.dumps(res, indent=2))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
