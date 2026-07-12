#!/usr/bin/env python3
"""Exact rational certification of the Erdos minimum-overlap UPPER bound.

Construction: agent "Hyra" (EinsteinArena, problem erdos-min-overlap),
an n=2400 step vector f in [0,1]^n. By default the vector is read from
../certs/erdos_hyra_current.json (relative to this script); pass an explicit
path as the first argument to override.
Our contribution: exact evaluation + the continuum-equals-discrete lemma.

Verifier semantics being replicated EXACTLY (in rational arithmetic):
    f <- f * ( (n/2) / sum(f) )          # sum-normalization rescale
    score = max_k np.correlate(f, 1-f, "full")[k] * 2/n

Mapping to White (arXiv 2201.05704, eq. (1.1) / Sec. 2): with h = 2/n and
F the step function on [-1,1] taking value f_i on [-1+ih, -1+(i+1)h),
G = 1-F on [-1,1] and F = G = 0 outside [-1,1]:
    integral F = h * sum f_i = (2/n)(n/2) = 1,
    M(x) = int F(t) G(t+x) dt,  M(mh) = h * sum_i f_i g_{i+m}.
By the piecewise-linearity lemma (M is continuous piecewise linear with
breakpoints in h*Z, M(+-2)=0), sup_{x in R} M(x) = max_{|m|<=n-1} M(mh)
= the verifier score. Hence mu <= exact verifier score (Swinnerton-Dyer /
White (1.1): mu = inf_f sup_x M(x), so any admissible f upper-bounds mu).

Everything below is exact integer arithmetic: floats are converted to
rationals via their bit patterns (Fraction(float) is exact), the rescale is
done exactly, and all 2n-1 correlation values are exact integers over S^2.
"""

import json
import os
import sys
from fractions import Fraction
from math import gcd
from operator import mul

try:
    from math import sumprod  # Python >= 3.12, exact for ints
    def dot(a, b):
        return sumprod(a, b)
except ImportError:
    def dot(a, b):
        return sum(map(mul, a, b))

N = 2400
_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_JSON_PATH = os.path.join(_HERE, os.pardir, "certs", "erdos_hyra_current.json")
BOARD_VALUE = 0.3808669097979876
DIGITS = 25


def main():
    if not __debug__:
        raise SystemExit("refusing to run under -O: validation is load-bearing")
    json_path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_JSON_PATH
    vals = json.load(open(json_path))["values"]
    assert len(vals) == N, f"expected n={N}, got {len(vals)}"

    # ---- exact rationals from the float bit patterns ----
    fr = [Fraction(v) for v in vals]  # exact: float -> dyadic rational
    assert all(0 <= x <= 1 for x in fr), "raw vector leaves [0,1]"

    # common power-of-two denominator: fr_i = A_i / 2^L
    Ls = []
    for x in fr:
        q = x.denominator
        l = q.bit_length() - 1
        assert q == (1 << l), "float denominator not a power of two?!"
        Ls.append(l)
    L = max(Ls)
    A = [x.numerator << (L - l) for x, l in zip(fr, Ls)]
    S = sum(A)  # raw sum = S / 2^L  (exact)

    raw_sum = Fraction(S, 1 << L)
    target = Fraction(N, 2)  # = 1200
    dev = raw_sum - target
    print(f"n = {N}, common denominator 2^{L}")
    print(f"raw sum        = {float(raw_sum):.17g}  (exact {raw_sum.numerator}/{raw_sum.denominator})")
    print(f"raw sum - n/2  = {float(dev):.6e}  ({'>0: rescale factor < 1' if dev > 0 else '<0: rescale factor > 1' if dev < 0 else 'exact'})")

    # ---- exact verifier rescale:  f_i' = f_i * (n/2)/sum(f) = 1200*A_i / S ----
    half_n = N // 2  # 1200
    P = [half_n * a for a in A]      # numerator of f_i' over denominator S
    # range check AFTER rescale: 0 <= P_i <= S  <=>  f_i' in [0,1]
    over = [i for i, p in enumerate(P) if p > S]
    under = [i for i, p in enumerate(P) if p < 0]
    print(f"[0,1] check after exact rescale: {len(over)} entries > 1, {len(under)} entries < 0")
    if over:
        worst = max(over, key=lambda i: P[i])
        exc = Fraction(P[worst] - S, S)
        print(f"  worst excess: f[{worst}] - 1 = {float(exc):.3e}")
    assert not over and not under, "rescaled vector leaves [0,1] -- bound NOT certified"

    Q = [S - p for p in P]           # numerator of g_i = 1 - f_i' over S

    # ---- all 2n-1 correlations exactly:  c_m = (sum_i P_i * Q_{i+m}) / S^2 ----
    # m ranges over -(n-1) .. n-1;  lags |m| >= n give empty overlap (= 0).
    corr = {}
    best_m, best_t = None, -1
    for m in range(-(N - 1), N):
        lo, hi = max(0, -m), min(N, N - m)
        t = dot(P[lo:hi], Q[lo + m:hi + m])
        corr[m] = t
        if t > best_t:
            best_t, best_m = t, m
    assert len(corr) == 2 * N - 1
    assert all(t >= 0 for t in corr.values())

    # ---- exact score = (2/n) * max_m c_m = max_t / (1200 * S^2) ----
    score = Fraction(2 * best_t, N * S * S)  # Fraction auto-reduces
    print(f"\nargmax lag m* = {best_m}  (x* = m*·h = {2*best_m}/{N} = {2*best_m/N:+.6f})")
    print(f"exact score  = {score.numerator}\n             / {score.denominator}")

    # 25-digit decimal, truncated and rounded UP
    p, q = score.numerator, score.denominator
    scaled = p * 10**DIGITS
    trunc = scaled // q
    ceil_ = -((-scaled) // q)
    def fmt(ip):
        s = str(ip).rjust(DIGITS + 1, "0")
        return s[:-DIGITS] + "." + s[-DIGITS:]
    print(f"decimal (truncated, {DIGITS} digits) = {fmt(trunc)}")
    print(f"decimal (ROUNDED UP, {DIGITS} digits) = {fmt(ceil_)}")
    print(f"RIGOROUS upper bound:  mu <= {fmt(ceil_)}  =  {ceil_}/10^{DIGITS}  (exact rational)")
    print(f"float(score) = {float(score):.17g}")
    print(f"board value  = {BOARD_VALUE:.17g}   match to float precision: {float(score) == BOARD_VALUE}")

    # ---- cross-check against numpy's own computation & indexing convention ----
    try:
        import numpy as np
    except ImportError:
        print("numpy unavailable -- skipping float cross-check")
        return
    f = np.array(vals, dtype=np.float64)
    f = f * ((N / 2) / f.sum())
    g = 1.0 - f
    c = np.correlate(f, g, "full") * 2.0 / N
    print(f"\nnumpy float verifier score = {c.max():.17g}")
    # pin the indexing convention: np.correlate(f, g, 'full')[k]  <->  which lag m?
    # test both orientations element-wise against the exact values
    exact_arr = [2 * corr[m] / (N * S * S) for m in range(-(N - 1), N)]
    exact_f = np.array([float(Fraction(2 * corr[m], N * S * S)) for m in range(-(N - 1), N)])
    fwd = np.max(np.abs(c - exact_f))
    rev = np.max(np.abs(c - exact_f[::-1]))
    print(f"element-wise |numpy - exact| : same-order max {fwd:.3e}, reversed max {rev:.3e}")
    conv = "c[k] = sum_i f_i g_(i+m), m = k-(n-1)" if fwd < rev else \
           "c[k] = sum_i f_i g_(i+m), m = (n-1)-k"
    print(f"numpy 'full' convention pinned: {conv}; max discrepancy {min(fwd, rev):.3e} (float rounding only)")
    assert min(fwd, rev) < 1e-12
    # 'full' mode length check: covers ALL lags with nonempty overlap
    assert len(c) == 2 * N - 1


if __name__ == "__main__":
    main()
