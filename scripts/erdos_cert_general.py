#!/usr/bin/env python3
"""General exact-rational certification of an Erdos min-overlap construction.

Replicates the board verifier EXACTLY in rational arithmetic:
  f <- f * ((n/2)/sum(f))            # sum-normalization
  score = max_k correlate(f, 1-f, 'full')[k] * (2/n)
and (piecewise-linear lemma: M continuous PL, breakpoints in hZ, M(+-2)=0)
=> mu <= exact score for any admissible f in [0,1]^n with the rescale in [0,1].

Usage: erdos_cert_general.py <values.json>   (json has {"values":[...]})
Prints the exact rational score, its 25-digit round-UP, and float cross-check.
"""
import json, sys
from fractions import Fraction
try:
    from math import sumprod
    def dot(a, b): return sumprod(a, b)
except ImportError:
    from operator import mul
    def dot(a, b): return sum(map(mul, a, b))

DIGITS = 25

def certify(vals):
    n = len(vals)
    fr = [Fraction(v) for v in vals]
    assert all(0 <= x <= 1 for x in fr), "raw vector leaves [0,1]"
    # common power-of-two denominator
    Ls = []
    for x in fr:
        q = x.denominator
        l = q.bit_length() - 1
        assert q == (1 << l), "float denom not power of two"
        Ls.append(l)
    L = max(Ls)
    A = [x.numerator << (L - l) for x, l in zip(fr, Ls)]
    S = sum(A)
    half_n = Fraction(n, 2)
    # exact rescale f_i' = f_i*(n/2)/sum = (n/2)*A_i / S  ; keep as numerator over S (times 2 for odd n)
    # use rationals directly to stay general for odd n
    P = [Fraction(n, 2) * a / S for a in A]   # f_i' exact in [0,1]
    over = [i for i, p in enumerate(P) if p > 1]
    under = [i for i, p in enumerate(P) if p < 0]
    assert not over and not under, f"rescaled leaves [0,1]: {len(over)}>1 {len(under)}<0"
    Q = [1 - p for p in P]                      # g_i = 1 - f_i'
    best_t, best_m = Fraction(-1), None
    for m in range(-(n - 1), n):
        lo, hi = max(0, -m), min(n, n - m)
        t = sum(P[i] * Q[i + m] for i in range(lo, hi))
        if t > best_t:
            best_t, best_m = t, m
    score = Fraction(2, n) * best_t
    return score, best_m, n

def _fmt(ip):
    s = str(ip).rjust(DIGITS + 1, "0")
    return s[:-DIGITS] + "." + s[-DIGITS:]

def fmt_down(score):
    p, q = score.numerator, score.denominator
    return _fmt((p * 10**DIGITS) // q)

def fmt_up(score):
    p, q = score.numerator, score.denominator
    return _fmt(-((-(p * 10**DIGITS)) // q))

def main():
    vals = json.load(open(sys.argv[1]))["values"]
    n = len(vals)
    fr = [Fraction(v) for v in vals]
    raw = sum(fr)                         # exact rational
    # exact post-rescale sum is (n/2)/sum(f) * sum(f) = n/2 identically
    print(f"n={n}  raw sum(f_i) = {float(raw):.12g}"
          f"{' (= n/2 exactly, admissible as given)' if raw == Fraction(n, 2) else f' (rescaled to n/2={Fraction(n,2)} exactly before scoring)'}")
    score, m, n = certify(vals)
    print(f"argmax lag m*={m}  (x* = 2*m*/n = {Fraction(2*m, n)} = {float(Fraction(2*m, n)):+.6f})")
    print(f"exact score Q = {score.numerator}/{score.denominator}")
    print(f"float(Q)      = {float(score):.17g}")
    print(f"decimal enclosure ({DIGITS} digits):")
    print(f"    {fmt_down(score)}  <=  Q  <  {fmt_up(score)}")
    print(f"RIGOROUS upper bound:  mu <= Q < {fmt_up(score)}")
    # float cross-check
    try:
        import numpy as np
        f = np.array(vals, float); f = f * ((n/2)/f.sum()); g = 1-f
        C = np.correlate(f, g, "full").max() * 2.0/n
        print(f"numpy float score = {C:.17g}  match={abs(C-float(score))<1e-12}")
    except ImportError:
        print("numpy unavailable -- float cross-check skipped (plays no role in the proof)")

if __name__ == "__main__":
    main()
