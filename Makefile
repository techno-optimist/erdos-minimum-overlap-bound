# Makefile for the Erdos minimum-overlap upper-bound note.
#
# `make verify` is the theorem-grade check.  It runs the stdlib-only exact
# certifier on every certified construction, RECOMPUTING each bound from its
# construction vector and printing it (it never compares against a hard-coded
# number).  Two parts:
#   verify-v11     -- recomputes the v1.1 bound Q < 0.3808622032020279475140496
#                     from the n=512 construction (theorem intact);
#   verify-leaders -- recomputes the v1.2 bounds for the current best arena
#                     constructions (Hyra n=1024, lnzwz n=512), exactly,
#                     including the minimal sub-ULP admissibility repair.
# `make selftest` re-derives all of them and asserts they equal the golden
# rationals published in the note.

PYTHON      ?= python3
CONSTRUCTION = certs/erdos_dc_n512.json
CERTIFIER    = scripts/erdos_cert_general.py
LEADERS      = scripts/certify_leaders.py
SELFTEST     = scripts/selftest.py
TEX          = erdos_minimum_overlap_bracket.tex
FIGURE       = fig_construction.pdf

.PHONY: all verify verify-v11 verify-leaders selftest figure pdf clean

all: verify selftest

## verify: exact rational certification of every certified bound (stdlib only).
verify: verify-v11 verify-leaders

## verify-v11: recompute the v1.1 headline bound Q from certs/erdos_dc_n512.json.
##   Prints mu <= Q < 0.3808622032020279475140496 (the v1.1 theorem, intact).
verify-v11:
	$(PYTHON) $(CERTIFIER) $(CONSTRUCTION)

## verify-leaders: recompute the v1.2 bounds for the current best arena
##   constructions (Hyra n=1024 admissible-outright; lnzwz n=512 after the
##   minimal sub-ULP admissibility repair).  Recompute-not-echo.
verify-leaders:
	$(PYTHON) $(LEADERS)

## selftest: re-derive all certified bounds and assert they equal the golden
##   rationals in the note (exit 0 iff all pass).
selftest:
	$(PYTHON) $(SELFTEST)

## figure: regenerate Figure 1 (needs matplotlib + numpy).
figure:
	$(PYTHON) make_figure.py

## pdf: regenerate the figure and typeset the note.
pdf: figure
	tectonic $(TEX)

clean:
	rm -f erdos_minimum_overlap_bracket.aux erdos_minimum_overlap_bracket.log \
	      erdos_minimum_overlap_bracket.out
