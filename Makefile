# Makefile for the Erdos minimum-overlap upper-bound note.
#
# `make verify` is the theorem-grade check: it runs the stdlib-only exact
# certifier on the n=512 construction, RECOMPUTING Q from the construction
# vector and printing it (it does NOT compare against a hard-coded Q). A
# passing run therefore establishes mu <= Q < 0.3808622032020279475140496,
# not merely that some pre-stored number was echoed back.

PYTHON      ?= python3
CONSTRUCTION = certs/erdos_dc_n512.json
CERTIFIER    = scripts/erdos_cert_general.py
TEX          = erdos_minimum_overlap_bracket.tex
FIGURE       = fig_construction.pdf

.PHONY: all verify figure pdf clean

all: verify

## verify: exact rational certification of the upper bound (stdlib only).
##   Recomputes Q = 117871142698558740618278313/309485009821345068724781056
##   from certs/erdos_dc_n512.json and prints mu <= Q < 0.3808622032020279475140496.
verify:
	$(PYTHON) $(CERTIFIER) $(CONSTRUCTION)

## figure: regenerate Figure 1 (needs matplotlib + numpy).
figure:
	$(PYTHON) make_figure.py

## pdf: regenerate the figure and typeset the note.
pdf: figure
	tectonic $(TEX)

clean:
	rm -f erdos_minimum_overlap_bracket.aux erdos_minimum_overlap_bracket.log \
	      erdos_minimum_overlap_bracket.out
