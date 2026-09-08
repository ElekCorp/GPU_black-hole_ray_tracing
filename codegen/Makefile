# Metric code generation and verification.
#
# The generated headers are committed, so building the ray tracer needs no
# Python at all.  This Makefile is only for regenerating them, which also
# re-runs every proof obligation: a header is written only if its obligations
# hold, so a checked-in header is one that passed.
#
#   make            regenerate every metric (runs all proofs; slow)
#   make fast       skip the Gappa round-off stage
#   make check      regenerate into a scratch dir and diff against the committed
#                   headers - use this in CI to catch a stale or hand-edited header
#   make venv       create the Python environment

PYTHON  ?= $(HOME)/.venv-bhcodegen/bin/python
OUTDIR  ?= ../OpenACC_cli/generated
METRICS ?= minkowski_cartesian schwarzschild kerr_newman

all:
	$(PYTHON) -m bhgen.cli $(METRICS) -o $(OUTDIR) --report verification-report.json

fast:
	$(PYTHON) -m bhgen.cli $(METRICS) -o $(OUTDIR) --no-roundoff

check:
	@rm -rf /tmp/bhgen-check && mkdir -p /tmp/bhgen-check
	$(PYTHON) -m bhgen.cli $(METRICS) -o /tmp/bhgen-check --no-roundoff >/dev/null
	@for m in $(METRICS); do \
	  diff -q $(OUTDIR)/bh_$$m_metric.h /tmp/bhgen-check/bh_$$m_metric.h >/dev/null 2>&1 \
	    || diff -u $(OUTDIR)/bh_$$m'_metric.h' /tmp/bhgen-check/bh_$$m'_metric.h' | head -40; \
	done
	@echo "committed headers match a fresh generation"

venv:
	python3 -m venv $(HOME)/.venv-bhcodegen
	$(HOME)/.venv-bhcodegen/bin/pip install -q sympy mpmath numpy
	@echo "Gappa is also required for 'make all': apt-get install gappa"

.PHONY: all fast check venv
