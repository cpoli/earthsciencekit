# Contributing

Install `python -m pip install -e '.[dev]'`. Use Python 3.10 or newer.

For every public model, document the governing expression, every input/output
unit, physical conventions and applicability limits. Validate mathematical bounds
and finite real inputs. Test independent known values and at least one meaningful
identity or limit. Add a runnable example and update the guide/history when the
model extends a scientific concept. Avoid runtime formula evaluation or hidden
network requests.

Run pytest with production coverage, Ruff lint/format and Sphinx warnings-as-errors
before proposing a change. CI enforces 100% statement and branch coverage and a
whole-package mypy check, runs the docstring `Examples` doctests, and
`earthsciencekit/tests/test_catalog.py` keeps `model_catalog.json` in step
with the code. These are release checks, not proof of scientific accuracy.
Give each new algebraic model a one-line doctested example. Keep version
changes explicit and update CHANGELOG.md. Releases are cut by pushing a
`vX.Y.Z` tag (see `.github/workflows/release.yml`).
