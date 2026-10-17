# Contributing

Install `python -m pip install -e '.[dev]'`. Use Python 3.10 or newer.

For every public model, document the governing expression, every input/output
unit, physical conventions and applicability limits. Validate mathematical bounds
and finite real inputs. Test independent known values and at least one meaningful
identity or limit. Add a runnable example and update the guide/history when the
model extends a scientific concept. Avoid runtime formula evaluation or hidden
network requests.

Run pytest with production coverage, Ruff lint/format and Sphinx warnings-as-errors
before proposing a change. CI enforces 100% statement coverage; this is a release
check, not proof of scientific accuracy. Keep version changes explicit and update
CHANGELOG.md. Releases are cut by pushing a `vX.Y.Z` tag (see
`.github/workflows/release.yml`).
