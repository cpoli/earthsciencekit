# EarthScienceKit 0.2.0 verification

## Verified release scope

The release contains seven domain subpackages: plate_tectonics (23 models),
seismology (23), atmosphere (22), climate (24), oceanography (24),
geophysics (31) and volcanology (24). These comprise 162 algebraic models and
nine dynamic models. Shared integration, results, visualization, units and
constants support them.

## Results

Verified on 2026-10-17 with Python 3.14 before the first push.

- Workspace: 2003 tests pass with warnings treated as errors.
- Production statement coverage: 100% (1579/1579 statements, zero misses).
- Distribution metadata and earthsciencekit.__version__ both equal 0.2.0.
- Ruff lint/format checks pass. Shared results/integration Mypy checks pass with
  check-untyped-defs, ignore-missing-imports and silent dependency following.
- Sphinx HTML builds from a clean gallery with warnings treated as errors; all
  114 scripts execute, and each has a JupyterLite button whose notebook starts
  with `%pip install earthsciencekit`.
- Sphinx doctests: nine tests pass with zero failures.
- Source distribution and wheel build and pass `twine check --strict`.
- The code cells of notebooks/quickstart.ipynb run against the checkout.

## Scientific evidence

Fixed-value checks cover every algebraic model. Physical bounds, invalid
real/finite inputs, broadcast arrays and cross-parameter contradictions are
exercised. Each domain's test_physics.py checks identities and limits, such as Euler-vector
tangency, acoustic interface energy conservation, moment-magnitude inversion and
shallow/deep wave limits.

Dynamic checks in each domain's test_physics.py compare sampled histories with
analytic responses, conserved quantities or limiting cases. Results are tested for
copied read-only arrays, labels, summaries, figure axes and real Pillow GIF export.

The 2003 count includes parameterized invalid-input cases; it is not 2003
distinct scientific benchmarks. Guides and primary university/agency references explain
physical scope. Bounds cannot establish whether a model applies to observations.

## Coverage policy

All production modules, private validation, shared tools and plotting utilities
are included. Only test modules are omitted. No production coverage-ignore
pragmas were introduced. fail_under=100 is enforced in pyproject.toml and CI.
This is statement coverage, not branch coverage or a proof of scientific accuracy.

## Reproduce

```bash
python -m pip install -e '.[dev]'
MPLBACKEND=Agg python -m pytest --cov=earthsciencekit --cov-report=term-missing
python -m ruff check .
python -m ruff format --check .
python -m mypy earthsciencekit/results.py earthsciencekit/integrators --ignore-missing-imports --follow-imports=silent --check-untyped-defs
MPLBACKEND=Agg python -m sphinx -E -W --keep-going -b html docs/source docs/_build/html
python -m sphinx -W -b doctest docs/source docs/_build/doctest
python -m build
```

Local verification used Python 3.14. The CI workflow runs Python 3.10-3.13;
those interpreters were not run locally. Documentation structure follows
MathematicsKit reference commit 7c075588a3208f5fdaed9cb7dfd41bf0cd01e9e4 and
the EngineeringKit adaptation: Subpackages, History, Examples and API.
