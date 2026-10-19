# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

EarthScienceKit (v0.3.0, on PyPI) is an educational Python library of inspectable Earth-system models: nine domain subpackages (`plate_tectonics`, `seismology`, `atmosphere`, `climate`, `oceanography`, `geophysics`, `volcanology`, `hydrology`, `glaciology`) holding 237 algebraic models and eleven sampled dynamic models; 0.2.0 shipped the first seven and 0.3.0 added `hydrology` and `glaciology`. It follows the "MathematicsKit family" structure (see `docs/README.md` for the reference commit). Runtime deps are only NumPy, SciPy and Matplotlib; Python >= 3.10.

## Commands

```bash
python -m pip install -e '.[dev]'

# Full test suite with the coverage gate (fail_under=100, enforced in CI)
MPLBACKEND=Agg python -m pytest --cov=earthsciencekit --cov-report=term-missing

# Single domain / single test
MPLBACKEND=Agg python -m pytest earthsciencekit/climate/tests/test_models.py
MPLBACKEND=Agg python -m pytest earthsciencekit/climate/tests/test_models.py -k absorbed_solar

python -m ruff check .
python -m ruff format --check .
python -m mypy earthsciencekit/results.py earthsciencekit/integrators --ignore-missing-imports --follow-imports=silent --check-untyped-defs

# Docs: warnings are errors and all gallery scripts are executed
cd docs && make html
make -C docs doctest

python -m build
```

CI (`.github/workflows/ci.yml`) runs pytest+coverage on Python 3.10–3.13, ruff check/format, the mypy check above and `build` + `twine check`; `docs.yml` builds the Sphinx HTML and doctests. Pushing a `vX.Y.Z` tag runs `release.yml`, which checks the tag against `__version__`, `CITATION.cff` and a `## [X.Y.Z] - date` CHANGELOG section, then publishes to PyPI (trusted publishing), creates the GitHub Release and deploys the docs to `gh-pages`. After Zenodo archives the release, add its version DOI (never the concept DOI) to `CITATION.cff` and the README badge. Pytest runs with `filterwarnings=["error"]`, so any warning fails a test. Ruff line length is 140; `E741` is ignored to allow engineering symbols like `I` and `l`.

## Architecture

**Family layout**: each domain `earthsciencekit/<domain>/` has `systems/<topic>.py` (the models, grouped by topic, each module with its own `__all__`), `core/base.py` (its result type, the shared `TimeSeries`), `visualizers/` (its plots, the shared `plot_series`/`animate_series`) and `tests/` (`test_models.py` for golden values and contracts, `test_physics.py` for identities, limits and dynamic checks). The domain `__init__.py` only re-exports the systems functions with the public `__all__`, so `from earthsciencekit.<domain> import name` is the stable API; add a new model to a systems module, its `__all__`, the domain `__init__.py` imports and `__all__`. `earthsciencekit/tests/` holds the cross-cutting tests: `test_core.py` (validation, integrators, results, units, visualizers, the layout contract) and `test_relations.py` (cross-parameter rejections). Pytest `testpaths` is `earthsciencekit`.

**Algebraic model pattern** — every public function follows the same shape, and new ones should match it exactly:
- Decorated with `@model` from `_validation.py`, which runs the body under `np.errstate(raise)`, converts `FloatingPointError` to `ValueError`, and rejects nonfinite results.
- Each input passes through `quantity(name, value, lower=..., upper=..., strict_lower=..., strict_upper=...)`, which rejects booleans, complex, strings, empty and nonfinite input, and returns a float array (scalars stay 0-d, arrays broadcast). Cross-parameter constraints use `require(condition, message)`.
- NumPy-style docstring with a fixed template: one-line title, a "why" sentence, Parameters giving unit + constraint kind, Returns with output unit, and Notes containing the applicability limit and `Equation: ``...```.
- Name added to the module's `__all__`.
- Tests: a `test_<name>_reference` golden value checked with `rtol=1e-11`, and a parametrized `test_<name>_rejects_invalid` over `[nan, inf, -inf, [], "bad", True, 1j]`.

`model_catalog.json` at the repo root holds the per-model spec (name, title, params with unit/kind, expr, reference args/expected, why, limit, output_unit) for each domain, and mirrors the docstrings and golden tests. Nothing in the repo reads it at runtime or build time, so keep it in sync by hand when you add or change a model.

**Dynamic models** (`simulate_*`-style functions, one per domain) take scalar parameters plus a time vector in seconds validated by `time_grid`, then either call `integrators.integrate_ode` (a validated wrapper over `scipy.integrate.solve_ivp`) or compute samples directly. In both cases they return `results.TimeSeries`, a frozen dataclass whose `time`/`values` arrays are copied and set read-only, with `.summary()`, `.plot()` and `.animate()` (lazy-imported from `visualizers`). They are appended to `__all__` after definition and use `@safe_arithmetic` rather than `@model` where the return value is not an array.

**Conventions**: SI units throughout; angles and latitude in radians. Conversions are explicit through `units.convert(value, from_unit, to_unit)`, never implicit. Shared physical constants live in `earthsciencekit.constants`. Runtime formula evaluation and network access are not allowed.

**Docs**: `docs/source/conf.py` owns the `SUBPACKAGES` registry. At `builder-inited` it generates the domain hubs, card grids and nav into `docs/source/_generated/`, and Sphinx-Gallery executes `examples/<domain>/plot_*.py` into `docs/source/api/gallery/<domain>/`. Both output directories are gitignored, and any failing example aborts the build. The handwritten pages are `docs/source/guides/<domain>.rst` (equations, conventions, checks), `docs/source/history/<domain>_breakthroughs.rst`, `docs/source/api/<domain>.rst` and `docs/source/limitations.rst`. Adding a domain means updating `SUBPACKAGES` as well as adding those pages and an examples directory.

## Contribution expectations (from CONTRIBUTING.md)

For every public model, document the governing expression, all units, conventions and applicability limits. Test an independent known value plus at least one meaningful identity or limit, and add a runnable example. Update the domain guide and history when the model extends a scientific concept. Production statement coverage must stay at 100% without `# pragma: no cover`, and test files are the only exclusion. Record version changes in `CHANGELOG.md`, and keep `__version__` in `earthsciencekit/__init__.py` in sync with `pyproject.toml` (`tests/test_core.py::test_version` asserts it). `VERIFICATION.md` records exact release results (test counts, coverage), so update it when you re-verify a release.

# Output Guidelines
- **Be Concise:** Provide direct code and answers first. Omit setup text, fluff, and conversational responses.
- **Code Generation:** Output only updated code blocks, functions, or unified diffs. Never output full files unless instructed.
- **Explanations:** Limit inline code comments. Explain high-level logic in 1–2 bullet points after the code block.
- **Format:** Use bulleted lists and tables for comparisons instead of paragraph blocks.