# Documentation architecture

The structure follows MathematicsKit at reference commit
7c075588a3208f5fdaed9cb7dfd41bf0cd01e9e4, using the EngineeringKit adaptation:
four entry points, grouped domain cards, generated hubs and History/Examples/API
navigation. conf.py owns the five-domain registry. Generated files are excluded
from source control and recreated on build.

Handwritten source lives in source/history, source/guides and source/api.
Gallery scripts live in ../examples/<domain>; Sphinx-Gallery writes
source/api/gallery/<domain>, including thumbnails, figures, source and notebooks.
Failures abort the build. API links and guide equations expose the same models.
Each example has a JupyterLite button (jupyterlite-sphinx) whose notebook starts
with `%pip install earthsciencekit`, so it runs the PyPI release in the browser.

Build from the project root (or run `make html` inside docs/):

```bash
python -m pip install -e '.[dev]'
cd docs && make html
```

Open docs/build/html/index.html. The offline ZIP starts at index.html.
