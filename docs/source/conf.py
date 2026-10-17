"""MathematicsKit-style domain hubs, cross-links and executed example galleries."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from earthsciencekit import __version__

project = "earthsciencekit"
author = "Charles Poli"
release = __version__
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinx.ext.doctest",
    "sphinx.ext.viewcode",
    "sphinx_autodoc_typehints",
    "sphinx_gallery.gen_gallery",
    "sphinx_design",
    "jupyterlite_sphinx",
    "sphinx_sitemap",
    "sphinxext.opengraph",
]
napoleon_google_docstring = False
napoleon_numpy_docstring = True
autodoc_default_options = {"members": True, "undoc-members": False, "show-inheritance": True}
autodoc_typehints = "description"
templates_path = ["_templates"]
html_static_path = ["_static"]
html_css_files = ["earthsciencekit.css"]
html_theme = "pydata_sphinx_theme"
# sphinx_gallery_conf holds a function (the JupyterLite hook), so Sphinx cannot pickle it.
suppress_warnings = ["config.cache"]

# Search-engine metadata: sphinx-sitemap writes sitemap.xml with absolute
# URLs from html_baseurl, and sphinxext-opengraph adds a <meta
# name="description"> plus og:* link-preview tags from each page's text.
# html_baseurl also gives every page a <link rel="canonical">.
html_baseurl = "https://cpoli.github.io/earthsciencekit/"
html_title = project
sitemap_url_scheme = "{link}"
sitemap_excludes = [
    "search.html",
    "genindex.html",
    "py-modindex.html",
    "*sg_execution_times.html",
    "_modules/*",
    # Include-only fragments (card grids, nav strips), not standalone pages.
    "_generated/grid_*",
    "_generated/vars.html",
    "_generated/nav/*",
]
ogp_site_url = html_baseurl
ogp_site_name = project
ogp_image = html_baseurl + "_static/images/readme_hero.png"
ogp_social_cards = {"enable": False}
html_theme_options = {
    "github_url": "https://github.com/cpoli/earthsciencekit",
    "icon_links": [
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/earthsciencekit/",
            "icon": "fa-brands fa-python",
            "type": "fontawesome",
        },
    ],
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "show_toc_level": 2,
    "navigation_with_keys": True,
    "navigation_depth": 2,
}
exclude_patterns = ["_build", "jupyterlite_contents", "_contents", "_generated/nav/*", "_generated/grid_*.rst", "_generated/vars.rst"]

SUBPACKAGES = [
    {
        "name": "plate_tectonics",
        "category": "Solid Earth",
        "blurb": "Euler plate motion, spreading, subduction, isostasy, flexure and lithosphere cooling.",
    },
    {
        "name": "seismology",
        "category": "Solid Earth",
        "blurb": "Elastic wave speeds, reflection, seismic moment, source pulses and damped response.",
    },
    {
        "name": "atmosphere",
        "category": "Fluid Earth and Climate",
        "blurb": "Hydrostatic balance, dry thermodynamics, moisture, rotation and parcel ascent.",
    },
    {
        "name": "climate",
        "category": "Fluid Earth and Climate",
        "blurb": "Radiative balance, forcing, heat storage, carbon baselines and climate adjustment.",
    },
    {
        "name": "oceanography",
        "category": "Fluid Earth and Climate",
        "blurb": "Tides, waves, seawater density, geostrophy, the Ekman spiral, wind-driven gyres and thermohaline boxes.",
    },
    {
        "name": "volcanology",
        "category": "Solid Earth",
        "blurb": (
            "Radiometric dating, isotope and trace-element geochemistry, magma viscosity and ascent, "
            "ground deformation and eruption columns."
        ),
    },
    {
        "name": "geophysics",
        "category": "Solid Earth",
        "blurb": "Gravity anomalies, the figure of the Earth, the geoid, the geomagnetic field, paleomagnetism and the dynamo.",
    },
]


def _jupyterlite_install_cell(notebook_content, notebook_filename):
    """Prepend a ``%pip install`` cell so the Pyodide kernel fetches earthsciencekit from PyPI."""
    notebook_content["cells"].insert(
        0,
        {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": ["%pip install earthsciencekit"]},
    )


sphinx_gallery_conf = {
    "examples_dirs": [f"../../examples/{s['name']}" for s in SUBPACKAGES],
    "gallery_dirs": [f"api/gallery/{s['name']}" for s in SUBPACKAGES],
    "filename_pattern": r"/plot_",
    "download_all_examples": True,
    "within_subsection_order": "FileNameSortKey",
    "backreferences_dir": "gen_modules/backreferences",
    "doc_module": ("earthsciencekit",),
    "abort_on_example_error": True,
    "matplotlib_animations": (True, "jshtml"),
    # "Launch JupyterLite" button: runs the notebook in the browser on Pyodide,
    # which has NumPy, SciPy and Matplotlib, with no server.
    "jupyterlite": {"use_jupyter_lab": True, "notebook_modification_function": _jupyterlite_install_cell},
}


def _cards(entries):
    text = ".. grid:: 1 2 3 3\n   :gutter: 2\n\n"
    for title, target, description in entries:
        text += f"   .. grid-item-card:: {title}\n      :link: {target}\n      :link-type: doc\n\n      {description}\n\n"
    return text


def _generate_domain_docs(app):
    generated = Path(app.srcdir) / "_generated"

    def write(relative, text):
        target = generated / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text)

    write("vars.rst", f".. |num_subpackages| replace:: {len(SUBPACKAGES)}\n")
    targets = {
        "hub": lambda n: f"/_generated/subpackages/{n}",
        "history": lambda n: f"/history/{n}_breakthroughs",
        "examples": lambda n: f"/api/gallery/{n}/index",
        "api": lambda n: f"/api/{n}",
    }
    categories = list(dict.fromkeys(s["category"] for s in SUBPACKAGES))
    for kind, target in targets.items():
        chunks = []
        for category in categories:
            chunks.append(f".. rubric:: {category}\n\n")
            chunks.append(
                _cards([(f"earthsciencekit.{s['name']}", target(s["name"]), s["blurb"]) for s in SUBPACKAGES if s["category"] == category])
            )
        write(f"grid_{kind}.rst", "\n".join(chunks))
    for domain in SUBPACKAGES:
        name = domain["name"]
        title = f"earthsciencekit.{name}"
        entries = [
            ("History", targets["history"](name), "Foundational ideas linked to the implementation."),
            ("Examples", targets["examples"](name), "Executed scripts, plots, source and notebook downloads."),
            ("API reference", targets["api"](name), "Public classes, functions and result objects."),
        ]
        write(
            f"subpackages/{name}.rst",
            f"{title}\n{'=' * len(title)}\n\n{domain['blurb']}\n\n"
            + _cards(entries)
            + f"\n:doc:`Units, assumptions and numerical validation </guides/{name}>`.\n",
        )
        nav = [
            f":doc:`Domain hub <{targets['hub'](name)}>`",
            f":doc:`History <{targets['history'](name)}>`",
            f":doc:`Examples <{targets['examples'](name)}>`",
            f":doc:`API reference <{targets['api'](name)}>`",
        ]
        write(f"nav/{name}.rst", ".. container:: subpkg-nav\n\n   " + " · ".join(nav) + "\n")


def setup(app):
    app.connect("builder-inited", _generate_domain_docs)
