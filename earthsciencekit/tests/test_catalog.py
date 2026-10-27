"""model_catalog.json agrees with the code: names per domain and golden values."""

import ast
import importlib
import json
import operator
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pytest

import earthsciencekit
from earthsciencekit import constants
from earthsciencekit._validation import model

CATALOG_PATH = Path(__file__).resolve().parents[2] / "model_catalog.json"
# The catalog lives at the repository root and is not shipped in the sdist.
pytestmark = pytest.mark.skipif(not CATALOG_PATH.exists(), reason="model_catalog.json is only in the repository")
CATALOG = json.loads(CATALOG_PATH.read_text()) if CATALOG_PATH.exists() else {}
MODEL_CODE = model(lambda: None).__code__
OPERATORS: dict[type, Callable[..., Any]] = {
    ast.Add: operator.add,
    ast.Sub: operator.sub,
    ast.Mult: operator.mul,
    ast.Div: operator.truediv,
    ast.Pow: operator.pow,
    ast.USub: operator.neg,
    ast.UAdd: operator.pos,
}
NUMPY_NAMES = {"pi", "e", "sqrt", "log", "log10", "exp", "sin", "cos", "tan", "arctan", "sinh", "cosh", "tanh"}


def evaluate(text):
    """Evaluate a catalog ``expected`` string: numbers, arithmetic, ``np.<name>`` and package constants only."""

    def visit(node):
        if isinstance(node, ast.Expression):
            return visit(node.body)
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)) and not isinstance(node.value, bool):
            return node.value
        if isinstance(node, ast.List):
            return [visit(item) for item in node.elts]
        if isinstance(node, ast.BinOp) and type(node.op) in OPERATORS:
            return OPERATORS[type(node.op)](visit(node.left), visit(node.right))
        if isinstance(node, ast.UnaryOp) and type(node.op) in OPERATORS:
            return OPERATORS[type(node.op)](visit(node.operand))
        if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name) and node.value.id == "np" and node.attr in NUMPY_NAMES:
            return getattr(np, node.attr)
        if isinstance(node, ast.Name) and node.id.isupper() and hasattr(constants, node.id):
            return getattr(constants, node.id)
        if isinstance(node, ast.Call) and not node.keywords:
            function = visit(node.func)
            assert callable(function)
            return function(*(visit(argument) for argument in node.args))
        raise ValueError(f"unsupported catalog expression {text!r}")

    return np.asarray(visit(ast.parse(text, mode="eval")), dtype=float)


ENTRIES = [pytest.param(domain, entry, id=f"{domain}.{entry['name']}") for domain, entries in CATALOG.items() for entry in entries]


@pytest.mark.parametrize("domain", list(CATALOG))
def test_catalog_names_match_algebraic_all(domain):
    module = importlib.import_module(f"earthsciencekit.{domain}")
    algebraic = [name for name in module.__all__ if getattr(module, name).__code__ is MODEL_CODE]
    names = [entry["name"] for entry in CATALOG[domain]]
    assert len(names) == len(set(names))
    assert sorted(names) == sorted(algebraic)


def test_catalog_covers_every_domain():
    assert sorted(CATALOG) == sorted(earthsciencekit.__all__)


@pytest.mark.parametrize(("domain", "entry"), ENTRIES)
def test_catalog_reference_reproduces_expected(domain, entry):
    function = getattr(importlib.import_module(f"earthsciencekit.{domain}"), entry["name"])
    assert [parameter["name"] for parameter in entry["params"]] == list(function.__wrapped__.__code__.co_varnames[: len(entry["params"])])
    np.testing.assert_allclose(function(*entry["args"]), evaluate(entry["expected"]), rtol=1e-11, atol=1e-12)
