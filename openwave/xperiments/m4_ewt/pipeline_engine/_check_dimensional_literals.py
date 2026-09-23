#!/usr/bin/env python3
"""
Scan physics/*.py for hardcoded dimensional literals.

Every dimensional constant the engine needs comes from a UnitSystem
feature. A file that assigns a numeric literal to a known dimensional
name freezes that value to whichever unit system the author had in mind
at the time, and the value silently stops matching the unit system the
pipeline actually runs under.

The linter catches three patterns:

    c = 1.0                     plain assignment
    c: float = 1.0              annotated assignment
    def f(c=1.0): ...           argument default

Assignments whose value is an expression (units.c, a formula, a call)
are not flagged. The rule is about literals, not about using the name.

Exempt files: units.py (defines the constants), __init__.py, and any
file whose name starts with `_test_` or `test_`.

Usage:
    python -m openwave.xperiments.m4_ewt.pipeline_engine._check_dimensional_literals

Exit code 0 if clean, 1 if any offending assignment is found.
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path

# Names that must never be assigned a numeric literal in physics/.
# Kept in sync with the UnitSystem contract in physics/units.py.
_FORBIDDEN = frozenset(
    {
        "c",
        "wavelength",
        "dx",
        "dt",
        "gamma",
        "rho_0",
        "A_pi",
        "eps_M",
        "N_geom",
        "X_eff",
        "N_nu_eff",
    }
)

_EXEMPT_FILES = frozenset({"units.py", "__init__.py"})
_EXEMPT_PREFIXES = ("_test_", "test_")


def _is_numeric_literal(node: ast.AST) -> bool:
    """
    True for a bare int or float literal, False for bool (which is a
    subclass of int in Python) and for every expression.
    """
    if not isinstance(node, ast.Constant):
        return False
    if isinstance(node.value, bool):
        return False
    return isinstance(node.value, (int, float))


def _check_assignment(node: ast.AST, path: Path) -> list[str]:
    issues: list[str] = []

    if isinstance(node, ast.Assign):
        for target in node.targets:
            if isinstance(target, ast.Name) and target.id in _FORBIDDEN:
                if _is_numeric_literal(node.value):
                    issues.append(
                        f"{path}:{node.lineno}: assign to '{target.id}' "
                        f"= {ast.unparse(node.value)}"
                    )

    elif isinstance(node, ast.AnnAssign):
        target = node.target
        if (
            isinstance(target, ast.Name)
            and target.id in _FORBIDDEN
            and node.value is not None
            and _is_numeric_literal(node.value)
        ):
            issues.append(
                f"{path}:{node.lineno}: annotated assign to "
                f"'{target.id}' = {ast.unparse(node.value)}"
            )

    elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
        args = node.args

        # Positional args (posonly + normal) share the same defaults list,
        # which applies to the last len(defaults) of them.
        positional = args.posonlyargs + args.args
        n_with_default = len(args.defaults)
        defaulted = positional[len(positional) - n_with_default :] if n_with_default else []
        for arg, default in zip(defaulted, args.defaults):
            if arg.arg in _FORBIDDEN and _is_numeric_literal(default):
                issues.append(
                    f"{path}:{node.lineno}: default for argument "
                    f"'{arg.arg}' = {ast.unparse(default)}"
                )

        # Keyword-only args have a parallel defaults list where None
        # means "no default".
        for arg, default in zip(args.kwonlyargs, args.kw_defaults):
            if arg.arg in _FORBIDDEN and default is not None and _is_numeric_literal(default):
                issues.append(
                    f"{path}:{node.lineno}: default for argument "
                    f"'{arg.arg}' = {ast.unparse(default)}"
                )

    return issues


def _check_file(path: Path) -> list[str]:
    try:
        source = path.read_text(encoding="utf-8")
    except UnicodeDecodeError as e:
        return [f"{path}: cannot decode: {e}"]

    try:
        tree = ast.parse(source)
    except SyntaxError as e:
        return [f"{path}: syntax error: {e}"]

    issues: list[str] = []
    for node in ast.walk(tree):
        issues.extend(_check_assignment(node, path))
    return issues


def main() -> int:
    physics_dir = Path(__file__).resolve().parent / "physics"
    if not physics_dir.exists():
        print(
            f"physics directory not found at {physics_dir}",
            file=sys.stderr,
        )
        return 1

    issues: list[str] = []
    checked = 0
    for path in sorted(physics_dir.glob("*.py")):
        name = path.name
        if name in _EXEMPT_FILES:
            continue
        if name.startswith(_EXEMPT_PREFIXES):
            continue
        checked += 1
        issues.extend(_check_file(path))

    if issues:
        print(
            f"Hardcoded dimensional literals found "
            f"({len(issues)} issue(s) in {checked} file(s)):",
            file=sys.stderr,
        )
        for issue in issues:
            print(f"  {issue}", file=sys.stderr)
        return 1

    print(f"physics/: no hardcoded dimensional literals ({checked} files checked)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
