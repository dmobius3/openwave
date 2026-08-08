#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Validate every xparameter module in M4:
  1. Compiles (py_compile) – catches syntax errors without importing.
  2. Imports and checks the structure of XPARAMETERS.

Usage (from any directory):
    python openwave/xperiments/m4_ewt/research/scripts/validate_xparameters.py

Or as a module (preferred – avoids sys.path manipulation):
    cd <project_root>
    python -m openwave.xperiments.m4_ewt.research.scripts.validate_xparameters
"""

from __future__ import annotations

import importlib
import py_compile
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Project-root discovery
# ---------------------------------------------------------------------------
# File lives at: <root>/openwave/xperiments/m4_ewt/research/scripts/<file>
#   parents[0] = scripts/
#   parents[1] = research/
#   parents[2] = m4_ewt/
#   parents[3] = xperiments/
#   parents[4] = openwave/          ← the package
#   parents[5] = <project_root>/    ← must be on sys.path
_THIS_FILE = Path(__file__).resolve()
_PROJECT_ROOT = _THIS_FILE.parents[5]  # contains the `openwave` package

if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

# Base import path for xparameter modules
_IMPORT_BASE = "openwave.xperiments.m4_ewt.xparameters"

# Required top-level keys in every XPARAMETERS dict
_REQUIRED_KEYS = ("engine",)

# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------


def discover_xparameter_modules() -> list[tuple[str, Path]]:
    """Return sorted (stem, path) tuples for all xparameter modules."""
    xp_dir = _THIS_FILE.parents[2] / "xparameters"

    if not xp_dir.exists():
        print(f"[ERROR] xparameters directory not found: {xp_dir}")
        sys.exit(1)

    modules: list[tuple[str, Path]] = []
    for f in sorted(xp_dir.glob("*.py")):
        if f.name.startswith("__"):
            continue
        try:
            content = f.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        if "XPARAMETERS" in content:
            modules.append((f.stem, f))

    return modules


# ---------------------------------------------------------------------------
# Individual checks
# ---------------------------------------------------------------------------


def check_compile(path: Path) -> tuple[bool, str]:
    """Syntax-check the file without importing it."""
    try:
        py_compile.compile(str(path), doraise=True)
        return True, "OK"
    except py_compile.PyCompileError as exc:
        return False, f"SyntaxError: {exc}"


def check_import(stem: str) -> tuple[bool, str]:
    """Import the module and verify XPARAMETERS structure."""
    fqn = f"{_IMPORT_BASE}.{stem}"
    try:
        mod = importlib.import_module(fqn)
    except ModuleNotFoundError as exc:
        missing = str(exc)
        if "openwave" in missing:
            return False, (
                f"ModuleNotFoundError: {exc}  "
                f"[hint: run from project root or add {_PROJECT_ROOT} to PYTHONPATH]"
            )
        return False, f"ModuleNotFoundError: {exc}"
    except Exception as exc:  # noqa: BLE001
        return False, f"ImportError: {type(exc).__name__}: {exc}"

    # Structure checks
    xp = getattr(mod, "XPARAMETERS", None)
    if xp is None:
        return False, "XPARAMETERS not found in module"
    if not isinstance(xp, dict):
        return False, f"XPARAMETERS is {type(xp).__name__}, expected dict"
    missing_keys = [k for k in _REQUIRED_KEYS if k not in xp]
    if missing_keys:
        return False, f"XPARAMETERS missing keys: {missing_keys}"

    return True, "OK"


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------


def main() -> None:
    modules = discover_xparameter_modules()
    if not modules:
        print("No xparameter modules found.")
        sys.exit(1)

    col_w = max(len(s) for s, _ in modules) + 2
    passed = 0

    print(f"\n{'Module':<{col_w}}  {'Compile':<10}  {'Import'}")
    print("-" * (col_w + 50))

    for stem, path in modules:
        ok_compile, msg_compile = check_compile(path)

        if ok_compile:
            ok_import, msg_import = check_import(stem)
        else:
            ok_import, msg_import = False, "skipped"

        ok = ok_compile and ok_import
        status = "PASS" if ok else "FAIL"

        compile_cell = "OK" if ok_compile else f"FAIL ({msg_compile})"
        import_cell = msg_import if ok_import else f"FAIL: {msg_import}"

        print(f"  {stem:<{col_w}} [{status}]  compile: {compile_cell:<10}  import: {import_cell}")
        if ok:
            passed += 1

    print(f"\n{'─'*60}")
    print(f"Result: {passed}/{len(modules)} passed")
    if passed < len(modules):
        print(f"        {len(modules)-passed} failed")
        print(f"\nProject root used: {_PROJECT_ROOT}")
        print("If imports fail, verify that all __init__.py files exist:")
        for part in [
            "openwave",
            "openwave/xperiments",
            "openwave/xperiments/m4_ewt",
            "openwave/xperiments/m4_ewt/xparameters",
        ]:
            p = _PROJECT_ROOT / part / "__init__.py"
            status = "✓" if p.exists() else "✗ MISSING"
            print(f"  {status}  {p.relative_to(_PROJECT_ROOT)}")
    print()
    sys.exit(0 if passed == len(modules) else 1)


if __name__ == "__main__":
    main()
