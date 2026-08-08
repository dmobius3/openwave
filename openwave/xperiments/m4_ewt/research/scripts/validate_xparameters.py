#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Validate every xparameter module in M4:
  1. Compiles (py_compile) – catches syntax errors without importing.
  2. Imports and checks the structure of XPARAMETERS.

Usage:
    python research/scripts/validate_xparameters.py
"""

import sys
import importlib
import py_compile
from pathlib import Path

# ---------------------------------------------------------------------------
# Discovery (mirrors _launcher.py non-recursive glob)
# ---------------------------------------------------------------------------

def discover_xparameter_modules():
    """Return sorted list of (stem, path) tuples for modules under xparameters/
    that define XPARAMETERS."""
    xp_dir = Path(__file__).resolve().parents[2] / "xparameters"
    if not xp_dir.exists():
        print(f"xparameters directory not found: {xp_dir}")
        sys.exit(1)

    modules = []
    for f in sorted(xp_dir.glob("*.py")):
        if f.name == "__init__.py":
            continue
        # Read safely with UTF-8
        try:
            content = f.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            continue
        if "XPARAMETERS" in content:
            modules.append((f.stem, f))
    return modules


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate_compile(path):
    """Return (ok, message) after attempting py_compile."""
    try:
        py_compile.compile(str(path), doraise=True)
        return True, "compiles"
    except py_compile.PyCompileError as exc:
        return False, f"compile error: {exc}"


def validate_import_and_structure(name):
    """Import *name* and verify it exposes a well-formed XPARAMETERS dict."""
    try:
        mod = importlib.import_module(
            f"openwave.xperiments.m4_ewt.xparameters.{name}"
        )
        xp = mod.XPARAMETERS
    except Exception as exc:
        return False, f"import failed: {exc}"

    if not isinstance(xp, dict):
        return False, "XPARAMETERS is not a dict"
    if "engine" not in xp:
        return False, "missing 'engine' key"

    return True, "OK"


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

def main():
    modules = discover_xparameter_modules()
    if not modules:
        print("No xparameter modules found.")
        sys.exit(1)

    passed = 0
    for stem, path in modules:
        # 1. Compile check
        ok_compile, msg_compile = validate_compile(path)

        # 2. Import + structure check (only if compile passed)
        if ok_compile:
            ok_import, msg_import = validate_import_and_structure(stem)
        else:
            ok_import, msg_import = False, "skipped (compile failed)"

        ok = ok_compile and ok_import
        status = "PASS" if ok else "FAIL"
        print(f"  {stem:60s} {status}  compile: {msg_compile:20s}  import: {msg_import}")
        if ok:
            passed += 1

    print(f"\n{passed}/{len(modules)} xparameter modules passed")
    sys.exit(0 if passed == len(modules) else 1)


if __name__ == "__main__":
    main()