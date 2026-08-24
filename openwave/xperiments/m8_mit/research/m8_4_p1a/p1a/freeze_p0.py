"""P1A.0: Freeze P0.

The qualified p0/ package and its manifest are immutable inputs.
Verify the manifest hashes on every run. Any edit to qualified P0
code reopens P0.
"""

import hashlib
import os


P0_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "p0")

FROZEN_MANIFEST = {
    "__init__.py":         "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
    "algebra.py":          "2a27f19a8eb70923",
    "bundle_operator.py":  "5513f7092c448518",
    "cloud.py":            "0926c44ab0c6df92",
    "frozen_tolerances.py":"3f78c5711f1ac76b",
    "group.py":            "1409d25cbdcdbcfd",
    "qualify.py":          "af944887c1d6164c",
    "rbffd.py":            "40ca9d68b0d65f07",
    "regression_tests.py": "f772fbe005fb5d2f",
    "representations.py":  "d1439829d79d5e9b",
}


def compute_manifest(directory):
    manifest = {}
    for fname in sorted(os.listdir(directory)):
        if not fname.endswith('.py'):
            continue
        path = os.path.join(directory, fname)
        with open(path, 'rb') as f:
            data = f.read()
        manifest[fname] = hashlib.sha256(data).hexdigest()
    return manifest


def verify_p0_frozen():
    """Verify every P0 file matches the frozen manifest from the qualification note.

    Returns (pass, details).
    """
    current = compute_manifest(P0_DIR)
    details = {}
    all_pass = True

    expected_files = set(FROZEN_MANIFEST.keys())
    actual_files = set(current.keys())

    if expected_files != actual_files:
        details["file_set_mismatch"] = {
            "missing": sorted(expected_files - actual_files),
            "extra": sorted(actual_files - expected_files),
        }
        all_pass = False

    for fname in sorted(expected_files & actual_files):
        expected_prefix = FROZEN_MANIFEST[fname]
        actual_hash = current[fname]
        match = actual_hash.startswith(expected_prefix)
        details[fname] = {
            "expected_prefix": expected_prefix,
            "actual": actual_hash[:16] + "...",
            "match": match,
        }
        if not match:
            all_pass = False

    return all_pass, details
