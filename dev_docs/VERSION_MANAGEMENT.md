# Version Management

## Overview

OpenWave uses a **single source of truth** approach for version management that works correctly with both regular installations and editable (development) installations.

## How It Works

### Version Definition

The version is defined in two places:

1. **Source code**: `openwave/__init__.py` contains `__version__ = "X.Y.Z"`
2. **Build metadata**: `pyproject.toml` contains `version = "X.Y.Z"`

Both files should always have the **same version number**.

### Version Access

The codebase accesses the version from the source code directly:

```python
from openwave import __version__
```

This approach ensures:

- Developers with editable installs (`pip install -e .`) see the current version immediately
- No need to reinstall after version bumps
- Works correctly with both development and production installations
- Fallback to metadata if `__version__` is not available

## For Developers

### Updating the Version

When bumping the version, update **both** files:

1. Edit `openwave/__init__.py`:

   ```python
   __version__ = "26.12"  # The previous release plus one
   ```

1. Edit `pyproject.toml`:

   ```toml
   version = "26.12"  # Update this to match
   ```

### Why This Works with Editable Installs

- **Editable install** (`pip install -e .`): Python imports directly from your source directory, so changes to `__init__.py` are immediately visible
- **Regular install** (`pip install .`): The version in `__init__.py` is copied during installation, and the metadata is generated from `pyproject.toml`

Both installation methods will see the same version number.

## Version Numbering

OpenWave uses **calendar versioning (CalVer)** in the form `YY.N`: the year the release was cut,
then the release count in that year. `26.11` is the eleventh-numbered release of 2026, and the
first release of 2027 is `27.1`: the count restarts each January. This is the scheme pip uses
(`24.0`, `24.1`, ..., `25.0`). Adopted 2026-10-05, replacing `YY.M.D` (the release date), which
had replaced SemVer at `1.6.10` on 2026-09-07.

### Why not SemVer

A SemVer number is a backward-compatibility promise addressed to a dependency resolver. This
project has no resolver to make one to: it is not published to PyPI, nothing pins a version
range against it, and with no CI there is nothing that could check such a promise even in
principle. `MINOR` versus `PATCH` was a judgment call made on every release for no reader.

### Why `YY.N` replaced `YY.M.D`

The date form needed a fourth component (`26.9.7.1`) for a second release on the same day, and a
letter there (`26.9.7b`) would have been read by PEP 440 as a beta of a release that already
shipped. A count has no same-day case: every release is the previous one plus one.

The tradeoff, stated plainly: the number no longer carries the day. The year still tells a reader
of a dated finding, method note or run record roughly where the engine stood, and the exact date
is the date of the `v<version>` git tag and its GitHub release. CalVer still says nothing about
compatibility, which was never verified anyway.

### Rules

| Rule | Reason |
| --- | --- |
| Each release is the previous one plus one: `26.11`, `26.12`, `26.13` | No judgment call, and no same-day collision |
| The first release of a new year is `YY.1`: `27.1` | The year in front would mean nothing if the count ran on forever (`27.143` would look like 143 releases in 2027) |
| No zero padding: `26.11`, never `26.011` | PEP 440 strips leading zeros, so the padded form installs as `26.11` and stops matching its own tag |
| No `v` in `__init__.py` or `pyproject.toml` | PEP 440 strips it too. The `v` prefix belongs on the git tag alone, which is where the existing convention already puts it |
| No letter suffix | ⚠️ PEP 440 reads `26.11b` as a **beta** of `26.11`, which sorts BEFORE it and which `pip install` skips by default |
| The count is releases (tags), not commits | Work between releases keeps the last released version |

### Ordering

`YY.N` sorts correctly under PEP 440 because each component is compared as an integer, not as
text: `26.9 < 26.11 < 26.100 < 27.1`. Each switch moved forward, never back: `1.6.10 < 26.9.7`,
and the first `YY.N`, `26.11`, is greater than every `26.10.x` and than the last dated release,
`26.9.22` (11 > 9 and 11 > 10 in the second component). Returning to either earlier scheme would
require a version decrease, which no installer would select.

## Implementation Details

### Files Using Version

The following files display the version to users:

- `openwave/i_o/cli.py`: CLI menu headers (lines 159-164, 223-229)
- `openwave/i_o/render.py`: Window title (lines 20-29)

All files use the same pattern:

```python
try:
    from openwave import __version__
    pkg_version = __version__
except ImportError:
    # Fallback to metadata if __version__ not available
    from importlib.metadata import version
    pkg_version = version("OPENWAVE")
```

### Why Not Use `importlib.metadata.version()` Only?

The `importlib.metadata.version()` function reads from package metadata installed by pip. This metadata is only updated when you reinstall the package:

- With editable installs, metadata is created once during `pip install -e .`
- Subsequent code changes (including version bumps) don't update the metadata
- Developers would need to run `pip install -e .` after every version bump

By reading from `__version__` in the source code, we avoid this issue entirely.

## Alternative Approaches

### setuptools-scm (Not Used)

An alternative approach is to use `setuptools-scm` to derive versions from git tags. We chose not to use this because:

- Adds complexity and dependencies
- Requires proper git tagging discipline
- Can be confusing when working with uncommitted changes
- Simple dual-file approach is more transparent and explicit

## When to Bump Version

### Recommended Workflow: Bump BEFORE Creating Tag/Release

The best practice is to bump the version **before** creating git tags and GitHub releases:

1. **Update the version** in your source code (`__init__.py` and `pyproject.toml`)
1. **Commit the version bump** with a clear message
1. **Create a git tag** matching that version
1. **Create a GitHub release** from that tag
1. **Publish to PyPI** (if applicable) using that tagged version

#### Why This Order?

This ensures:

- The tag points to code that actually contains that version number
- Users installing from that tag get the correct version
- Clear history: `git log` shows when each version was created
- The git tag and the package version are synchronized

#### Typical Workflow Example

```bash
# 1. Make your changes and test them
git add .
git commit -m "Add new feature X"

# 2. Bump the version: the previous release plus one, or YY.1 for the first release of a year
git tag --list 'v*' --sort=-v:refname | head -1   # the previous release, e.g. v26.11
# Edit: openwave/__init__.py → __version__ = "26.12"
# Edit: pyproject.toml → version = "26.12"

# 3. Commit the version bump
git add openwave/__init__.py pyproject.toml
git commit -m "Bump version to 26.12"

# 4. Create a git tag
git tag -a v26.12 -m "Release version 26.12"

# 5. Push everything
git push origin main
git push origin v26.12

# 6. Create GitHub release (via UI or gh cli)
gh release create v26.12 --title "v26.12" --notes "Release notes here"

# 7. (Optional) Publish to PyPI
python -m build
python -m twine upload dist/*
```

#### Alternative: Separate Version Bump Commit

Some teams prefer a dedicated "version bump" commit at the end of a release cycle:

```bash
# After all feature work is done:
git commit -m "Implement feature X"
git commit -m "Fix bug Y"
git commit -m "Update docs"

# Then bump version as last commit before tag
# Edit version files to the next number...
git commit -m "Bump version to 26.12"
git tag -a v26.12 -m "Release v26.12"
```

### What NOT to Do

- **Don't bump version AFTER creating the tag**: The tag would point to old version number
- **Don't commit version bumps on every commit**: Creates noise in git history
- **Don't leave version bumps uncommitted**: Other developers won't see the new version

### Which Number to Use

The previous release plus one, read from the tags:

```bash
git tag --list 'v*' --sort=-v:refname | head -1   # v26.11, so the next is 26.12
```

The only other case is the first release of a year: `YY.1` (`27.1`), whatever the last count of
the year before was.

### Pre-release Versions

⚠️ SemVer-style suffixes such as `0.2.0-dev` or `0.2.0-alpha.1` are NOT canonical PEP 440 and
silently mutate on install. Python rewrites them:

| Written | What pip records |
| --- | --- |
| `26.12-dev` | `26.12.dev0` |
| `26.12-alpha.1` | `26.12a1` |
| `26.12-beta.1` | `26.12b1` |
| `26.12-rc.1` | `26.12rc1` |

Write the canonical form directly if a pre-release is ever needed:

```python
__version__ = "26.12.dev0"   # Development (ongoing work)
__version__ = "26.12a1"      # Alpha (early testing)
__version__ = "26.12b1"      # Beta (feature complete)
__version__ = "26.12rc1"     # Release candidate (final testing)
__version__ = "26.12"        # Final release
```

⚠️ All four sort BEFORE `26.12`, and `pip install` skips them unless asked for with
`--pre`. That is what they are for. It is also why a hotfix must never be numbered `26.12b`:
it would be treated as a beta of a release that already shipped. A hotfix is simply the next number.

In practice a counted scheme has little use for these. The working version between releases
is simply the last released version until the day a new one is cut.

### Automation Options

Version bumping needs no tool now that the number is a count, and `bumpver update --patch`
style commands no longer map onto anything. If any automation is added later, the useful
target is checking that the two files agree and that the tag matches, not computing the number:

```bash
# Verify the two version strings are in sync before tagging
python -c "import tomllib,re,sys; \
t=tomllib.load(open('pyproject.toml','rb'))['project']['version']; \
i=re.search(r'__version__ = \"([^\"]+)\"',open('openwave/__init__.py').read()).group(1); \
sys.exit(0 if t==i else f'version mismatch: pyproject {t} vs __init__ {i}')"
```

## Best Practices Summary

1. Always update both `__init__.py` and `pyproject.toml` together
1. Bump version BEFORE creating git tags and releases
1. Use the pattern shown above when accessing version in code
1. Keep version numbers synchronized between source and build config
1. Use `YY.N`, the previous release plus one (`YY.1` in a new year), with no zero padding and no `v` prefix in the files
1. Create dedicated version bump commits
1. Document version changes in commit messages and release notes
