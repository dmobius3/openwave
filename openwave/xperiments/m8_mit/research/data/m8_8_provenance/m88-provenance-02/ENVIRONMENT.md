# Environment record and rerun instructions

Delivered out of band per § 4.2 check 1 and check 2, so these bytes stay clear of the
clean-room base the room opens from.

| | |
| --- | --- |
| provenance ID | `M88-CONSTR-02` |
| provenance class | `derived` (§ 4.2) |
| construction packet SHA-256 | `df00c0222f98c481eb56b882cd867a6c3a4f8604b8633e81dec0cce1f8460a06` |
| derivation source SHA-256 | `8a3a1c87f54372a446356a5c2a5ece4d9b4ba7a32367ef129b8baf18b44733f6` |
| M8.5-A group packet SHA-256 | `e3b0c945bbbb15b4549fa641234c9461062c2337b3d1e372af621b614d4883a9` |
| certificate prime | `p = 10**9 + 7`, the one prime every rank in the published certificate used |
| built under | Python 3.13.13, NumPy 2.5.0, macOS (darwin) |

Compare hashes in FULL. The group-packet hash begins `e3b0c9…`, which eyeball-collides with
the SHA-256 of the empty string (`e3b0c442…`); a glance-comparison can confuse the two in
either direction, including the case where a pipeline hashes a missing file and gets the
empty-string value.

Every command below was executed as written against this archive before being documented.
Paths in the frozen sources are cwd-relative, so each recipe states its working directory.
The rebuild steps fail with `FileNotFoundError` from any other directory, which is a
property of the frozen bytes, not damage; `certify.py` happens to also run from the archive
root, where both packets sit.

## The `source_content_sha256` concatenation rule, exactly

SHA-256 over the concatenation of these six files, **in this order**, raw bytes, nothing
between them:

1. `build/qphi.py`
2. `build/complex.py`
3. `build/kernel.py`
4. `build/sat.py`
5. `build/search_sym.py`
6. `build/build_packet_v2.py`

From the archive root:

```
cat build/qphi.py build/complex.py build/kernel.py build/sat.py build/search_sym.py build/build_packet_v2.py | shasum -a 256
```

Editing any of the six, including whitespace or a comment, invalidates the published hash.

## Independent verification, from the archive root

`recover.py` is an author-supplied self-check and CANNOT serve as the audit: it runs inside
the artifact it verifies, and its `RECOVERY VERIFIED` line is informational only. Verify
independently instead; all three steps run from the archive root:

```
# 1. every manifest-listed file: hash and byte count
python3 - <<'EOF'
import hashlib, re, pathlib
man = pathlib.Path('MANIFEST.md').read_text(); bad = n = 0
for m in re.finditer(r'\|(?:\s*\d+\s*\|)?\s*`([^`]+)` \| `([0-9a-f]{64})` \| (\d+) \|', man):
    p = pathlib.Path(m.group(1)); n += 1
    ok = hashlib.sha256(p.read_bytes()).hexdigest() == m.group(2) and p.stat().st_size == int(m.group(3))
    bad += not ok
    if not ok: print('MISMATCH', m.group(1))
print(f'{n - bad} of {n} match')
EOF

# 2. the frozen six, which alone define source_content_sha256
cat build/qphi.py build/complex.py build/kernel.py build/sat.py build/search_sym.py build/build_packet_v2.py | shasum -a 256

# 3. the packet itself
shasum -a 256 m8_8_construction_packet.json
```

Expected: `19 of 19 match` for step 1 (`MANIFEST.md` lists every file but itself); then
`8a3a1c87f54372a446356a5c2a5ece4d9b4ba7a32367ef129b8baf18b44733f6`; then
`df00c0222f98c481eb56b882cd867a6c3a4f8604b8633e81dec0cce1f8460a06`. If any recomputation
disagrees with any file in this archive, believe the recomputation.

## Rerunning the maintainer audit, from the archive root

```
python3 build/m8_8_packet_audit.py --packet m8_8_construction_packet.json --group-packet m8_5a_packet.json --no-write --mutation-tests
```

`--no-write` matters: the script is the maintainer's, shipped byte-identical to their `main`
so that a rerun is a rerun of their gate rather than of a fork of it, and without the flag it
tries to write its result into a `data/` directory this archive does not have and dies on a
`FileNotFoundError` after the checks have already passed. The transcript of exactly this
invocation is `build/audit_A1_A11_output.txt`.

## Rebuilding the packet from scratch

Work in a COPY of the archive, never the original tree: the rebuild rewrites `m3_v2.json`
and the packet in its working directory (deterministically and byte-identically, but a
pristine tree is the point of having a manifest). The frozen sources read both packets from
their working directory, so copy them in first:

```
cp -R . /tmp/m88-work && cd /tmp/m88-work
cp m8_5a_packet.json m8_8_construction_packet.json build/ && cd build
mv m3_v2.json m3_v2.json.shipped
python3 search_sym.py | tee search.log
grep -q "FOUND: basis 6" search.log
cmp m3_v2.json m3_v2.json.shipped
python3 build_packet_v2.py
cmp m8_8_construction_packet.json ../m8_8_construction_packet.json && echo BYTE-IDENTICAL
python3 certify.py        # prints the verdict; add --emit to regenerate the JSON
```

The `mv` and the `grep` are load-bearing, not ceremony. `search_sym.py` exits 0 even when it
accepts nothing, and the shipped `m3_v2.json` would otherwise still be sitting in `build/`,
so a silently failed search would rebuild the packet from the stale file and the final `cmp`
would print `BYTE-IDENTICAL` anyway. Setting the shipped file aside forces the search to
reproduce it, the `grep` fails if the search accepted nothing, and the first `cmp` fails if
it accepted anything else. Both failure paths were exercised against a sabotaged search
before this recipe was documented.

`certify.py` writes nothing without `--emit`, precisely so a verification rerun cannot
perturb a manifest-listed file.

**Set `PYTHONDONTWRITEBYTECODE=1`, or expect `__pycache__`.** Python writes bytecode beside
any module it imports, so running anything creates `build/__pycache__/*.pyc`. Those are
Python's artifacts, not this archive's, and they are deliberately absent from `MANIFEST.md`.

## The container command

```
docker run --rm -e PYTHONDONTWRITEBYTECODE=1 -v "$PWD:/w" -w /w python:3.13-slim sh -c \
  "pip install --quiet numpy==2.5.0 && python recover.py"
```

`recover.py` copies everything to a temporary directory first and re-hashes the full
manifest on exit, so it cannot mutate the tree; it also runs the maintainer audit's mutation
suite by default (`--skip-mutations` to drop that step, at the cost of never demonstrating
the gates can fail).

## Reproducibility, stated exactly

The derivation is deterministic end to end, and rerunning `search_sym.py` exactly as frozen
reproduces the canonical packet byte for byte, which is the § 4.2 check 2 obligation.

The frozen file contains `random.Random(20260803)`. That header seed is inert for the
accepted result, for two reasons readable from the frozen code rather than established by
editing it: the 119 kernel-basis vectors are queued as candidates BEFORE the RNG generates
any sparse candidate, and acceptance happens at basis vector 6, candidate 7 of 6119, inside
that deterministic prefix, so the header seed feeds only candidates that are never reached;
and the minor-sampling RNG inside the acceptance test is seeded per candidate by its index,
`random.Random(tried)`, independent of the header seed. Do not paraphrase this as
"seed-free": the seed exists, it is frozen into the hashed bytes, and it is provably inert,
which is a different and stronger statement. One field to read correctly: `m3_v2.json`
records `"seed": 20260803` as a literal from the frozen source, so it documents the constant
in the file, not a measurement of the run that produced the output.

## Why `certify.py` travels with the archive

It is the certificate generator, and it is the only file in this derivation whose int64
arithmetic can approach overflow, so it carries the guard for its own elimination
(`(p-1)^2` must stay below `2^63`; past that, products wrap silently and the routine
returns a wrong rank rather than raising). Stated carefully, because an earlier version of
this sentence overclaimed in the other direction: the frozen `search_sym.py` also contains
an int64 modular screen, but its primes are frozen at 37 and below, where every entry is
reduced below the prime and products stay under 37², sixteen orders of magnitude inside the
bound, so it cannot overflow as frozen. The frozen `sat.py` is pure-Python
arbitrary-precision arithmetic and cannot overflow at any prime. The one place a
user-visible prime meets int64 is `certify.py`'s `rank_p` at `p = 10⁹ + 7`. An earlier version of this
file claimed the guard lived in `certify.py` because frozen `sat.py` could not absorb it;
that rationale was false and is corrected here. What remains true, and general: a hash-pinned
source set cannot take even a strictly safe fix without invalidating itself, so repairs land
in files outside the frozen set, and those files must then travel with the archive.

## What is in here and what is not

Included: the frozen six, `certify.py`, the recovery wrapper, the construction packet, the
M8.5-A group packet, the derivation inputs `pres.json` (the generator choice, echoed by the
packet's `abstract_generators`) and `m3_v2.json` (the accepted generator, reproduced by the
rebuild recipe), the maintainer audit script byte-identical to their `main`, its full A1-A11
transcript including the mutation suite, the saturation certificate for the shipped packet,
and the author-side construction-audit account.

Not included: any reproduced quantity. No torsion value, ratio, sector product or decimal
rendering of any reproduced quantity appears anywhere in this archive.
