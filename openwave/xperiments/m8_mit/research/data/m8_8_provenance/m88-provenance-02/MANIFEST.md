# Manifest

**20 files in this archive. Only the SIX below define `source_content_sha256`.**
Everything else travels with them so the derivation can be rerun and audited, and none of it
is covered by that hash. The distinction is the point of this file.

Recompute rather than read. This manifest is author-supplied like everything else here.

**`m8_8_reproduction_protocol.md` here is a SNAPSHOT, not a pin.** It records the protocol as
of delivery. The protocol is still an open PR and its § 11 pins land later by design, so
further edits are expected and do NOT invalidate this archive or require a rebuild. If you
need the live text, read the PR. Every other file here is either frozen (§ 1) or an input to
the derivation, and those do not drift.

| | |
| --- | --- |
| provenance ID | `M88-CONSTR-02` |
| construction packet SHA-256 | `df00c0222f98c481eb56b882cd867a6c3a4f8604b8633e81dec0cce1f8460a06` |
| derivation source SHA-256 | `8a3a1c87f54372a446356a5c2a5ece4d9b4ba7a32367ef129b8baf18b44733f6` |
| M8.5-A group packet SHA-256 | `e3b0c945bbbb15b4549fa641234c9461062c2337b3d1e372af621b614d4883a9` |
| frozen prime | `p = 10**9 + 7` |

## 1. The frozen six, and ONLY these, define `source_content_sha256`

That statement is about the HASH, and it should not be read as a derivation closure: the
derivation also consumes two inputs outside the hashed code, the M8.5-A group packet and
`pres.json` (the generator choice). Both are closed elsewhere rather than left loose: the
group packet by `group_packet_sha256` inside the construction packet, and `pres.json` by the
packet's own `abstract_generators`, which echo it and which the maintainer audit validates
against the independently rebuilt closure (A3). Vary either input and the packet itself
changes.

SHA-256 over the concatenation of exactly these files, **in this order**, raw bytes, nothing
between them. Editing any of them, including whitespace or a comment, invalidates the
published hash.

```
cd build && cat qphi.py complex.py kernel.py sat.py search_sym.py build_packet_v2.py | shasum -a 256
```

| # | path | sha256 | bytes |
| --- | --- | --- | --- |
| 1 | `build/qphi.py` | `6819952e076b409183c451efa6397d5d7628ab958bb6a5dd9f3ad17e1bfd8775` | 2124 |
| 2 | `build/complex.py` | `344fbea1d6e07005825428b05e0e7db1f84420ab3278b7f1acf3685ec39979f9` | 1691 |
| 3 | `build/kernel.py` | `966c799e8a7825fdf89ef6ba62eb49aff14b0c43b8bc0eff6ab4d0790feac26e` | 1397 |
| 4 | `build/sat.py` | `d057815902f027f6de3803edeac084d62135aabe4452f133048df1e2e37259e2` | 3500 |
| 5 | `build/search_sym.py` | `7fc8bf4db0ae6fd1b8bb37624da2d190cbf13cb58ec25c244fd74ff50654a85e` | 2848 |
| 6 | `build/build_packet_v2.py` | `0698958be6a200a8f5189ccc07b2b073a4f67b93da3ebddcb453f98bfb512fc6` | 2306 |

## 2. Execution and recovery, NOT covered by the source hash

Needed to rerun and verify. `certify.py` in particular must ship: it is the certificate
generator, and the only file whose int64 arithmetic can approach overflow, so it carries the
guard for its own elimination. (The frozen `search_sym.py` has an int64 screen too, but at
frozen primes of at most 37 it cannot overflow; frozen `sat.py` is pure-Python arbitrary
precision and cannot overflow at any prime.)

| path | sha256 | bytes |
| --- | --- | --- |
| `build/certify.py` | `72950ee24c258a289848e0b5d4b0d228937f6f77726aa2a18d9fb4ab74984c2e` | 8421 |
| `build/m3_v2.json` | `74329e7b55cbd55fb288f0738dead159f2fe691aaabb65b827513aa69b40f48e` | 239 |
| `build/pres.json` | `224ee4b9ef3a1f20aeaadf44eb0b40548a5ed80cce88ebb567742fa009e38473` | 25 |
| `m8_5a_packet.json` | `e3b0c945bbbb15b4549fa641234c9461062c2337b3d1e372af621b614d4883a9` | 565 |
| `m8_8_construction_packet.json` | `df00c0222f98c481eb56b882cd867a6c3a4f8604b8633e81dec0cce1f8460a06` | 2657 |
| `recover.py` | `4c6e0c38699bad03cfcdafdcf69104f173cb92a44977b152402015f1b6bb2160` | 6755 |

## 3. Contextual audit material, NOT covered by the source hash

| path | sha256 | bytes |
| --- | --- | --- |
| `ENVIRONMENT.md` | `1653b1ec3d3cbca3e18a091ab40798c181c5dfc681a87240951270f1a53ab01a` | 9413 |
| `README.md` | `ab0ddee5c708af97bc64677fb8875e72a6ae6e356ead957d648997fb407e7dd6` | 4366 |
| `build/audit_A1_A11_output.txt` | `94b3e421d2e290a4fcb113484c5a2d2470f3ed1c085b3c2da28802c0f43d005e` | 3503 |
| `build/author_side/construction_audit.md` | `69a8ddd17e949d7ff5925e3216bab0251d19e68826e3290c944303d57beb6fd9` | 16355 |
| `build/m8_8_packet_audit.py` | `4d33adc79a00e27a3ba318de84695adeaa143c74acfc81a84609bef8fb708664` | 46985 |
| `build/saturation_certificate.json` | `4f1b43991e49d1c5cfc07ca3e5519c92cfb15a0ec7c804356fd66c493cb56e66` | 2184 |
| `m8_8_reproduction_protocol.md` | `5543ece4ea7483e2670a8d7be68ce3c1b631b7fcd699a9a6d57295213a19c2da` | 50268 |

## Candidate ordering, and why the seed does not matter

`search_sym.py` initializes an RNG, but the ordered candidate list begins with the 119
saturated kernel-basis vectors and accepts **basis vector 6**, candidate 7 of 6119, inside
that deterministic prefix, so the header seed feeds only candidates that are never reached;
the minor-sampling RNG inside the acceptance test is seeded per candidate by index,
independent of the header seed. Rerunning the file exactly as frozen reproduces the
canonical packet byte for byte. See `ENVIRONMENT.md` for the full statement; do not
paraphrase this as "seed-free".

## Container command

```
docker run --rm -e PYTHONDONTWRITEBYTECODE=1 -v "$PWD:/w" -w /w python:3.13-slim sh -c \
  "pip install --quiet numpy==2.5.0 && python recover.py"
```
