"""Route (a) group-closure battery (Addendum 12.3 requalification).

The coverage this battery adds is exactly what the M85B-ADJ-04 STOP proved
missing: route (a)'s stencil path exercised on an effective group that is
NOT sign-exact closed. Every check prints PASS or FAIL and the process
exits nonzero on any FAIL. The mutation arm rebuilds the PRE-repair
sign-exact lookup locally and requires it to fail, so the repair is
load-bearing, and the truncation arm requires refusal, so lookup did not
become nearest-match.
"""
import sys

for _p in ("production", "pilot", "gates"):
    sys.path.insert(0, _p)

import numpy as np                                          # noqa: E402
from route_a_twosided import ONE, close_pairs, cloud, act_matrix, qmul  # noqa: E402
from equivariant_stencils import group_multiplication_table  # noqa: E402
import route_a_producer as route_a                          # noqa: E402
import route_b_producer as route_b                          # noqa: E402
import step3_runner                                          # noqa: E402

fails = 0
def ck(label, ok, detail=""):
    global fails
    print(f"  {'PASS' if ok else 'FAIL'}  {label}" + (f"   {detail}" if detail else ""))
    if not ok:
        fails += 1

mk = lambda t: np.array([np.cos(t), np.sin(t), 0.0, 0.0])
def lens_gens(q, s1, s2):
    return [(mk(np.pi * (s1 + s2) / q), mk(np.pi * (s1 - s2) / q))]

print("route (a) group-closure battery")

# 1. homogeneous regression
p31 = close_pairs(lens_gens(3, 1, 1))
try:
    group_multiplication_table(p31)
    ck("homogeneous L(3,1): table closes (regression)", len(p31) == 3)
except ValueError as e:
    ck("homogeneous L(3,1): table closes (regression)", False, str(e))

# 2. binary polyhedral regression: 2T, which contains -1 natively
h = 0.5
tet = [(np.array([h, h, h, h]), ONE.copy()), (np.array([0.0, 1.0, 0.0, 0.0]), ONE.copy())]
p2t = close_pairs(tet)
try:
    group_multiplication_table(p2t)
    ck("binary polyhedral 2T: table closes (regression)", len(p2t) == 24,
       f"order {len(p2t)}")
except ValueError as e:
    ck("binary polyhedral 2T: table closes (regression)", False, str(e))

# 3+4. the previously failing regime, with ground truth from the ACTION,
# which is sign-free: T[i][j] must index the element whose SO(4) matrix
# equals the matrix product. This is correctness, not mere non-crashing.
p712 = close_pairs(lens_gens(7, 1, 2))
try:
    T = group_multiplication_table(p712)
except ValueError as e:
    # the pre-repair state reaches here: report it as a scored FAIL rather
    # than dying, so a reverted repair produces a readable battery, not a
    # traceback. The exit code is nonzero either way.
    T = None
    ck("inhomogeneous L(7;1,2): every table entry matches the SO(4) ground truth",
       False, f"table could not be built: {e}")
    ck("central-partner products actually occur and resolve", False, "table absent")
if T is not None:
    mats = [act_matrix(u, v) for u, v in p712]
    good = all(np.abs(mats[i] @ mats[j] - mats[T[i][j]]).max() < 1e-8
               for i in range(len(p712)) for j in range(len(p712)))
    wraps = sum(1 for i in range(len(p712)) for j in range(len(p712))
                if not (np.abs(qmul(p712[i][0], p712[j][0]) - p712[T[i][j]][0]).max() < 1e-9))
    ck("inhomogeneous L(7;1,2): every table entry matches the SO(4) ground truth", good)
    ck("central-partner products actually occur and resolve", wraps > 0,
       f"{wraps} of {len(p712)**2} products land on the central partner")

# 5. mutation: the PRE-repair sign-exact lookup must FAIL here
def sign_exact_table(pairs, tol=1e-9):
    def same(a, b):
        return abs(a[0] - b[0]).max() < tol and abs(a[1] - b[1]).max() < tol
    def idx_of(p):
        for i, q in enumerate(pairs):
            if same(p, q):
                return i
        raise ValueError("product left the supplied element list")
    return [[idx_of((qmul(pairs[i][0], pairs[j][0]), qmul(pairs[i][1], pairs[j][1])))
             for j in range(len(pairs))] for i in range(len(pairs))]
try:
    sign_exact_table(p712)
    ck("mutation: removing central equivalence makes L(7;1,2) fail", False,
       "sign-exact lookup unexpectedly succeeded")
except ValueError:
    ck("mutation: removing central equivalence makes L(7;1,2) fail", True)

# 6. genuinely absent product still refused
try:
    group_multiplication_table(p712[:-1])
    ck("truncated element list is refused (lookup is not nearest-match)", False)
except ValueError:
    ck("truncated element list is refused (lookup is not nearest-match)", True)

# 7. cloud cardinality: effective order exactly, no coincident points
rng = np.random.default_rng(20260816)
S = rng.normal(size=(12, 4)); S /= np.linalg.norm(S, axis=1, keepdims=True)
X, oid, gid, M = cloud(S, p712)
from scipy.spatial.distance import pdist
ck("cloud cardinality is seeds x effective order, with no coincident points",
   len(X) == 12 * len(p712) and pdist(X).min() > 1e-6,
   f"{len(X)} nodes, min separation {pdist(X).min():.2e}")

# 8. route (a) END TO END on the previously failing regime, frozen stencil
# configuration at a reduced ladder rung, validated route-locally by the
# frozen step-3 runner (which includes the full-band census)
try:
    art_a = route_a.produce(p712, "requal-a-L712", "m8_5b-v1-requal",
                            "SYN-L712-REQUAL", seeds=60, adjudication=False)
    probs = step3_runner._validate_one(art_a, p712, "route a")
    ck("route (a) traverses the stencil path end to end on L(7;1,2)",
       not probs, "; ".join(probs[:2]))
except Exception as e:
    # the pre-repair state cannot reach the end of this path; report it as a
    # scored FAIL so the whole battery stays readable under a reverted repair
    ck("route (a) traverses the stencil path end to end on L(7;1,2)", False,
       f"{type(e).__name__}: {e}")

# 9a. route (b) is UNCHANGED by this repair, proven at RUNTIME: a fresh
# subprocess executes the same produce() path 9b tests, then asks whether
# the repaired module ever entered sys.modules. Import-time closure alone
# would not establish this, since a lazy import inside produce() would
# evade it; running the path first makes the sentence true of what was
# actually executed.
import subprocess as _sp
_probe_code = (
    "import sys\n"
    "sys.path[:0] = ['production', 'pilot', 'gates']\n"
    "import numpy as np\n"
    "import route_b_producer\n"
    "mk = lambda t: np.array([np.cos(t), np.sin(t), 0.0, 0.0])\n"
    "gens = [(mk(np.pi * 3 / 7), mk(-np.pi / 7))]\n"
    "route_b_producer.produce(gens, 'requal-b-9a', 'm8_5b-v1-requal',\n"
    "                         'SYN-L712-REQUAL', adjudication=False)\n"
    "print('equivariant_stencils' in sys.modules)\n"
)
probe = _sp.run([sys.executable, "-c", _probe_code],
                capture_output=True, text=True, cwd=".")
ck("repaired module never enters route (b)'s runtime through produce() "
   "(runtime dependency closure)",
   probe.stdout.strip().endswith("False"),
   probe.stdout.strip()[-24:] or probe.stderr[:80])

# 9b. and its scalar records equal an INDEPENDENT character-sum
# recomputation performed here with this battery's own arithmetic, so the
# baseline is mathematics, not a rerun of the same code.
art_b = route_b.produce(lens_gens(7, 1, 2), "requal-b-L712", "m8_5b-v1-requal",
                        "SYN-L712-REQUAL", adjudication=False)
p712b = close_pairs(lens_gens(7, 1, 2))
def chi(quat, n):
    ang = np.arctan2(np.linalg.norm(quat[1:]), quat[0])
    if abs(np.sin(ang)) < 1e-12:
        return (n + 1) * (np.cos(ang) ** n if n else 1.0)
    return np.sin((n + 1) * ang) / np.sin(ang)
sc = [r for r in art_b["records"] if r["sector_id"] == "scalar"
      and r["quotient_multiplicity"] is not None]
ok_b = all(abs(r["quotient_multiplicity"]
               - sum(chi(u, r["harmonic_level"]) * chi(v, r["harmonic_level"])
                     for u, v in p712b) / len(p712b)) < 1e-9 for r in sc)
ck("route (b) scalar records equal the independent character-sum baseline",
   ok_b and len(sc) > 0, f"{len(sc)} scalar cells recomputed independently")

print(f"\n{('BATTERY PASS' if fails == 0 else 'BATTERY FAIL')}: {fails} failures")
sys.exit(1 if fails else 0)
