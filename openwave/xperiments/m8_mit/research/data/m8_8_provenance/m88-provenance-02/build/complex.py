"""Fox calculus over Z[2I], and the based 3-complex of S^3/2I."""
import sys; sys.path.insert(0,'.')
from qphi import *
from collections import defaultdict

P='m8_5a_packet.json'
GENS,_ = parse_packet(P)
G = close_group(GENS)
ORDER, IDS = canonical_ids(G)
KEY = lambda e: tuple((c[0],c[1]) for c in e)
N = len(ORDER)
def gid(x): return IDS[KEY(x)]
MUL = [[gid(hmul(ORDER[a], ORDER[b])) for b in range(N)] for a in range(N)]
INV = [next(b for b in range(N) if MUL[a][b] == gid(HONE)) for a in range(N)]
E = gid(HONE)

# --- group-ring elements: dict id -> int coefficient ---
def rz(*terms):
    d = defaultdict(int)
    for c, g in terms: d[g] += c
    return {g: c for g, c in d.items() if c}
def radd(x, y):
    d = defaultdict(int)
    for g, c in x.items(): d[g] += c
    for g, c in y.items(): d[g] += c
    return {g: c for g, c in d.items() if c}
def rneg(x): return {g: -c for g, c in x.items()}
def rmul(x, y):
    d = defaultdict(int)
    for g1, c1 in x.items():
        for g2, c2 in y.items(): d[MUL[g1][g2]] += c1*c2
    return {g: c for g, c in d.items() if c}
def rconj(x): return {INV[g]: c for g, c in x.items()}
ZR, UNIT = {}, {E: 1}

# --- Fox calculus on words, letters as (gen_index, +/-1) ---
def word_eval(w, gens):
    x = E
    for i, e in w: x = MUL[x][gens[i] if e > 0 else INV[gens[i]]]
    return x
def fox(w, k, gens):
    """d(w)/d(g_k) in Z[G], w a list of (gen_index, exponent sign)."""
    out, pref = ZR, E
    for i, e in w:
        if e > 0:
            if i == k: out = radd(out, {pref: 1})
            pref = MUL[pref][gens[i]]
        else:
            pref = MUL[pref][INV[gens[i]]]
            if i == k: out = radd(out, {pref: -1})
    return out
