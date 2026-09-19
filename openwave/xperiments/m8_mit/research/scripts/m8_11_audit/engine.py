"""Exact section calculus (the P-trick): functions on SU(2) with values in C^7 as sums of x^T D^J(g) Y."""
import pickle
from fractions import Fraction
from mq import MQ, ZERO, ONE, mq
from su2 import couple, theta, cg

with open("out/group.pkl", "rb") as f:
    GRP = pickle.load(f)
P = GRP["P"]  # d -> 7x7 exact projector


def Pcols(d):
    M = P[d]
    return [[M[k][a] for k in range(7)] for a in range(7)]


def block(u, d):
    """Phi~ = u^T D^3 P as {J: [(x, Ycols)]}."""
    return {3: [(list(u), Pcols(d))]}


def vdot(x, y):
    """x^dag y."""
    s = ZERO
    for a, b in zip(x, y):
        if not a.is_zero() and not b.is_zero():
            s = s + a.conj() * b
    return s


def vadd(x, y):
    return [a + b for a, b in zip(x, y)]


def vscale(x, s):
    return [a * s for a in x]


def trYY(Y1, Y2):
    return sum((vdot(c1, c2) for c1, c2 in zip(Y1, Y2)), ZERO)


def inner(f, h):
    """<f, h> = int sum_a conj(f_a) h_a."""
    s = ZERO
    for J in f:
        if J not in h:
            continue
        for x, Y in f[J]:
            for x2, Y2 in h[J]:
                a = vdot(x, x2)
                if a.is_zero():
                    continue
                b = trYY(Y, Y2)
                if b.is_zero():
                    continue
                s = s + a * b * Fraction(1, 2 * J + 1)
    return s


def fadd(*fs):
    out = {}
    for f in fs:
        for J, terms in f.items():
            out.setdefault(J, []).extend(terms)
    return out


def fscale(f, s):
    """multiply function by scalar s (applied to the left vector)."""
    return {J: [(vscale(x, s), Y) for x, Y in terms] for J, terms in f.items()}


def trilinear(f1, f2, f3, levels=None):
    """T(f1,f2,f3)_a = sum_b conj(f1_b) f2_b f3_a; restricted to output levels (set of J) if given."""
    out = {}
    for j1, t1 in f1.items():
        for j2, t2 in f2.items():
            for j3, t3 in f3.items():
                for x1, Y1 in t1:
                    tx1 = theta(x1, j1)
                    tY1 = [theta(c, j1) for c in Y1]
                    for x2, Y2 in t2:
                        for L in range(abs(j1 - j2), j1 + j2 + 1):
                            Js = [J for J in range(abs(L - j3), L + j3 + 1) if levels is None or J in levels]
                            if not Js:
                                continue
                            X12 = couple(tx1, j1, x2, j2, L)
                            if all(c.is_zero() for c in X12):
                                continue
                            W = [ZERO] * (2 * L + 1)
                            for b in range(7):
                                W = vadd(W, couple(tY1[b], j1, Y2[b], j2, L))
                            if all(c.is_zero() for c in W):
                                continue
                            for x3, Y3 in t3:
                                for J in range(abs(L - j3), L + j3 + 1):
                                    if levels is not None and J not in levels:
                                        continue
                                    A = couple(X12, L, x3, j3, J)
                                    if all(c.is_zero() for c in A):
                                        continue
                                    Z = [couple(W, L, Y3[a], j3, J) for a in range(7)]
                                    if all(c.is_zero() for col in Z for c in col):
                                        continue
                                    out.setdefault(J, []).append((A, Z))
    return out


def compress(f):
    """Rewrite each level as a sum over a linearly independent set of right factors (exact)."""
    out = {}
    for J, terms in f.items():
        basis = []   # list of (flatZ, Z, x_accum, pivot_index, reduced_vector)
        red = []     # reduced echelon rows: (pivot, row) in the coordinates of flattened Z
        coeffs = []  # for each basis element: accumulated x
        # represent each Z in terms of chosen basis Zs via exact elimination
        chosen = []  # list of (Z, flat)
        ech = []     # list of (pivot, row(list), combination over chosen (list of MQ))
        xs = []
        for x, Z in terms:
            flat = [c for col in Z for c in col]
            comb = [ZERO] * len(chosen)
            row = list(flat)
            for piv, erow, ecomb in ech:
                if not row[piv].is_zero():
                    fct = row[piv] / erow[piv]
                    row = [a - fct * b for a, b in zip(row, erow)]
                    comb = [a + fct * b for a, b in zip(comb + [ZERO] * (len(ecomb) - len(comb)), ecomb)]
            nz = [i for i, c in enumerate(row) if not c.is_zero()]
            if nz:
                # new independent direction: row = flat - sum comb_k chosen_k
                chosen.append((Z, flat))
                xs.append(list(x))
                newcomb = [-c for c in comb] + [ONE]
                ech.append((nz[0], row, newcomb))
                # ech rows expressed over chosen: row = sum newcomb_k chosen_k ; extend old combos
                ech = [(p, r, c + [ZERO] * (len(chosen) - len(c))) for p, r, c in ech]
            else:
                # flat = sum_k comb_k chosen_k
                for k, c in enumerate(comb):
                    if not c.is_zero():
                        xs[k] = vadd(xs[k], vscale(x, c))
        out[J] = [(xs[k], chosen[k][0]) for k in range(len(chosen))]
    return out


def block_fibre(f, d):
    """fibre vector of the level-3 part (right factors must be multiples of P)."""
    Pm = P[d]
    Pflat = [Pm[k][a] for a in range(7) for k in range(7)]
    tot = [ZERO] * 7
    for x, Y in f.get(3, []):
        flat = [c for col in Y for c in col]
        c = sum((Pm[k][a] * Y[a][k] for a in range(7) for k in range(7)), ZERO) * Fraction(1, d)  # tr(P Y)/d (P Hermitian: tr(P^dag Y))
        # c = tr(P^dag Y)/d with P^dag=P ; check exact proportionality
        for a1, a2 in zip(flat, Pflat):
            if not (a1 - c * a2).is_zero():
                raise AssertionError("level-3 right factor not proportional to P")
        tot = vadd(tot, vscale(x, c))
    return tot


def norm2(f):
    return inner(f, f)


def N(f, levels=None):
    return trilinear(f, f, f, levels)


def DN(f, h, levels=None):
    return fadd(trilinear(h, f, f, levels), trilinear(f, h, f, levels), trilinear(f, f, h, levels))
