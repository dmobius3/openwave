"""Standalone rung-3b theorem evaluator: p-form multiplicities of a lens space.

PROVENANCE.  This implements `THEOREM_TRANSCRIPTION.md`, which is the
implementation contract, NOT the papers directly.  The papers are provenance
sources and are retained under `../refs/` with their measured SHA-256:

    theorem      Lauret, arXiv:1604.02471v5, Thm 2.1 / Lem 3.2 / eq (3.1) / Thm 3.3
                 6b20b886450706d72467121560662fda6f88b319323930f3270ed31741fcb244
    generator    Lauret, Miatello, Rossetti, arXiv:1311.7167v4, eq (3.2)
                 db17fc71f629c23e3c2f674baadbd6da19ddf81cecf592090eb27f2f667ad702

Two corrections are applied, and they are DIFFERENT IN KIND.  Each must redden a
different check; a single passing tower comparison establishes at most one.

    CORRECTION A, group-action convention.  Lauret v5 section 2 prints the
    generator's final rotation block with s_1 repeated.  LMR v4 eq (3.2) carries
    the sequence through to s_n.  The LMR form is used.

    CORRECTION B, spectral-index correspondence.  Theorem 3.3 AS PRINTED is
    inconsistent with Theorem 2.1, with its own proof and eq (3.1), and with the
    paper's own 0-spectrum specialization.  The shifted correspondence is used:
        mult(lambda_{k,p-1}) = M_Gamma(k, p)
        mult(lambda_{k,p})   = M_Gamma(k, p+1)
    No claim is made about the typeset version of record, which was not fetched.

ISOLATION (protocol section 4.2 step 3).  This module is standalone.  It imports
nothing from route (a) or route (b), and no route may import it.  It is not a
shared library.  Its only case-specific inputs are the parameters (q, s) supplied
at evaluation time; nothing case-specific appears in this file.

ARITHMETIC.  Exact integers throughout.  No floating point anywhere in the
multiplicity path.
"""

from math import comb as _math_comb

__all__ = [
    "paper_binom",
    "generator_block_parameters",
    "in_congruence_lattice",
    "lattice_counts",
    "weight_multiplicity",
    "M_Gamma",
    "p_form_spectrum",
]

FORMULA_PROVENANCE = (
    "Lauret arXiv:1604.02471v5 Thm 2.1 / Lem 3.2 / eq (3.1) / Thm 3.3; "
    "generator convention from LMR arXiv:1311.7167v4 eq (3.2) in place of the "
    "repeated-final-parameter form printed in Lauret section 2; "
    "Theorem 3.3 branch mapping shifted per THEOREM_TRANSCRIPTION.md section 4B"
)


# --- layer 0: the paper's binomial convention --------------------------------

def paper_binom(b, a):
    """C(b, a) under the convention stated in Lauret v5, Notation 3.1:

        C(b, a) = 0   if a < 0 or b < a
                = the usual binomial coefficient otherwise.

    This is NOT math.comb, which raises ValueError on negative arguments.  The
    sums below generate negative upper arguments routinely, and the paper's
    convention is that those terms vanish.  Scattering try/except or explicit
    range guards through the finite sums would let Python's exception semantics
    stand in for the mathematics; one named primitive keeps the sums a literal
    transcription.
    """
    if a < 0 or b < a:
        return 0
    return _math_comb(b, a)


# --- layer 1: generator and action construction (CORRECTION A lives here) ----

def generator_block_parameters(s):
    """The parameter carried by each 2x2 rotation block of gamma, in order.

    gamma = diag( R(2*pi*s_1/q), ..., R(2*pi*s_n/q) )   [LMR v4 eq (3.2)]

    Block j carries s_j.  The LAST block carries s_n, not s_1.  Lauret v5
    section 2 prints s_1 in the final block; that printed form is CORRECTION A
    and is mutation 1.
    """
    return tuple(s)


def in_congruence_lattice(a, q, blocks):
    """mu = sum a_j eps_j lies in L_Gamma iff sum a_j s_j = 0 (mod q).

    Lauret v5 eq (2.5); LMR v4 eq (3.3).  Derived from gamma^mu = 1, so it
    consumes the BLOCK PARAMETERS rather than s directly: a wrong final block
    changes the lattice, which is how correction A becomes observable here.
    """
    return sum(aj * sj for aj, sj in zip(a, blocks)) % q == 0


# --- layer 2: congruence-lattice census --------------------------------------

def _vectors_with_one_norm(n, k):
    """Every a in Z^n with sum |a_j| = k, enumerated exactly."""
    if n == 0:
        if k == 0:
            yield ()
        return
    for first in range(-k, k + 1):
        for rest in _vectors_with_one_norm(n - 1, k - abs(first)):
            yield (first,) + rest


def one_norm(a):
    """||mu||_1 = sum |a_j|, Lauret v5 eq (2.6).  Named so the census mutation
    can replace it without touching anything else."""
    return sum(abs(aj) for aj in a)


def zero_count(a):
    """Z(mu) = #{ j : a_j = 0 }, Lauret v5 eq (2.6).  Named for the same reason:
    Z(mu) feeds the Lemma 3.2 kernel through ell, so a miscount here is a
    lattice-census defect and must redden a different check than a wrong
    generator or a wrong branch index."""
    return sum(1 for aj in a if aj == 0)


def lattice_counts(q, s, kmax, blocks=None):
    """N_L(k, ell) for 0 <= k <= kmax.

    Lauret v5 eq (2.6) and (2.7):
        ||mu||_1 = sum |a_j|,   Z(mu) = #{ j : a_j = 0 }
        N_L(k, ell) = #{ mu in L : ||mu||_1 = k, Z(mu) = ell }
    """
    if blocks is None:
        blocks = generator_block_parameters(s)
    n = len(s)
    counts = {}
    for k in range(kmax + 1):
        for a in _vectors_with_one_norm(n, k):
            if not in_congruence_lattice(a, q, blocks):
                continue
            counts[(one_norm(a), zero_count(a))] = (
                counts.get((one_norm(a), zero_count(a)), 0) + 1)
    return counts


# --- layer 3: Lemma 3.2, the weight multiplicity -----------------------------

def _lemma_3_2_kernel(r, p, ell, n, binom=paper_binom):
    """The inner multi-sum of Lemma 3.2, as a function of r, p, Z(mu)=ell, n.

    Lemma 3.2 depends on mu only through ||mu||_1 and Z(mu), which is why the
    same kernel serves eq (3.1) with ell in place of Z(mu).  Transcribed term
    for term; the j/t ordering difference between the two printings is addition
    and does not change the value.
    """
    total = 0
    for j in range(1, p + 1):
        sign = -1 if (j - 1) % 2 else 1
        for t in range(0, (p - j) // 2 + 1):
            c1 = binom(n - p + j + 2 * t, t)
            if c1 == 0:
                continue
            for beta in range(0, p - j - 2 * t + 1):
                c2 = 2 ** (p - j - 2 * t - beta)
                c3 = binom(n - ell, beta)
                c4 = binom(ell, p - j - 2 * t - beta)
                if c3 == 0 or c4 == 0:
                    continue
                for alpha in range(0, beta + 1):
                    c5 = binom(beta, alpha)
                    if c5 == 0:
                        continue
                    for i in range(0, j):
                        c6 = binom(r - i - p + j + alpha + t + n - 2, n - 2)
                        total += sign * c1 * c2 * c3 * c4 * c5 * c6
    return total


def weight_multiplicity(k, p, one_norm, zeros, n, binom=paper_binom):
    """m_{pi_{k,p}}(mu), Lemma 3.2.

    r(mu) = (k + p - ||mu||_1)/2 must be a non-negative integer, else 0.
    """
    num = k + p - one_norm
    if num < 0 or num % 2 != 0:
        return 0
    return _lemma_3_2_kernel(num // 2, p, zeros, n, binom=binom)


# --- layer 4: M_Gamma(k, p), eq (3.1) ----------------------------------------

def M_Gamma(k, p, q, s, counts=None, binom=paper_binom, blocks=None):
    """M_Gamma(k, p), Lauret v5 eq (3.1).

    Equals dim V^Gamma_{pi_{k-1, p}}; that identity is the primary evidence for
    correction B and is recorded in THEOREM_TRANSCRIPTION.md section 5.
    """
    n = len(s)
    top = k - 1 + p
    if counts is None:
        counts = lattice_counts(q, s, max(top, 0), blocks=blocks)
    total = 0
    for ell in range(0, n + 1):
        for r in range(0, top // 2 + 1):
            cnt = counts.get((top - 2 * r, ell), 0)
            if cnt == 0:
                continue
            total += cnt * _lemma_3_2_kernel(r, p, ell, n, binom=binom)
    return total


# --- layer 5: Theorem 3.3 branch mapping (CORRECTION B lives here) -----------

def _branch_indices(p, mapping):
    """Which M_Gamma index each spectral branch reads.

    corrected   mult(lambda_{k,p-1}) = M_Gamma(k, p)     mult(lambda_{k,p}) = M_Gamma(k, p+1)
    as_printed  mult(lambda_{k,p-1}) = M_Gamma(k, p-1)   mult(lambda_{k,p}) = M_Gamma(k, p)

    `as_printed` exists ONLY as mutation 2.  It is the mistake an implementer
    working straight from the paper would make, and at Gamma = 1 it returns
    M_Gamma(k, 0) = 0 on the lower branch, because eq (3.1)'s j-sum is empty at
    p = 0.  It must fail the Gamma = 1 gate.
    """
    if mapping == "corrected":
        return p, p + 1
    if mapping == "as_printed":
        return p - 1, p
    raise ValueError(f"unknown branch mapping {mapping!r}")


def eigenvalue(k, p, n):
    """lambda_{k,p}, Lauret v5 eq (2.2)."""
    if p == -1:
        return 0
    return (k + p) * (k + 2 * n - 2 - p)


# --- layer 6: the reported multiplicities ------------------------------------

def p_form_spectrum(p, q, s, kmax, mapping="corrected", binom=paper_binom,
                    blocks=None):
    """The p-spectrum of Gamma\\S^(2n-1) for k = 1..kmax.

    Returns a list of records, one per k, each carrying BOTH branches with the
    eigenvalue, the M_Gamma index actually read, and the multiplicity.  The
    index is reported rather than hidden so that a reader can see which branch
    mapping produced the number without rerunning anything.
    """
    n = len(s)
    if not 0 <= p <= n - 1:
        raise ValueError(f"Theorem 3.3 requires 0 <= p <= n-1; got p={p}, n={n}")
    if blocks is None:
        blocks = generator_block_parameters(s)
    lo_idx, hi_idx = _branch_indices(p, mapping)
    top = kmax - 1 + max(lo_idx, hi_idx)
    counts = lattice_counts(q, s, max(top, 0), blocks=blocks)

    out = []
    for k in range(1, kmax + 1):
        lo = M_Gamma(k, lo_idx, q, s, counts=counts, binom=binom, blocks=blocks)
        hi = M_Gamma(k, hi_idx, q, s, counts=counts, binom=binom, blocks=blocks)
        out.append({
            "k": k,
            "lower_branch": {
                "eigenvalue": eigenvalue(k, p - 1, n),
                "M_Gamma_index": lo_idx,
                "multiplicity": lo,
            },
            "upper_branch": {
                "eigenvalue": eigenvalue(k, p, n),
                "M_Gamma_index": hi_idx,
                "multiplicity": hi,
            },
        })
    return out
