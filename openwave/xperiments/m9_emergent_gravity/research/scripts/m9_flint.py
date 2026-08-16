"""ARB/Flint backend for *future* Peschel. Do not replay Papers 70--73.

Priority 1: same math as mpmath at dps=80.
Priority 2: compiled ARB eigensolve when a new script needs it.
"""

from __future__ import annotations

try:
    import flint

    HAVE_FLINT = True
except ImportError:
    flint = None
    HAVE_FLINT = False

DPS = 80


def set_dps(dps: int = DPS) -> None:
    if not HAVE_FLINT:
        raise RuntimeError("python-flint is not installed")
    flint.ctx.dps = dps


def peschel_s(block) -> float:
    """von Neumann entropy of a real symmetric correlator block.

    ``block`` is a square sequence of sequences (mpmath matrix,
    list of lists, or arb_mat). Returns a Python float of the
    mid-point; the work is ARB at ``flint.ctx.dps``.
    """
    if not HAVE_FLINT:
        raise RuntimeError("python-flint is not installed")
    set_dps()
    n = block.rows if hasattr(block, "rows") else len(block)
    mat = flint.arb_mat(n, n)
    for i in range(n):
        for j in range(n):
            mat[i, j] = flint.arb(str(block[i, j] if hasattr(block, "rows") else block[i][j]))
    ev = mat.eig()
    # eig may return acb; take real mid
    acc = flint.arb(0)
    clip = flint.arb("1e-30")
    one = flint.arb(1)
    for k in range(n):
        w = ev[k]
        if hasattr(w, "real"):
            w = w.real
        if w < clip:
            w = clip
        if w > one - clip:
            w = one - clip
        acc -= w * w.log() + (one - w) * (one - w).log()
    return float(acc)
