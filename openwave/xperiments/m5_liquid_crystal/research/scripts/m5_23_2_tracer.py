"""M5.23.2 arm (1) — the general disclination-line TRACER (research instrument).

Finds defect lines in a live M field instead of assuming them from the seed
geometry: rods, ring cords, and split-vortex loops become MEASURED objects
(count, closure, extent), the closure instrument the M5.21.10 census could
not provide (its fixed-radii biaxial census was cut-sensitive at the grade
that mattered: 1-cell filament doublets reconnecting past the r-cut).

DETECTION (the Stage D criterion, measured at the M5.23 close): the
disclination core is an exact UNIAXIAL ESCAPE of the spatial spectrum —
with eigenvalues sorted ascending (l1 <= l2 <= l3) the minor split
s = l2 - l1 is 0.000 on the rod axis vs ~ delta (0.265 measured) in the
biaxial bulk. The tracer thresholds s below a SELF-CALIBRATED bulk scale
(median over the domain) — no hardcoded arena radii (the M5.21.10 audit
requirement (c)).

SCOPE NOTE: the criterion is meaningful on BIAXIAL-family states (far field
D = diag(1, delta, 0), split = delta). A uniaxial far field (1, d, d) has
split 0 EVERYWHERE and would light up whole-domain: trace() warns when the
bulk scale collapses.

ASSEMBLY (audit requirement (a)): connected components of the detection
mask under 26-connectivity — lines by connectivity, NOT blob components
inside a radius cut.

CLOSURE (audit requirement (b)): explicit, per component:
  - "boundary"     — touches the domain edge (run-out; in-box half-lines)
  - "closed-loop"  — Euler characteristic chi = 0 of the voxel cubical
                     complex (a solid torus has chi 0; first Betti b1 = 1)
  - "open"         — chi = 1 (contractible: a rod / filament segment)
  - "chi=<n>"      — anything else (multi-handle / cavity; inspect by hand)

Consumers: the m5_23_2 selftest (T gates), the M5.21.10 C-pair two-loop
read (m5_23_2_a_cpair_read.py), the M5.8.8 rod-localization check (future).
Record: research/tasks/m5_23_2_task_details.md.
"""

import numpy as np
from scipy import ndimage

STRUCT26 = np.ones((3, 3, 3), dtype=bool)


def split_map(M):
    """The uniaxial-escape order parameter s = l2 - l1 (ascending eigvalsh).

    M: (..., 3, 3) spatial blocks, or (..., 4, 4) (spatial block sliced).
    Returns s with the leading grid shape.
    """
    if M.shape[-1] == 4:
        M = M[..., 1:, 1:]
    ev = np.linalg.eigvalsh(M)  # ascending
    return ev[..., 1] - ev[..., 0]


def euler_characteristic(mask):
    """chi = V - E + F - C of the cubical complex of the occupied voxels.

    Each voxel is a closed unit cube; shared vertices/edges/faces counted
    once (union complex). Ball -> 1, solid torus -> 0.
    """
    occ = np.ascontiguousarray(mask, dtype=bool)
    p = np.pad(occ, 1)
    C = int(occ.sum())
    # faces: exists if either adjacent cube (along its normal) is occupied
    F = 0
    for ax in range(3):
        a = np.swapaxes(p, 0, ax)
        F += int((a[:-1] | a[1:]).sum())
    # edges: exists if any of the 4 cubes around it is occupied
    E = 0
    for ax in range(3):  # edge direction ax; neighbors span the other two axes
        a = np.moveaxis(p, ax, 2)  # edge axis last
        e = a[:-1, :-1] | a[:-1, 1:] | a[1:, :-1] | a[1:, 1:]
        E += int(e.sum())
    # vertices: exists if any of the 8 cubes around it is occupied
    v = (
        p[:-1, :-1, :-1]
        | p[:-1, :-1, 1:]
        | p[:-1, 1:, :-1]
        | p[:-1, 1:, 1:]
        | p[1:, :-1, :-1]
        | p[1:, :-1, 1:]
        | p[1:, 1:, :-1]
        | p[1:, 1:, 1:]
    )
    V = int(v.sum())
    return V - E + F - C


def trace(M, h=1.0, thr=None, thr_frac=0.35, min_size=1):
    """Trace defect lines in a spatial M field.

    Args:
        M: (n, n, n, 3, 3) or (n, n, n, 4, 4) field array.
        h: grid spacing (research units) for the length estimates.
        thr: absolute split threshold; None -> thr_frac * bulk (median).
        thr_frac: threshold as a fraction of the self-calibrated bulk scale.
        min_size: drop components smaller than this many voxels.

    Returns dict:
        bulk_split, thr, warn, split (the map), mask,
        split_min_on_lines, lines = [{id, n_vox, verdict, chi,
        touches_boundary, centroid, bbox, extent_h, length_est_h}]
        sorted by n_vox descending.
    """
    s = split_map(np.asarray(M, dtype=np.float64))
    bulk = float(np.median(s))
    warn = []
    if bulk < 1e-3:
        warn.append(
            f"bulk split scale ~ 0 ({bulk:.2e}): the far field is not "
            "biaxial-family; the uniaxial-escape criterion is not meaningful here"
        )
    t = float(thr) if thr is not None else thr_frac * bulk
    mask = s < t
    labels, n_comp = ndimage.label(mask, structure=STRUCT26)
    lines = []
    n = s.shape[0]
    smin_global = float("inf")
    for lid in range(1, n_comp + 1):
        comp = labels == lid
        n_vox = int(comp.sum())
        if n_vox < min_size:
            continue
        idx = np.argwhere(comp)
        touches = bool(
            (idx.min() == 0)
            or (idx[:, 0].max() == n - 1)
            or (idx[:, 1].max() == s.shape[1] - 1)
            or (idx[:, 2].max() == s.shape[2] - 1)
        )
        chi = euler_characteristic(comp)
        if touches:
            verdict = "boundary"
        elif chi == 0:
            verdict = "closed-loop"
        elif chi == 1:
            verdict = "open"
        else:
            verdict = f"chi={chi}"
        lo, hi = idx.min(axis=0), idx.max(axis=0)
        extent = float(np.linalg.norm((hi - lo).astype(float)) * h)
        smin = float(s[comp].min())
        smin_global = min(smin_global, smin)
        lines.append(
            {
                "id": lid,
                "n_vox": n_vox,
                "verdict": verdict,
                "chi": int(chi),
                "touches_boundary": touches,
                "centroid": [float(c * h) for c in idx.mean(axis=0)],
                "bbox": [lo.tolist(), hi.tolist()],
                "extent_h": extent,
                "length_est_h": float(n_vox * h),  # thin-line upper estimate
                "split_min": smin,
            }
        )
    lines.sort(key=lambda d: -d["n_vox"])
    return {
        "bulk_split": bulk,
        "thr": t,
        "warn": warn,
        "split": s,
        "mask": mask,
        "split_min_on_lines": smin_global if lines else float("nan"),
        "lines": lines,
    }


if __name__ == "__main__":
    import argparse
    import os

    ap = argparse.ArgumentParser(description="Trace defect lines in a saved M npz")
    ap.add_argument("npz")
    ap.add_argument("--key", default="M")
    ap.add_argument("--h", type=float, default=None, help="grid unit (default: npz 'h' or 1)")
    ap.add_argument("--thr-frac", type=float, default=0.35)
    args = ap.parse_args()
    z = np.load(args.npz)
    h_arg = args.h if args.h is not None else float(z["h"]) if "h" in z.files else 1.0
    r = trace(z[args.key], h=h_arg, thr_frac=args.thr_frac)
    print(f"{os.path.basename(args.npz)}[{args.key}]: bulk {r['bulk_split']:.4f}, thr {r['thr']:.4f}")
    for w in r["warn"]:
        print(f"  WARNING: {w}")
    for ln in r["lines"]:
        print(
            f"  line {ln['id']}: {ln['n_vox']} vox, {ln['verdict']} (chi {ln['chi']}), "
            f"extent {ln['extent_h']:.1f} h-units, centroid {np.round(ln['centroid'], 1)}"
        )
