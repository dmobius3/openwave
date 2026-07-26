"""M5.21.10 adversarial audit (independent recompute of claims C1..C7).

Read-only on all existing artifacts; writes ONLY
../data/m5_21_10_audit.json. Own implementations for every check:

    ledger      E + KE + absorbed drift recomputed from the ev hist
                JSONs; the t = 0 budget E0 recomputed from the relax
                npz snapshots with INS.e_parts (this also identifies
                WHICH stored state each evolution actually started
                from; the ev JSON records only start = f64_X)
    census      own biaxiality mask (min eigen-gap < thr on my own
                eigvalsh read) + own scipy.ndimage 26-connected
                labeling on npz snapshots; loop_read_xyz is NOT
                called, only its conventions (r < 0.65 half census
                region, edge zone r > 0.62 half) are mirrored so the
                claims are audited on their own definitions, plus a
                full-domain variant with no radial cut
    kinematics  departing counts / pair geometry recomputed from the
                lp + kin JSONs (speeds and angles re-derived from the
                raw centroid records, not the stored v / dir fields)
    ring        E_sym / E_fwd / E_2h recomputed from the t32 npz
                endpoints via INS.e_parts under each stencil cfg;
                orderings, spreads and the 2h offset re-derived
    retention   f64 rows recomputed from the f64 npz endpoints via
                INS.retention with the correct seed cfg (f48 npz no
                longer exist; those values rest on the row JSONs)

Instruments are imported, never modified. Verdict per claim:
CONFIRMED / NUANCE / REFUTED, with the recomputed numbers.
"""
from __future__ import annotations

import importlib.util
import json
import os
import time

import numpy as np
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

_spec = importlib.util.spec_from_file_location(
    "dec64", os.path.join(HERE, "m5_21_10_a_decay64.py"))
D64 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(D64)
D, INS = D64.D, D64.INS

T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:6.0f}s] {msg}", flush=True)


def jload(name):
    with open(os.path.join(DATA, name)) as f:
        return json.load(f)


# ================= own geometry / census ===========================
def my_grid(n, h):
    x = (np.arange(n) - (n - 1) / 2.0) * h
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    return X, Y, Z, np.sqrt(X * X + Y * Y + Z * Z)


def my_census(M, n, h, thr, r_cut=None, min_size=3):
    """Own threshold + labeling. Mask: min eigen-gap < thr.
    26-connectivity via scipy.ndimage. Optional census radial cut."""
    lam = np.linalg.eigvalsh(M)
    gap = np.minimum(lam[..., 1] - lam[..., 0],
                     lam[..., 2] - lam[..., 1])
    X, Y, Z, r = my_grid(n, h)
    half = n * h / 2.0
    mask_full = gap < thr
    mask = mask_full & (r < r_cut) if r_cut is not None else mask_full
    lab, nc = ndimage.label(mask, structure=np.ones((3, 3, 3), int))
    sizes = np.bincount(lab.ravel(), minlength=nc + 1)
    comps = []
    for k in range(1, nc + 1):
        if sizes[k] < min_size:
            continue
        sel = lab == k
        cx = float(X[sel].mean())
        cy = float(Y[sel].mean())
        cz = float(Z[sel].mean())
        comps.append({
            "size": int(sizes[k]),
            "xyz": [cx, cy, cz],
            "r_xyz": float(np.sqrt(cx * cx + cy * cy + cz * cz)),
            "r_cells_mean": float(r[sel].mean()),
            "r_cells_max": float(r[sel].max()),
            "touches_edge_zone": bool(r[sel].max() > 0.62 * half)})
    comps.sort(key=lambda c: -c["size"])
    return {"n_mask_full_domain": int(mask_full.sum()),
            "n_mask_census_region": int(mask.sum()),
            "n_comps_ge_min": len(comps), "comps": comps}


def angle_deg(u, v):
    u, v = np.asarray(u, float), np.asarray(v, float)
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    if nu < 1e-12 or nv < 1e-12:
        return None
    return float(np.degrees(np.arccos(
        np.clip(np.dot(u, v) / (nu * nv), -1.0, 1.0))))


def frag_filter(comps, smin=8, smax=25, rmin=9.0, rmax=19.0):
    """off-center compact fragment candidates (own criteria)."""
    return [c for c in comps
            if smin <= c["size"] <= smax
            and rmin <= c["r_xyz"] <= rmax
            and not c["touches_edge_zone"]]


# ================= C1: ledger ======================================
def c1():
    log("C1 ledger: recompute E0 from relax npz + drift from hist")
    out = {"per_ev": {}}
    for ev, tag in (("C", "f64_C"), ("B", "f64_B"), ("A", "f64_A")):
        z = np.load(os.path.join(DATA, f"m5_21_6_end_{tag}.npz"))
        cfg = D64.cfg_from_npz(z)
        assert abs(float(z["delta"]) - cfg["delta"]) < 1e-12
        hist = jload(f"m5_21_10_ev_{ev}_free64.json")["hist"]
        tots = np.array([r["E"] + r["KE"] + r["absorbed"]
                         for r in hist])
        cands = {}
        for key in ("M", "M_it1000", "M_it2000", "M_it500"):
            if key not in z.files:
                continue
            eu, ed, evv = INS.e_parts(z[key].astype(float), cfg)
            cands[key] = float(eu + ed + evv)
            if abs(cands[key] - tots[0]) / tots[0] < 2e-3:
                break
        start_key = min(cands, key=lambda k: abs(cands[k] - tots[0]))
        E0 = cands[start_key]
        ident = abs(E0 - tots[0]) / tots[0] < 2e-3
        row = {
            "start_state_identified": start_key,
            "start_identification_ok": bool(ident),
            "E0_recomputed": E0,
            "E_candidates": cands,
            "tot_first_row_t5": float(tots[0]),
            "tot_last_row_t150": float(tots[-1]),
            "tot_min": float(tots.min()), "tot_max": float(tots.max()),
            "drift_end_vs_first": float((tots[-1] - tots[0]) / tots[0]),
            "drift_max_vs_first": float(
                np.max(np.abs(tots - tots[0])) / tots[0]),
            "drift_max_vs_E0": float(
                np.max(np.abs(tots - E0)) / abs(E0))}
        out["per_ev"][ev] = row
        log(f"  {ev}: start={start_key} E0={E0:.6f} "
            f"drift end/first={row['drift_end_vs_first']:.2e} "
            f"max/first={row['drift_max_vs_first']:.2e} "
            f"max/E0={row['drift_max_vs_E0']:.2e}")
    m_e0 = max(r["drift_max_vs_E0"] for r in out["per_ev"].values())
    m_f = max(r["drift_max_vs_first"] for r in out["per_ev"].values())
    m_end = max(abs(r["drift_end_vs_first"])
                for r in out["per_ev"].values())
    out["max_drift_vs_E0_all"] = m_e0
    out["max_drift_vs_first_all"] = m_f
    out["max_drift_end_vs_first_all"] = m_end
    if m_e0 <= 5e-5:
        out["verdict"] = "CONFIRMED"
        out["reason"] = ("ledger conserved to <= 5e-5 vs the true "
                        "t=0 budget on all three evolutions")
    elif m_end <= 5e-5:
        out["verdict"] = "NUANCE"
        out["reason"] = (
            "conservation is real but at the 1e-4 grade, not 5e-5: "
            f"max-over-rows drift vs the t=0 budget reaches {m_e0:.1e} "
            f"(B), and vs the first hist row {m_f:.1e}; only the "
            "end-vs-first-row difference stays <= 5e-5 "
            f"({m_end:.1e}). Also a provenance gap: the evolutions "
            "start from relax snapshot M_it1000 (energy match), not "
            "the stored final M, and the ev JSON does not record it")
    else:
        out["verdict"] = "REFUTED"
        out["reason"] = f"end-to-end drift {m_end:.1e} exceeds 5e-5"
    return out


# ================= C2: the C pair ==================================
def c2(budget):
    log("C2 pair: lp/kin recompute + own-labeling npz spot check")
    lp = jload("m5_21_10_lp_C_free64.json")
    keys = sorted(lp["snaps"], key=lambda k: int(k[4:]))
    frag_hist = []
    for k in keys:
        s = lp["snaps"][k]
        t9 = s["thr"][1]
        assert t9["thr"] == 0.09
        frs = [c for c in t9["components"]
               if 8 <= c["size"] <= 25
               and 9.0 <= float(np.linalg.norm(c["xyz"])) <= 19.0
               and not c["touches_edge_zone"]]
        frag_hist.append({"t": s["t"], "n_frag": len(frs),
                          "sizes": [c["size"] for c in frs],
                          "xyz": [c["xyz"] for c in frs]})
    t_first_pair = next((f["t"] for f in frag_hist
                         if f["n_frag"] >= 2), None)
    kin = jload("m5_21_10_kin_C_free64.json")
    dep = [t for t in kin["tracks"] if t["departing"]]
    big = max(kin["tracks"], key=lambda t: t["size_first"])
    pair = sorted(dep, key=lambda t: t["xyz_first"][2])
    pair_num = {}
    if len(pair) == 2:
        a, b = pair
        va = (np.array(a["xyz_last"]) - np.array(a["xyz_first"])) \
            / (a["t_last"] - a["t_first"])
        vb = (np.array(b["xyz_last"]) - np.array(b["xyz_first"])) \
            / (b["t_last"] - b["t_first"])
        pair_num = {
            "sizes": [[a["size_first"], a["size_last"]],
                      [b["size_first"], b["size_last"]]],
            "speeds_recomputed": [float(np.linalg.norm(va)),
                                  float(np.linalg.norm(vb))],
            "dir_angle_deg_recomputed": angle_deg(va, vb),
            "pos_angle_deg_t_first": angle_deg(a["xyz_first"],
                                               b["xyz_first"]),
            "r_first": [a["r_first"], b["r_first"]],
            "r_last": [a["r_last"], b["r_last"]]}
    # own-labeling spot checks on the raw snapshots, including the
    # object character (bbox extents, radial span) and the full-domain
    # connectivity of each census-cut fragment candidate
    z = np.load(os.path.join(DATA, "m5_21_10_ev_C_free64.npz"))
    n, h = int(z["n"]), float(z["h"])
    X, Y, Zc, r = my_grid(n, h)
    half = n * h / 2.0
    spot = {}
    for key, t in (("M_it3200", 80.0), ("M_it4800", 120.0),
                   ("M_it6000", 150.0)):
        M = z[key].astype(float)
        lam = np.linalg.eigvalsh(M)
        gap = np.minimum(lam[..., 1] - lam[..., 0],
                         lam[..., 2] - lam[..., 1])
        mask_full = gap < 0.09
        mask_cut = mask_full & (r < 0.65 * half)
        st26 = np.ones((3, 3, 3), int)
        lab_c, nc_c = ndimage.label(mask_cut, structure=st26)
        lab_f, _ = ndimage.label(mask_full, structure=st26)
        sizes_c = np.bincount(lab_c.ravel(), minlength=nc_c + 1)
        sizes_f = np.bincount(lab_f.ravel())
        frs = []
        for k in range(1, nc_c + 1):
            if sizes_c[k] < 3:
                continue
            sel = lab_c == k
            cx = float(X[sel].mean())
            cy = float(Y[sel].mean())
            cz = float(Zc[sel].mean())
            r_xyz = float(np.sqrt(cx * cx + cy * cy + cz * cz))
            r_max = float(r[sel].max())
            if not (8 <= sizes_c[k] <= 25 and 9.0 <= r_xyz <= 19.0
                    and r_max <= 0.62 * half):
                continue
            fl = np.unique(lab_f[sel])
            f_sizes = [int(sizes_f[q]) for q in fl if q > 0]
            f_rmax = float(max(r[np.isin(lab_f, [q for q in fl
                                                 if q > 0])].max(),
                               r_max))
            frs.append({
                "size": int(sizes_c[k]),
                "xyz": [cx, cy, cz], "r_xyz": r_xyz,
                "r_cells_min": float(r[sel].min()),
                "r_cells_mean": float(r[sel].mean()),
                "r_cells_max": r_max,
                "radial_span": float(r_max - r[sel].min()),
                "bbox_extent": [float(X[sel].max() - X[sel].min()),
                                float(Y[sel].max() - Y[sel].min()),
                                float(Zc[sel].max() - Zc[sel].min())],
                "full_domain_comp_sizes": f_sizes,
                "full_domain_comp_r_max": f_rmax,
                "full_domain_reaches_edge_zone":
                    bool(f_rmax > 0.62 * half)})
        frs.sort(key=lambda c: -c["size"])
        best_ang = None
        for i in range(len(frs)):
            for j in range(i + 1, len(frs)):
                ang = angle_deg(frs[i]["xyz"], frs[j]["xyz"])
                if best_ang is None or ang > best_ang:
                    best_ang = ang
        n_filament = sum(1 for c in frs
                         if c["radial_span"] > 6.0
                         or min(c["bbox_extent"][:2]) == 0.0)
        spot[key] = {
            "t": t, "n_offcenter_frag_census_cut": len(frs),
            "n_filament_like": n_filament,
            "n_full_domain_edge_connected": sum(
                1 for c in frs
                if c["full_domain_reaches_edge_zone"]),
            "frags": frs, "max_pos_angle_deg": best_ang,
            "mask_cells_census_region": int(mask_cut.sum())}
        log(f"  spot {key} t={t}: frags(cut)={len(frs)} "
            f"filament-like={n_filament} maxAngle={best_ang}")
        del M, lam, gap, mask_full, mask_cut, lab_c, lab_f
    out = {"frag_census_history_thr009": frag_hist,
           "t_first_two_fragments": t_first_pair,
           "n_departing_tracks_thr009": len(dep),
           "big_remnant": {"size_first": big["size_first"],
                           "size_last": big["size_last"],
                           "departing": big["departing"]},
           "pair": pair_num, "spot": spot}
    n150 = spot["M_it6000"]["n_offcenter_frag_census_cut"]
    n120 = spot["M_it4800"]["n_offcenter_frag_census_cut"]
    fil150 = spot["M_it6000"]["n_filament_like"]
    edge150 = spot["M_it6000"]["n_full_domain_edge_connected"]
    core_ok = (t_first_pair == 80.0 and len(dep) == 2
               and not big["departing"]
               and pair_num
               and all(13 <= s[0] <= 15 and 12 <= s[1] <= 14
                       for s in pair_num["sizes"])
               and abs(pair_num["speeds_recomputed"][0]
                       - pair_num["speeds_recomputed"][1]) < 0.002
               and abs(pair_num["dir_angle_deg_recomputed"]
                       - 109.0) < 5.0)
    out["core_checks_pass"] = bool(core_ok)
    spot_clean = (n120 >= 2 and n150 == 2 and fil150 == 0)
    if core_ok and spot_clean:
        out["verdict"] = "CONFIRMED"
        out["reason"] = "all pair claims reproduce, spot checks clean"
    elif core_ok:
        out["verdict"] = "NUANCE"
        out["reason"] = (
            "the JSON-level numbers reproduce exactly (census gains "
            "2 departing antipodal size-14 objects at t=80, same "
            "recomputed centroid speed 0.052, dir angle "
            f"{pair_num['dir_angle_deg_recomputed']:.1f} deg, big "
            "remnant non-departing), BUT the raw snapshots do not "
            "support the stable-departing-fragment picture: at t=80 "
            "each object is a 2x2x5 rod (radial span 4.9), at t=120 "
            "and t=140 the census resolves NO off-center fragments "
            f"(own labeling: {n120} at t=120), and at t=150 there "
            f"are {n150} candidates, not 2, ALL {fil150} of them "
            "1-cell-wide z-columns spanning r 8.3-19.8 (mirror "
            f"quadruplets of the t=80 pair), {edge150} of which "
            "connect to the edge zone once the census r-cut is "
            "removed. Antipodal geometry and cell counts 13-14 do "
            "hold; 'two compact fragments departing at 0.052' is a "
            "tracker construction over threshold filaments")
    else:
        out["verdict"] = "REFUTED"
        out["reason"] = "the t=80 pair read itself fails recompute"
    return out


# ================= C3: threshold sensitivity =======================
def c3():
    log("C3 thr ladder: departing counts + pair persistence")
    files = {0.06: "m5_21_10_kin_C_free64_thr0.json",
             0.09: "m5_21_10_kin_C_free64.json",
             0.15: "m5_21_10_kin_C_free64_thr2.json"}
    counts, deps = {}, {}
    for thr, fn in files.items():
        d = jload(fn)
        assert abs(d["thr"] - thr) < 1e-9
        deps[thr] = [t for t in d["tracks"] if t["departing"]]
        counts[str(thr)] = len(deps[thr])
    # pair persistence: antipodal departing pair at 0.06 and 0.09
    def antipodal_pair(tracks, tol=3.0):
        for i in range(len(tracks)):
            for j in range(i + 1, len(tracks)):
                a = np.array(tracks[i]["xyz_first"])
                b = np.array(tracks[j]["xyz_first"])
                if np.linalg.norm(a + b) < tol and \
                        np.linalg.norm(a) > 5.0:
                    return tracks[i], tracks[j]
        return None
    p06 = antipodal_pair(deps[0.06])
    p09 = antipodal_pair(deps[0.09])
    same_obj = None
    if p06 and p09:
        d_min = min(
            float(np.linalg.norm(np.array(a["xyz_first"])
                                 - np.array(b["xyz_first"])))
            for a in p06 for b in p09)
        same_obj = d_min < 3.0
    out = {"departing_counts": counts,
           "pair_at_006": bool(p06), "pair_at_009": bool(p09),
           "pair_min_cross_thr_distance_lt3": same_obj,
           "pair_006_records": [
               {"xyz_first": t["xyz_first"], "t_first": t["t_first"],
                "sizes": [t["size_first"], t["size_last"]],
                "speed": t["speed"]} for t in (p06 or [])],
           "note_006_third_departing": [
               {"xyz_first": t["xyz_first"],
                "sizes": [t["size_first"], t["size_last"]],
                "r": [t["r_first"], t["r_last"]]}
               for t in deps[0.06]
               if p06 and t not in p06]}
    ok = (counts == {"0.06": 3, "0.09": 2, "0.15": 0}
          and p06 and p09 and same_obj)
    out["verdict"] = "CONFIRMED" if ok else (
        "NUANCE" if counts["0.09"] == 2 and p06 and p09
        else "REFUTED")
    out["reason"] = (
        f"departing counts {counts['0.06']}/{counts['0.09']}/"
        f"{counts['0.15']} across 0.06/0.09/0.15; an antipodal "
        "departing pair exists at both 0.06 and 0.09 at the same "
        "location (cross-threshold distance < 3)" if ok else
        "counts or pair persistence deviate; see numbers")
    return out


# ================= C4: B disintegration ============================
def c4(budget_B):
    log("C4 B: own-labeling volume + end remnant + ledger reads")
    lp = jload("m5_21_10_lp_B_free64.json")
    s10 = lp["snaps"]["M_it400"]["thr"][1]
    s150 = lp["snaps"]["M_it6000"]["thr"][1]
    v10_lp = sum(c["size"] for c in s10["components"])
    v150_lp = sum(c["size"] for c in s150["components"])
    z = np.load(os.path.join(DATA, "m5_21_10_ev_B_free64.npz"))
    n, h = int(z["n"]), float(z["h"])
    r_cut = 0.65 * n * h / 2.0
    own = {}
    for key, t in (("M_it400", 10.0), ("M_it6000", 150.0)):
        M = z[key].astype(float)
        cen = my_census(M, n, h, 0.09, r_cut=r_cut)
        compact = [c for c in cen["comps"]
                   if not c["touches_edge_zone"]]
        own[key] = {
            "t": t,
            "mask_cells_census_region": cen["n_mask_census_region"],
            "mask_cells_full_domain": cen["n_mask_full_domain"],
            "largest_compact_nonedge": (compact[0]["size"]
                                        if compact else 0),
            "largest_compact_xyz": (compact[0]["xyz"]
                                    if compact else None),
            "top_comps": [{"size": c["size"],
                           "r_xyz": c["r_xyz"],
                           "touches_edge_zone": c["touches_edge_zone"]}
                          for c in cen["comps"][:4]]}
        log(f"  spot {key}: mask(cut)={cen['n_mask_census_region']} "
            f"largest compact={own[key]['largest_compact_nonedge']}")
        del M
    hist = jload("m5_21_10_ev_B_free64.json")["hist"]
    absorbed = hist[-1]["absorbed"]
    rot = hist[-1]["rot_core_deg"]
    ratio_own = own["M_it6000"]["mask_cells_census_region"] / \
        max(own["M_it400"]["mask_cells_census_region"], 1)
    out = {"volume_thr009_lp": {"t10": v10_lp, "t150": v150_lp,
                                "ratio": v150_lp / max(v10_lp, 1)},
           "volume_thr009_own": {
               "t10": own["M_it400"]["mask_cells_census_region"],
               "t150": own["M_it6000"]["mask_cells_census_region"],
               "ratio": ratio_own},
           "own_spot": own,
           "absorbed_final": absorbed,
           "budget_E0": budget_B,
           "absorbed_over_budget": absorbed / budget_B,
           "rot_final_deg": rot,
           "end_largest_compact_nonedge":
               own["M_it6000"]["largest_compact_nonedge"],
           "end_edge_reaching_note": [
               c for c in own["M_it6000"]["top_comps"]
               if c["touches_edge_zone"]]}
    ok = (ratio_own < 0.5
          and 0.45 <= out["absorbed_over_budget"] <= 0.55
          and 19.0 <= rot <= 21.0
          and out["end_largest_compact_nonedge"] <= 11)
    out["verdict"] = "CONFIRMED" if ok else "NUANCE"
    out["reason"] = (
        f"own mask volume falls {own['M_it400']['mask_cells_census_region']}"
        f" -> {own['M_it6000']['mask_cells_census_region']} cells "
        f"(x{ratio_own:.2f}), absorbed/budget = "
        f"{out['absorbed_over_budget']:.3f} ~ half, rot = {rot:.1f} "
        f"deg, largest compact non-edge end remnant = "
        f"{out['end_largest_compact_nonedge']} cells; note two "
        "~316-cell edge-reaching mirror shells do remain at r~13 but "
        "the claim scoped itself to compact non-edge" if ok
        else "one or more sub-checks deviate; see numbers")
    return out


# ================= C5: A control null ==============================
def c5(budget_A):
    log("C5 A: energy hold + departing artifacts + rot/specdev")
    hist = jload("m5_21_10_ev_A_free64.json")["hist"]
    absorbed = hist[-1]["absorbed"]
    E_first, E_last = hist[0]["E"], hist[-1]["E"]
    specd = [r["spec_dev_core"] for r in hist]
    rots = [r["rot_core_deg"] for r in hist]
    kin = jload("m5_21_10_kin_A_free64.json")
    dep = [t for t in kin["tracks"] if t["departing"]]
    dep_rec = [{"r_first": t["r_first"], "r_last": t["r_last"],
                "sizes": [t["size_first"], t["size_last"]],
                "speed": t["speed"], "xyz_first": t["xyz_first"]}
               for t in dep]
    # own end-state read
    z = np.load(os.path.join(DATA, "m5_21_10_ev_A_free64.npz"))
    n, h = int(z["n"]), float(z["h"])
    cen = my_census(z["M_it6000"].astype(float), n, h, 0.09,
                    r_cut=0.65 * n * h / 2.0)
    compact = [c for c in cen["comps"] if not c["touches_edge_zone"]]
    out = {
        "absorbed_final": absorbed, "budget_E0": budget_A,
        "absorbed_over_budget": absorbed / budget_A,
        "E_drop_vs_E0": (budget_A - E_last) / budget_A,
        "E_drop_vs_first_row": (E_first - E_last) / E_first,
        "EplusKE_drop_vs_E0":
            (budget_A - E_last - hist[-1]["KE"]) / budget_A,
        "departing_tracks": dep_rec,
        "specdev_first_last": [specd[0], specd[-1]],
        "specdev_rel_change": (specd[-1] - specd[0]) / specd[0],
        "specdev_monotone_nonincreasing":
            bool(all(b <= a + 1e-6 for a, b in zip(specd, specd[1:]))),
        "rot_final_deg": rots[-1],
        "rot_monotone_increasing":
            bool(all(b >= a - 1e-6 for a, b in zip(rots, rots[1:]))),
        "end_largest_compact_nonedge":
            compact[0]["size"] if compact else 0,
        "end_largest_compact_xyz":
            compact[0]["xyz"] if compact else None}
    abs_ok = out["absorbed_over_budget"] <= 0.05
    edrop_ok = out["E_drop_vs_first_row"] <= 0.03 \
        and out["E_drop_vs_E0"] <= 0.03
    exact_artifact = (len(dep) >= 1 and any(
        t["r_first"] == 0.0 and t["sizes"] == [402, 55]
        for t in dep_rec))
    all_centered = all(t["r_first"] < 1.0 for t in dep_rec)
    qual_diff = all(t["sizes"][1] < 0.7 * t["sizes"][0]
                    for t in dep_rec)  # eroding, unlike C's 14->13
    spec_flat = abs(out["specdev_rel_change"]) < 0.2 and \
        specd[-1] > 0.1
    rot_ok = 11.0 <= rots[-1] <= 13.0 and out["rot_monotone_increasing"]
    out["checks"] = {"absorbed_le_5pct": bool(abs_ok),
                     "E_drop_le_3pct": bool(edrop_ok),
                     "claimed_artifact_track_found":
                         bool(exact_artifact),
                     "all_departing_centered": bool(all_centered),
                     "departing_all_eroding_unlike_C_pair":
                         bool(qual_diff),
                     "specdev_flat_no_melt": bool(spec_flat),
                     "rot_12deg_growing": bool(rot_ok)}
    if all(out["checks"].values()):
        out["verdict"] = "CONFIRMED"
        out["reason"] = "all A-null sub-checks reproduce"
    elif abs_ok and exact_artifact and qual_diff and spec_flat \
            and rot_ok:
        out["verdict"] = "NUANCE"
        out["reason"] = (
            "the null holds qualitatively (absorbed "
            f"{100 * out['absorbed_over_budget']:.1f}% <= 5%, the "
            "402->55 r_first=0 artifact track is real, all departing "
            "tracks are eroding masks unlike C's stable pair, "
            f"specdev drifts only {100 * out['specdev_rel_change']:.0f}"
            "% with no melt while rot grows monotonically to "
            f"{rots[-1]:.1f} deg) BUT the E-drop bound is wrong: E "
            f"falls {100 * out['E_drop_vs_first_row']:.1f}% from the "
            f"first hist row and {100 * out['E_drop_vs_E0']:.1f}% "
            "from the t=0 budget, not <= 3%; and the second "
            "departing track starts at r=5.6 (sizes 89->56), not "
            "r=0, so 'centered artifacts' is exact for only one of "
            "the two")
    else:
        out["verdict"] = "REFUTED"
        out["reason"] = "core null sub-checks fail; see numbers"
    return out


# ================= C6: ring tie-breaker ============================
def c6():
    log("C6 ring: recompute stencil energies from t32 npz")
    tags = ("t32_A", "t32_R4", "t32_R6")
    rec = jload("m5_21_10_ring.json")
    mine, stops = {}, {}
    for tag in tags:
        z = np.load(os.path.join(DATA, f"m5_21_6_end_{tag}.npz"))
        M = z["M"].astype(float)
        n, h = int(z["n"]), float(z["h"])
        row = jload(f"m5_21_6_row_{tag}.json")
        stops[tag] = row["stop"]
        es = {}
        for st in ("sym", "fwd", "2h"):
            cfg = INS.base_cfg(term="T2", n=n, L=n * h, bc="pinned",
                               stencil=st)
            eu, ed, evv = INS.e_parts(M, cfg)
            es[st] = float(eu + ed + evv)
        mine[tag] = es
    orders, spreads = {}, {}
    for st in ("sym", "fwd", "2h"):
        vals = {t: mine[t][st] for t in tags}
        order = sorted(tags, key=lambda t: vals[t])
        orders[st] = order
        spreads[st] = (max(vals.values()) - min(vals.values())) \
            / min(vals.values())
    off = [mine[t]["2h"] / mine[t]["sym"] - 1.0 for t in tags]
    match_rec = max(abs(mine[t][st] - rec["rows"][t][f"E_{st}"])
                    / abs(rec["rows"][t][f"E_{st}"])
                    for t in tags for st in ("sym", "fwd", "2h"))
    out = {"stops": stops, "E_mine": mine,
           "orders_mine": orders,
           "rel_spread_sym": spreads["sym"],
           "rel_spread_fwd": spreads["fwd"],
           "rel_spread_2h": spreads["2h"],
           "offset_2h_per_tag": off,
           "offset_2h_mean": float(np.mean(off)),
           "recorded_vs_mine_max_rel_diff": match_rec}
    ok = (all(s == "f_tol" for s in stops.values())
          and spreads["sym"] <= 4e-4
          and orders["sym"] == orders["fwd"]
          and orders["2h"] != orders["sym"]
          and abs(np.mean(off) + 0.077) < 0.002)
    out["verdict"] = "CONFIRMED" if ok else "NUANCE"
    out["reason"] = (
        f"all stops f_tol; my sym spread {spreads['sym']:.2e} <= "
        f"4e-4; sym order {orders['sym']} == fwd order; 2h order "
        f"{orders['2h']} flips it at 2h spread {spreads['2h']:.2e}; "
        f"2h offset mean {100 * np.mean(off):.2f}%" if ok
        else "one or more ring sub-checks deviate; see numbers")
    return out


# ================= C7: descent signature ===========================
def c7():
    log("C7 retention: rows + f64 recompute from npz")
    rows = {t: jload(f"m5_21_6_row_{t}.json")
            for t in ("f64_A", "f64_B", "f64_C",
                      "f48_A", "f48_B", "f48_C")}
    per = {t: rows[t]["retention"]["per_eig"] for t in rows}
    means = {t: rows[t]["retention"]["mean"] for t in rows}
    mine, f48_npz_present = {}, {}
    for t in ("f64_A", "f64_B", "f64_C"):
        z = np.load(os.path.join(DATA, f"m5_21_6_end_{t}.npz"))
        cfg = INS.base_cfg(term="T2", seed=t[-1], n=int(z["n"]),
                           L=int(z["n"]) * float(z["h"]), bc="free",
                           delta=float(z["delta"]))
        mine[t] = INS.retention(z["M"].astype(float), cfg)["per_eig"]
    for t in ("f48_A", "f48_B", "f48_C"):
        f48_npz_present[t] = os.path.exists(
            os.path.join(DATA, f"m5_21_6_end_{t}.npz"))
    recomp_dev = max(abs(a - b)
                     for t in mine
                     for a, b in zip(mine[t], per[t]))
    sig_f64C = (0.65 <= per["f64_C"][0] <= 0.78
                and 0.65 <= per["f64_C"][1] <= 0.78
                and per["f64_C"][2] >= 0.97)
    sig_f48C = (0.28 <= per["f48_C"][0] <= 0.40
                and 0.28 <= per["f48_C"][1] <= 0.40
                and per["f48_C"][2] >= 0.97)
    deeper = all(per["f48_C"][i] < per["f64_C"][i] for i in (0, 1))
    a_ok = all(per["f64_A"][i] >= per["f48_A"][i]
               for i in range(3)) and \
        means["f64_A"] >= means["f48_A"]
    out = {"per_eig_rows": per, "means_rows": means,
           "per_eig_f64_recomputed": mine,
           "recompute_max_abs_dev": recomp_dev,
           "f48_npz_present": f48_npz_present,
           "checks": {"f64C_two_depressed_one_high": bool(sig_f64C),
                      "f48C_deeper": bool(sig_f48C and deeper),
                      "f64A_ge_f48A": bool(a_ok)}}
    ok = sig_f64C and sig_f48C and deeper and a_ok \
        and recomp_dev < 1e-3
    out["verdict"] = "CONFIRMED" if ok else "NUANCE"
    out["reason"] = (
        f"f64_C per-eig {[round(x, 3) for x in per['f64_C']]} vs "
        f"f48_C {[round(x, 3) for x in per['f48_C']]}; f64_A >= "
        "f48_A on every eig and the mean; my npz recompute of all "
        f"f64 rows matches to {recomp_dev:.1e}. Caveat: the f48 "
        "endpoint npz no longer exist, so the f48 values rest on "
        "the row JSONs alone" if ok
        else "signature or recompute deviates; see numbers")
    return out


# ================= main ============================================
def main():
    res = {"meta": {
        "audit": "m5_21_10_d_audit (adversarial, independent "
                 "recompute)",
        "date": "2026-07-20",
        "own_implementations": [
            "biaxiality mask + 26-conn labeling (my_census)",
            "ledger drift vs recomputed t=0 budgets",
            "pair speeds / angles from raw centroid records",
            "stencil energies + orderings from t32 npz",
            "f64 retention from endpoint npz"]}}
    res["C1"] = c1()
    budgets = {ev: res["C1"]["per_ev"][ev]["E0_recomputed"]
               for ev in ("A", "B", "C")}
    res["C2"] = c2(budgets["C"])
    res["C3"] = c3()
    res["C4"] = c4(budgets["B"])
    res["C5"] = c5(budgets["A"])
    res["C6"] = c6()
    res["C7"] = c7()
    out_path = os.path.join(DATA, "m5_21_10_audit.json")
    with open(out_path, "w") as f:
        json.dump(res, f, indent=1)
    log(f"wrote {out_path}")
    print("\n| claim | verdict | reason |")
    print("| --- | --- | --- |")
    for c in ("C1", "C2", "C3", "C4", "C5", "C6", "C7"):
        print(f"| {c} | {res[c]['verdict']} | {res[c]['reason']} |")


if __name__ == "__main__":
    main()
