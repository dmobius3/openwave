"""M5.32 R27-3: the stored-field reads the author asked for after R26 (the tube-masked moment,
the director's charge inside the carrier lines, the carriers' distance against the box, the
Coulomb calibration on a uniaxial exterior) and the delta arm continued.

EQUATIONS FIRST
---------------
(a) THE TUBE-MASKED MOMENT (the author's item 5, 2026-09-25 19:04 UTC). On each shell of radius
r the pair-gap deviation dev = (lambda_1 - lambda_0) - delta is fitted by least squares to
    dev(r, theta, phi) = a_0 + a_1 cos(theta) + a_2 P_2(cos theta) + b_1 sin(theta) cos(phi)
                         + c_1 sin(theta) sin(phi)
on the cells OUTSIDE tubes of radius r_t around each carrier line (the carriers located per
sphere by the R26-0 partition reader; the tube = every cell within r_t of the line from the
origin through the carrier). A gap-zero line piercing a sphere removes a fixed AREA, so an
unmasked l = 1 read of a z-asymmetric line pattern falls as r^-2 by geometry (the author's
footprint); the masked read removes it, a field dipole survives the mask. Validated in R27-0
(e) on a synthetic dipole and a synthetic footprint. The m = 0 content = a_1^2 / (a_1^2 + b_1^2
+ c_1^2). The footprint prediction: a_1(footprint) = -(delta / (4 pi)) sum_lines
Omega_line cos(theta_line) x 3 with Omega_line = pi r_t^2 / r^2 the solid angle of one tube
(the l = 1 projection of a zero patch of solid angle Omega at polar angle theta_0 on a field
of value -delta: a_1 = 3 (-delta) Omega cos(theta_0) / (4 pi)).
(b) THE CHARGE IN THE LINES (item 6b). The director's topological charge density rho (R26-3
charge_density, the hedgehog Jacobian divergence) split into the carrier tubes (radius 3 h
around each carrier line) and the rest: the fraction of |rho| inside, <r^2> of each part.
(c) THE CARRIERS' DISTANCE FROM THE AXIS against the box (18:51 UTC item 2): rho_from_axis of
each carrier on r 9 (and r 18 where readable) on the stored R25-2 rows and the census rows;
the capped rows labeled FALLING.
(d) THE COULOMB CALIBRATION (18:51 UTC item 4). E(> R) = sum of the certified energy density
    e = 4 h^3 sum_{i<j} <F_ij, F_ij>_eta + V4
over the free cells with r > R; a Coulomb exterior e ~ C / r^4 gives E(> R) = 4 pi C / R, so
c = 4 pi C is read from a fit of E(> R) to c (1 / R - 1 / R_max) on R 9 to 18 with the local
exponent d ln E(> R) / d ln R beside it (-1 for Coulomb, +1 for the biaxial halo's 1 / R^2
density). No femtometres: the charge and mass normalization is the author's.
(e) THE DELTA ARM (item 6a): the S1 rows at delta 0.1 and 0.03 continued from their R26-3
stage fields (8000 iterations) for 8000 more under R25-2's descent (n 32, L 48, cap 16000,
gate 1e-4), the kick if at the gate, C_int / C_orb / C_rigid at the end with the sensitivity
to the 8000-iteration field.

PRE-REGISTERED LABELS (the record's R27 PLANNING section)
--------------------------------------------------------
    MOMENT_FOOTPRINT         on {1,1,1,1}: the masked l = 1 coefficient below three times the
                             masked control's floor on r 6 to 18, OR its slope outside -2 +- 0.5,
                             OR the unmasked coefficient within 30 percent of the footprint
                             prediction
    MOMENT_FIELD_CANDIDATE   the masked coefficient keeps the slope -2 +- 0.5 on r 6 to 18,
                             exceeds three times the floor on every shell, m = 0 content > 2/3
    MOMENT_UNREAD            otherwise
    the delta arm: R25-2's gate labels per row; C_rigid quoted at the gate only.
    (b), (c), (d) are reported, never labeled.

Modes: run_arm [workers] | reads | collect | smoke. Output: data/m5_32_r27_3_reads.json (the
reads), data/m5_32_r27_3_arm.json (the arm rows), arrays in data/m5_32_r27_3/ (local, kept).
Regenerate: run_arm 2 about 2 to 3 h; reads about 10 min; smoke a minute.
"""

import importlib.util
import json
import multiprocessing as mp
import os
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_JSON = os.path.join(DATA, "m5_32_r27_3_reads.json")
ARM_JSON = os.path.join(DATA, "m5_32_r27_3_arm.json")
OUT_NPZ = os.path.join(DATA, "m5_32_r27_3")
R26_3_DIR = os.path.join(DATA, "m5_32_r26_3")
R26_3_JSON = os.path.join(DATA, "m5_32_r26_3_spin.json")
R26_4_DIR = os.path.join(DATA, "m5_32_r26_4")
R26_4_JSON = os.path.join(DATA, "m5_32_r26_4_census.json")
R25_2_DIR = os.path.join(DATA, "m5_32_r25_2")
R25_2_JSON = os.path.join(DATA, "m5_32_r25_2_charge.json")
T0 = time.time()
G = 8.0
W1S = 25.0
ARM_CAP = 16000
TUBE_R = (1.5, 3.0)  # in units of h: the mask radii of (a)
CHARGE_TUBE_H = 3.0
FIT = (6.0, 18.0)
FOOT_TOL = 0.30
FLOOR_FACTOR = 3.0


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


S3 = _load("m5_32_r26_3_spin", "m5_32_r26_3_spin.py")
F0, R25 = S3.F0, S3.R25
CS, R21, R20, B3, R0, W1 = R25.CS, R25.R21, R25.R20, R25.B3, R25.R0, R25.W1
# the delta arm's arrays land here (R25's run_job writes to its module's OUT_NPZ)
R25.OUT_NPZ = OUT_NPZ
R25.OUT_JSON = ARM_JSON.replace(".json", "_rows.json")


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


def cfg_pot(n, L, delta, w1s=W1S):
    return S3.cfg_pot(n, L, delta, w1s)


# ============================================================================
# (e) THE DELTA ARM
# ============================================================================
def _load_retry(f, tries=5):
    """RUN-TIME ADDITION (14:45 UTC): the first arm pool died on a transient CRC error reading a
    file it had just written (the file loads cleanly afterwards); retry with a pause."""
    for k in range(tries):
        try:
            return np.load(f, allow_pickle=True)
        except Exception:  # noqa: BLE001
            if k == tries - 1:
                raise
            time.sleep(3.0)


def jobs_arm():
    return [dict(kind="S1", delta=d, w1s=W1S, n=32, L=48.0, cap=ARM_CAP) for d in (0.1, 0.03)]


def run_arm_job(j):
    """continue from the R26-3 stage field (copied once into this rung's folder)."""
    tag = R25.job_tag(j)
    os.makedirs(OUT_NPZ, exist_ok=True)
    stage = os.path.join(OUT_NPZ, tag + "_stage.npz")
    src = os.path.join(R26_3_DIR, tag + "_stage.npz")
    if not os.path.exists(stage):
        shutil.copy(src, stage)
    Z0 = np.load(src, allow_pickle=True)
    done0 = int(Z0["done"])
    row = R25.run_job(j)
    row["continued_from_iters"] = done0
    if row.get("status") == "OK":
        cfg, p, pot = cfg_pot(j["n"], j["L"], j["delta"])
        M = _load_retry(os.path.join(OUT_NPZ, row["tag"] + ".npz"))["M"]
        row["spin_reads"] = S3.spin_reads(M, cfg, p, pot, j["delta"], row["tag"])
        M0 = Z0["M"]
        row["spin_reads_at_8000"] = S3.spin_reads(
            M0, cfg, p, pot, j["delta"], row["tag"] + "@8000"
        )
        g1, g0 = row["spin_reads"]["generator"], row["spin_reads_at_8000"]["generator"]
        row["stage_to_end_sensitivity"] = {
            k: abs(g1[k] - g0[k]) / max(abs(g0[k]), 1e-300)
            for k in ("internal", "orbital", "rigid")
        }
        row["own_reads"] = {
            "biaxiality": F0.biaxiality_reads(M, cfg, j["delta"]),
            "gap_tail": F0.gap_tail(M, cfg, j["delta"]),
            "partition_r9": F0.partition_reader(M, cfg, 9.0, j["delta"])["partition"],
            "partition_r18": F0.partition_reader(M, cfg, 18.0, j["delta"])["partition"],
        }
    return row


def load_arm():
    if os.path.exists(ARM_JSON):
        with open(ARM_JSON) as f:
            return json.load(f)
    return {"task": "M5.32 R27-3 (e) the delta arm", "rows": {}, "collect": {}}


def save_arm(J):
    tmp = ARM_JSON + ".tmp"
    with open(tmp, "w") as f:
        json.dump(J, f, indent=1, default=str)
    os.replace(tmp, ARM_JSON)


def run_arm(workers):
    os.makedirs(OUT_NPZ, exist_ok=True)
    J = load_arm()
    pending = [j for j in jobs_arm() if J["rows"].get(R25.job_tag(j), {}).get("status") != "OK"]
    log(f"arm pool: {len(pending)} jobs, {workers} workers")
    with ProcessPoolExecutor(max_workers=int(workers), mp_context=mp.get_context("spawn")) as ex:
        futs = [ex.submit(run_arm_job, j) for j in pending]
        for fut in as_completed(futs):
            row = fut.result()
            Jn = load_arm()
            Jn["rows"][row["tag"]] = row
            save_arm(Jn)
    log("pool done")


def collect_arm():
    J = load_arm()
    prev = {}
    if os.path.exists(R26_3_JSON):
        with open(R26_3_JSON) as f:
            prev = json.load(f).get("arm_rows", {})
    table = {}
    for tag, r in J["rows"].items():
        if r.get("status") != "OK":
            table[tag] = {"status": r.get("status"), "stop": r.get("stop")}
            continue
        g = r["spin_reads"]["generator"]
        table[tag] = {
            "delta": r["delta"],
            "iters": r["iters"],
            "gate_label": r["gate_label"],
            "kick_label": r["kick_label"],
            "E": r["E"],
            "E_at_8000": (prev.get(tag) or {}).get("E"),
            "fmax_end": r["chunks"][-1]["fmax_spatial"],
            "drop_end": r["chunks"][-1]["drop"],
            "C_int": g["internal"],
            "C_orb": g["orbital"],
            "C_rigid": g["rigid"],
            "C_rigid_at_8000": r["spin_reads_at_8000"]["generator"]["rigid"],
            "stage_to_end_sensitivity": r["stage_to_end_sensitivity"],
            "C_rigid_quotable": r["gate_label"] == "AT_GATE",
            "partition_r9": r["own_reads"]["partition_r9"],
        }
    # the sign of C_rigid against delta across 0.03, 0.1, 0.3 (the 0.3 row from R25-2 via R26-3)
    c03 = {}
    if os.path.exists(R26_3_JSON):
        with open(R26_3_JSON) as f:
            lad = json.load(f).get("ladder", {})
        r03 = lad.get("S1_d0.3_w25_n32_L48") or lad.get("h1.5_L48") or {}
        c03 = r03.get("end", r03) if isinstance(r03, dict) else {}
    J["collect"] = {
        "rows": table,
        "delta_0.3_reference": c03,
        "all_three_at_gate": all(
            t.get("gate_label") == "AT_GATE" for t in table.values() if "gate_label" in t
        )
        and len(table) == 2,
        "note": "C_rigid quoted at the gate only; the sign against delta reported, never labeled unless all three rows are at the gate",
    }
    save_arm(J)
    return J["collect"]


# ============================================================================
# the reads (a) to (d) on the stored fields
# ============================================================================
F = _load("m5_32_r27_0_form", "m5_32_r27_0_form.py")
CENSUS_TAGS = {p: f"P{p}_d0.3_w25_n32_L48" for p in ("1_1_1_1", "2_1_1", "2_2", "3_1", "4")}
R25_2_ROWS = [
    ("S1_d0.3_w25_n32_L48", 32, 48.0),
    ("S1_d0.3_w25_n48_L72", 48, 72.0),
    ("S1_d0.3_w25_n64_L96", 64, 96.0),
    ("S1_d0.3_w25_n48_L48", 48, 48.0),
    ("S1_d0.3_w25_n64_L64", 64, 64.0),
]
SPHERES = (3.0, 4.5, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0)


def carriers_by_sphere(M, cfg, delta):
    """the carriers (theta, phi, half-units, rho) on every readable sphere."""
    out = {}
    for R in SPHERES:
        if R > 0.5 * cfg["L"] - 1.6:
            continue
        rd = F0.partition_reader(M, cfg, R, delta)
        out[R] = {
            "partition": rd["partition"],
            "carriers": [
                (c["theta"], c["phi"], c["half_units"], c["rho_from_axis"], c["z"])
                for c in rd["carriers"]
            ],
            "in_pin": rd["in_pin"],
        }
    return out


def ball_mask(cfg, points, r_t):
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = B3.coords(n, h)
    mask = np.zeros((n, n, n), dtype=bool)
    for q in points:
        mask |= (X - q[0]) ** 2 + (Y - q[1]) ** 2 + (Z - q[2]) ** 2 < r_t**2
    return mask


def moment_read(M, cfg, delta, r_t_h, carriers):
    """(a): the masked and unmasked l = 1 reads with the footprint from the masked-out content."""
    n, h = cfg["n"], cfg["h"]
    r_t = r_t_h * h
    lam = np.linalg.eigvalsh(M[..., 1:, 1:])
    dev = (lam[..., 1] - lam[..., 0]) - delta
    # the straight lines through the carriers read at r 9 (the plan's mask) ...
    c9 = carriers.get(9.0, {}).get("carriers", [])
    dirs = F.line_dirs([(th, ph) for th, ph, hu, rho, z in c9])
    line_mask = F.tube_mask(cfg, dirs, r_t) if dirs else np.zeros((n, n, n), dtype=bool)
    # ... and the balls around every carrier point read on the spheres (the lines curve)
    pts = []
    for R, rec in carriers.items():
        if rec["in_pin"]:
            continue
        for th, ph, hu, rho, z in rec["carriers"]:
            pts.append((R * np.sin(th) * np.cos(ph), R * np.sin(th) * np.sin(ph), R * np.cos(th)))
    balls = ball_mask(cfg, pts, r_t) | (
        line_mask & (np.sqrt(sum(a * a for a in B3.coords(n, h))) < 6.0)
    )
    out = {
        "r_t": r_t,
        "lines_from_r9": [[np.degrees(th), np.degrees(ph), hu] for th, ph, hu, rho, z in c9],
    }
    un = F.masked_harmonics(dev, cfg, np.zeros((n, n, n), dtype=bool))
    out["unmasked"] = un
    for name, mask in (("line", line_mask), ("balls", balls)):
        rd = F.masked_harmonics(dev, cfg, mask)
        # the footprint: the l = 1 projection of the masked-out content on each shell
        X, Y, Z = B3.coords(n, h)
        r = np.sqrt(X * X + Y * Y + Z * Z)
        cth = Z / np.maximum(r, 1e-12)
        pin = B3.pin_shell(n, h)
        foot = []
        for q in rd["shells"]:
            sh = (np.abs(r - q["r"]) < 0.75) & ~pin & mask
            dOmega = h**3 / (q["r"] ** 2 * 1.5)
            foot.append(float(3.0 / (4.0 * np.pi) * np.sum(dev[sh] * cth[sh]) * dOmega))
        rd["footprint_a1_from_masked_content"] = foot
        rd["unmasked_a1"] = [q["a1"] for q in un["shells"]][: len(foot)]
        rd["masked_plus_footprint_over_unmasked"] = [
            (q["a1"] + f_) / u if abs(u) > 1e-300 else None
            for q, f_, u in zip(rd["shells"], foot, rd["unmasked_a1"])
        ]
        out[name] = rd
    return out


def moment_label(rd_balls, rd_un, fit=(6.0, 18.0)):
    """the pre-registered rule on {1,1,1,1}: FOOTPRINT / FIELD_CANDIDATE / UNREAD."""
    win = [q for q in rd_balls["shells"] if fit[0] <= q["r"] <= fit[1]]
    winu = [q for q in rd_un["shells"] if fit[0] <= q["r"] <= fit[1]]
    if len(win) < 3:
        return "MOMENT_UNREAD", {"reason": "fewer than three shells"}
    floor = [q["residual_rms"] / np.sqrt(max(q["cells"], 1)) for q in win]
    above = [abs(q["a1"]) > FLOOR_FACTOR * f_ for q, f_ in zip(win, floor)]
    slope = rd_balls["slope_a1"]
    m0 = [q["m0_content"] for q in win]
    foot = rd_balls["footprint_a1_from_masked_content"][: len(rd_balls["shells"])]
    footw = [f_ for q, f_ in zip(rd_balls["shells"], foot) if fit[0] <= q["r"] <= fit[1]]
    un_agrees = [
        abs(u["a1"] - f_) < FOOT_TOL * abs(u["a1"])
        for u, f_ in zip(winu, footw)
        if abs(u["a1"]) > 0
    ]
    info = {
        "masked_a1": [q["a1"] for q in win],
        "floor_per_shell": floor,
        "above_3x_floor_per_shell": above,
        "masked_slope": slope,
        "m0_content_per_shell": m0,
        "unmasked_a1": [q["a1"] for q in winu],
        "footprint_a1": footw,
        "unmasked_within_30pct_of_footprint_per_shell": un_agrees,
    }
    slope_in = slope is not None and abs(slope + 2.0) <= 0.5
    if (not all(above)) or (not slope_in) or (un_agrees and all(un_agrees)):
        return "MOMENT_FOOTPRINT", info
    if all(above) and slope_in and min(m0) > 2.0 / 3.0:
        return "MOMENT_FIELD_CANDIDATE", info
    return "MOMENT_UNREAD", info


def load_census_row(part):
    f = os.path.join(R26_4_DIR, CENSUS_TAGS[part] + ".npz")
    return np.load(f)["M"] if os.path.exists(f) else None


def reads_a():
    J = json.load(open(R26_4_JSON))["rows"]
    out = {}
    cfg, p, pot = cfg_pot(32, 48.0, 0.3)
    for part in ("1_1_1_1", "3_1", "2_2", "2_1_1", "4"):
        M = load_census_row(part)
        if M is None:
            continue
        rw = J.get(CENSUS_TAGS[part], {})
        car = carriers_by_sphere(M, cfg, 0.3)
        rec = {
            "gate_label": rw.get("gate_label"),
            "kick_label": rw.get("kick_label"),
            "certified": rw.get("gate_label") == "AT_GATE"
            and rw.get("kick_label") in ("STABLE", "SADDLE"),
            "carriers_by_sphere": {f"R{R:g}": v for R, v in car.items()},
        }
        for rth in TUBE_R:
            rec[f"tube_{rth:g}h"] = moment_read(M, cfg, 0.3, rth, car)
        lab, info = moment_label(rec["tube_3h"]["balls"], rec["tube_3h"]["unmasked"])
        rec["label_by_rule"] = lab
        rec["label_info"] = info
        out[part] = rec
        log(
            f"(a) {part}: {lab} masked slope {info.get('masked_slope')} unmasked slope {rec['tube_3h']['unmasked']['slope_a1']}"
        )
    out["label"] = out.get("1_1_1_1", {}).get("label_by_rule", "MOMENT_UNREAD")
    out["note"] = (
        "the label is {1,1,1,1}'s by the pre-registered rule on the 3 h balls mask; the other branches reported"
    )
    return out


def reads_b():
    out = {"floor_from_R27_0_f": None, "rows": {}}
    J0 = json.load(open(os.path.join(DATA, "m5_32_r27_0_form.json")))
    f0 = J0["f_charge_control"]
    out["floor_from_R27_0_f"] = {
        "net_outer": f0["net_charge_outer"],
        "abs_outer": f0["fraction_outer_in_tubes"] + f0["fraction_outer_outside_tubes"],
    }
    Jc = json.load(open(R26_4_JSON))["rows"]
    Jr = json.load(open(R25_2_JSON))["rows"]
    items = []
    for tag, n, L in R25_2_ROWS:
        f = os.path.join(R25_2_DIR, tag + ".npz")
        if os.path.exists(f):
            items.append((tag, n, L, np.load(f)["M"], Jr.get(tag, {})))
    for part in ("1_1_1_1", "3_1", "2_2", "2_1_1", "4"):
        M = load_census_row(part)
        if M is not None:
            items.append((CENSUS_TAGS[part], 32, 48.0, M, Jc.get(CENSUS_TAGS[part], {})))
    for tag, n, L, M, rw in items:
        cfg, p, pot = cfg_pot(n, L, 0.3)
        car = carriers_by_sphere(M, cfg, 0.3)
        c9 = car.get(9.0, {}).get("carriers", [])
        dirs = F.line_dirs([(th, ph) for th, ph, hu, rho, z in c9])
        sp = F.charge_split(M, cfg, dirs, CHARGE_TUBE_H * cfg["h"])
        sp["carriers_r9"] = [
            [np.degrees(th), np.degrees(ph), hu, rho] for th, ph, hu, rho, z in c9
        ]
        sp["gate_label"] = rw.get("gate_label")
        sp["kick_label"] = rw.get("kick_label")
        sp["n_L_h"] = [n, L, cfg["h"]]
        sp["tubes_minus_floor"] = sp["fraction_outer_in_tubes"] - out["floor_from_R27_0_f"][
            "abs_outer"
        ] * sp["tube_cells_outer"] / max(int((~B3.pin_shell(n, cfg["h"])).sum()), 1)
        out["rows"][tag] = sp
        log(
            f"(b) {tag}: in tubes {sp['fraction_outer_in_tubes']:.3f} rest {sp['fraction_outer_outside_tubes']:.3f} core {sp['fraction_inside_core']:.3f} net outer {sp['net_charge_outer']:.3f}"
        )
    return out


def reads_c():
    out = {"rows": {}}
    Jr = json.load(open(R25_2_JSON))["rows"]
    Jc = json.load(open(R26_4_JSON))["rows"]
    items = []
    for tag, n, L in R25_2_ROWS:
        f = os.path.join(R25_2_DIR, tag + ".npz")
        if os.path.exists(f):
            items.append((tag, n, L, np.load(f)["M"], Jr.get(tag, {})))
    for part in ("1_1_1_1", "3_1", "2_2", "2_1_1", "4"):
        M = load_census_row(part)
        if M is not None:
            items.append((CENSUS_TAGS[part], 32, 48.0, M, Jc.get(CENSUS_TAGS[part], {})))
    for tag, n, L, M, rw in items:
        cfg, p, pot = cfg_pot(n, L, 0.3)
        car = carriers_by_sphere(M, cfg, 0.3)
        out["rows"][tag] = {
            "n_L_h": [n, L, cfg["h"]],
            "gate_label": rw.get("gate_label"),
            "kick_label": rw.get("kick_label"),
            "capped": rw.get("gate_label") != "AT_GATE",
            "carriers_by_sphere": {
                f"R{R:g}": {
                    "partition": v["partition"],
                    "rho_hu": [[rho, hu, z] for th, ph, hu, rho, z in v["carriers"]],
                }
                for R, v in car.items()
            },
        }
    out["note"] = (
        "the carriers' distance from the polar axis (rho) and their z per sphere; capped rows labeled; the box ladder is the S1 rows"
    )
    return out


def reads_d():
    out = {}
    J0 = json.load(open(os.path.join(DATA, "m5_32_r27_0_form.json")))["g_coulomb_uniaxial"]
    out["uniaxial_hedgehog"] = {
        k: J0[k]
        for k in (
            "C_sympy_pole",
            "c_4piC",
            "c_fit",
            "R_eff_fit",
            "c_fit_over_4piC",
            "local_exponent_E_gt_R",
            "E_total_free",
            "c_over_E_total",
            "shells_r_C_lattice",
        )
    }
    out["convention"] = (
        "e = 4 h^3 sum_{i<j} <F_ij, F_ij>_eta + V4 on M = n n^T (delta 0), E(> R) = 4 pi C / R = c / R, c = 32 pi; c / E in code lengths, no femtometres"
    )
    Jr = json.load(open(R25_2_JSON))["rows"]
    out["biaxial_rows"] = {}
    for tag, n, L in R25_2_ROWS:
        f = os.path.join(R25_2_DIR, tag + ".npz")
        if not os.path.exists(f):
            continue
        M = np.load(f)["M"]
        cfg, p, pot = cfg_pot(n, L, 0.3)
        e = R20.density(M, cfg, pot)
        h = cfg["h"]
        X, Y, Z = B3.coords(n, h)
        r = np.sqrt(X * X + Y * Y + Z * Z)
        pin = B3.pin_shell(n, h)
        rmax = 0.5 * L - 1.6
        Rs = np.arange(6.0, rmax, 1.5)
        Eout = np.array([float(e[(r > R_) & ~pin].sum()) for R_ in Rs])
        sel = (Rs >= 9.0) & (Rs <= min(18.0, rmax - 3.0))
        lg = np.polyfit(np.log(Rs[sel]), np.log(np.maximum(Eout[sel], 1e-300)), 1)
        lin = np.polyfit(Rs[sel], Eout[sel], 1)
        # a Coulomb exterior gives ONE c on every window: c(R) = E(> R) / (1 / R - 1 / rmax)
        c_of_R = [
            [float(R_), float(E_ / (1.0 / R_ - 1.0 / rmax))]
            for R_, E_ in zip(Rs, Eout)
            if R_ < rmax - 1.5
        ]
        # the shell energy density's own exponent: e(r) on shells |r - R_s| < 0.75, per unit volume
        sh = []
        for R_ in Rs:
            m = (np.abs(r - R_) < 0.75) & ~pin
            if m.sum() > 20:
                sh.append([float(R_), float(np.mean(e[m]) / h**3)])
        sha = np.array([q for q in sh if 9.0 <= q[0] <= min(18.0, rmax - 3.0)])
        pd = (
            np.polyfit(np.log(sha[:, 0]), np.log(np.maximum(sha[:, 1], 1e-300)), 1)
            if len(sha) >= 3
            else [None]
        )
        out["biaxial_rows"][tag] = {
            "n_L_h": [n, L, h],
            "gate_label": Jr.get(tag, {}).get("gate_label"),
            "E_total_free": float(e[~pin].sum()),
            "E_gt_R": [[float(a), float(b)] for a, b in zip(Rs, Eout)],
            "local_exponent_E_gt_R_on_9_18": float(lg[0]),
            "linear_slope_dE_dR_on_9_18": float(lin[0]),
            "c_of_R_if_coulomb": c_of_R,
            "c_drift_9_to_18": (
                [q[1] for q in c_of_R if q[0] == 9.0][0]
                and [q[1] for q in c_of_R if q[0] == 18.0][0]
                / [q[1] for q in c_of_R if q[0] == 9.0][0]
                if any(q[0] == 18.0 for q in c_of_R) and any(q[0] == 9.0 for q in c_of_R)
                else None
            ),
            "shell_density_r_e": sh,
            "shell_density_exponent_on_9_18": float(pd[0]) if pd[0] is not None else None,
            "note": "a Coulomb exterior has e ~ r^-4 and one c; the R23 halo has e ~ r^-2",
        }
    out["note"] = (
        "the biaxial exterior's E(> R) is linear in R (a 1 / R^2 halo density, R23), not c / R; c is read on the uniaxial hedgehog only"
    )
    return out


def run_reads():
    out = {"task": "M5.32 R27-3 the reads"}
    only = os.environ.get("R27_3_ONLY")
    prev = json.load(open(OUT_JSON)) if (only and os.path.exists(OUT_JSON)) else {}
    for k, fn in (
        ("a_moment", reads_a),
        ("b_charge_split", reads_b),
        ("c_carriers_vs_box", reads_c),
        ("d_coulomb", reads_d),
    ):
        if only and k[0] not in only and k in prev:
            out[k] = prev[k]
            continue
        t = time.time()
        try:
            out[k] = fn()
        except Exception as e:  # noqa: BLE001
            import traceback

            out[k] = {"error": repr(e), "traceback": traceback.format_exc()}
        out[k]["wall_s"] = round(time.time() - t, 1)
        log(f"{k} done ({out[k]['wall_s']} s)")
    out["labels"] = {"moment": out.get("a_moment", {}).get("label")}
    out["wall_s"] = round(time.time() - T0, 1)
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=1, default=str)
    return out


def smoke():
    """the masked reader against R26-0's gap_tail on one stored field (the unmasked a1 must agree)."""
    cfg, p, pot = cfg_pot(32, 48.0, 0.3)
    M = load_census_row("1_1_1_1")
    lam = np.linalg.eigvalsh(M[..., 1:, 1:])
    dev = (lam[..., 1] - lam[..., 0]) - 0.3
    un = F.masked_harmonics(dev, cfg, np.zeros((32, 32, 32), dtype=bool))
    gt = F0.gap_tail(M, cfg, 0.3)
    a1_new = {q["r"]: q["a1"] for q in un["shells"]}
    a1_old = {row[0]: row[2] for row in gt["shells_r_a0_a1_a2_maxdev_n"]}
    diffs = [
        abs(a1_new[r] - a1_old[r]) / max(abs(a1_old[r]), 1e-300)
        for r in a1_old
        if r in a1_new and 6.0 <= r <= 18.0
    ]
    res = {"max_rel_diff_a1_vs_gap_tail": float(max(diffs)), "PASS": bool(max(diffs) < 0.05)}
    print(json.dumps(res, indent=1))
    return res


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "smoke"
    if mode == "run_arm":
        run_arm(sys.argv[2] if len(sys.argv) > 2 else 2)
        print(json.dumps(collect_arm(), indent=1, default=str))
    elif mode == "collect_arm":
        print(json.dumps(collect_arm(), indent=1, default=str))
    elif mode == "reads":
        o = run_reads()
        print(json.dumps(o["labels"], indent=1))
    elif mode == "smoke":
        smoke()
    else:
        raise SystemExit(f"unknown mode {mode}")


if __name__ == "__main__":
    main()
