"""M5.32 R0 (b): baseline reproduction driver (comparison only, no physics).

Runs the certified record scripts UNMODIFIED on this machine, with their
DATA / PLOTS module globals redirected to ../data/m5_32_r0_b/ so the
stored reference JSONs in ../data/ are never overwritten, then compares
every numeric leaf of each reproduced JSON against its reference.

Usage (from research/scripts/):
    python m5_32_r0_b_baseline.py run <item>      # one item, writes
                                                  # ../data/m5_32_r0_b/<item>.meta.json
    python m5_32_r0_b_baseline.py compare         # ../data/m5_32_r0_baseline.json
    python m5_32_r0_b_baseline.py list

Items (item -> record script -> reference JSON):
    s16a   m5_21_16_a_symbolic.py    m5_21_16_symbolic.json   (CHAN/INV/IDENT symbolic side)
    s16b   m5_21_16_b_field.py       m5_21_16_field.json      (IDENT / INV / CHAN / DRESS)
    s15f   m5_21_15_f_guard.py       m5_21_15_guard.json      (fixed-J guard ladder)
    s15g   m5_21_15_g_fomnarrow.py   m5_21_15_fom_narrow.json (omega* = 0.59, E = +115.9)
    s15c   m5_21_15_c_fixedj.py      m5_21_15_fixedj.json     (fixed-J bridge, full)
    p4     m5_21_4_a_pair.py ladder it=120 n=32
                                     m5_21_4_ladder_it120.json (like-charge repulsion, 3 d)
    p17    m5_17_two_charge.py fixed m5_17_two_charge_fixed.json (Coulomb 1/d fit, both signs)
    g11    m5_21_11_g_controls.py    m5_21_11_garm_controls.json (vacuum null, identity, E_min)
    m14c   m5_21_14_c_minimize.py    m5_21_14_minimize.json   (analytic dressed family MIN)
    m14a   m5_21_14_a_symbolic.py    m5_21_14_symbolic.json   (T1 kin = -8 sum |Mdot3 v_i|^2)

Gate: every numeric leaf within 1e-3 relative (abs floor 1e-12 for
near-zero leaves such as invariance drifts), booleans / strings equal.
Keys runtime_s / wall_s / note / reading are excluded from the gate.
"""
from __future__ import annotations

import importlib.util
import json
import os
import platform
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT = os.path.join(DATA, "m5_32_r0_b")
PY = sys.executable

REL_TOL = 1e-3
ABS_FLOOR = 1e-12
SKIP_KEYS = {"runtime_s", "wall_s", "note", "reading", "statement"}

ITEMS = {
    "s16a": dict(script="m5_21_16_a_symbolic.py", ref="m5_21_16_symbolic.json",
                 inputs=[], call="main"),
    "s16b": dict(script="m5_21_16_b_field.py", ref="m5_21_16_field.json",
                 inputs=[], call="main"),
    "s15f": dict(script="m5_21_15_f_guard.py", ref="m5_21_15_guard.json",
                 inputs=[], call="main"),
    "s15g": dict(script="m5_21_15_g_fomnarrow.py", ref="m5_21_15_fom_narrow.json",
                 inputs=["m5_21_15_guard.json"], call="main"),
    "s15c": dict(script="m5_21_15_c_fixedj.py", ref="m5_21_15_fixedj.json",
                 inputs=["m5_21_14_minimize.json"], call="main"),
    "p4": dict(script="m5_21_4_a_pair.py", ref="m5_21_4_ladder_it120.json",
               inputs=[], call="ladder", kwargs={"it": 120, "n": 32}),
    "p17": dict(script="m5_17_two_charge.py", ref="m5_17_two_charge_fixed.json",
                inputs=[], call="run_fixed", argv=["m5_17_two_charge.py", "fixed"]),
    "g11": dict(script="m5_21_11_g_controls.py", ref="m5_21_11_garm_controls.json",
                inputs=[], call="main"),
    "m14c": dict(script="m5_21_14_c_minimize.py", ref="m5_21_14_minimize.json",
                 inputs=["m5_21_14_verify.json"], call="main"),
    "m14a": dict(script="m5_21_14_a_symbolic.py", ref="m5_21_14_symbolic.json",
                 inputs=[], call="main"),
}

# headline leaves per item (path into the JSON), reported explicitly
HEADLINES = {
    "s16a": ["N1_baseline_lead", "N2_variantA_lead", "N4_equals_variantA",
             "B1_H_equals_sum_inner_eta", "B2_Hflip_equals_sum_frobenius", "all_pass"],
    "s16b": ["IDENT.u_eta", "IDENT.u_flip", "IDENT.rel_diff",
             "INV.eta.so3_rot", "INV.eta.so13_boost", "INV.eta.so4_compact",
             "INV.flip.so3_rot", "INV.flip.so13_boost", "INV.flip.so4_compact",
             "CHAN.rows.clock_local.kin_eta", "CHAN.rows.boost_z.kin_eta",
             "CHAN.rows.boost_z.kin_flip", "CHAN.rows.boost_x.kin_eta",
             "CHAN.rows.boost_x.kin_flip", "DRESS.kin_base_eta"],
    "s15f": ["rows[7].E_total", "rows[7].omega_star", "rows[7].bound",
             "rows[1].E_total", "rows[1].omega_star"],
    "s15g": ["omega_min", "E_min_total", "interior_minimum", "E_min_positive"],
    "s15c": ["kin_base_clock", "E_base_u", "EJ[0].omega_star", "EJ[0].E_total",
             "EJ[1].omega_star", "EJ[1].E_total", "EJ[2].omega_star", "EJ[2].E_total"],
    "p4": ["rows[0].E", "rows[1].E", "rows[2].E", "rows[3].E", "rows[4].E",
           "rows[5].E", "rows[3].charge.far[0]", "E_single", "c2_selfcal"],
    "p17": ["E_single_same_box", "curves.likepair.fit_A", "curves.likepair.fit_E0",
            "curves.likepair.sign_ok", "curves.antipair.fit_A",
            "curves.antipair.sign_ok", "curves.likepair.A_over_prediction"],
    "g11": ["C1_vacuum_null[0].E_min", "C1_vacuum_null[1].E_min",
            "C2_field_identity[0].field_diff_max", "C2_field_identity[1].field_diff_max",
            "C2_field_identity[2].field_diff_max",
            "C3_record_match[0].E_min_reproduced", "C3_record_match[1].E_min_reproduced",
            "C3_record_match[2].E_min_reproduced", "all_pass"],
    "m14c": ["verdicts.E_corr_at_bstar", "verdicts.kin_corr_at_bstar",
             "rigid.E_rigid", "rigid.b_const_best", "g_flatness.g32.E_corr",
             "g_flatness.g32.kin_corr", "lattice_crosscheck.E_corr_lattice",
             "lattice_crosscheck.kin_corr_lattice", "DIAG.t1_reached"],
    "m14a": ["T1.t1_kin_str", "T1.t1_kin_identity", "T1.t1_static_identity",
             "BR[0].t1_u_value", "BR[0].t1_k_value", "BR[1].t1_u_value",
             "BR[1].t1_k_value", "all_green"],
}


def _versions():
    import numpy
    import scipy
    import sympy
    return {"python": sys.version.split()[0], "numpy": numpy.__version__,
            "scipy": scipy.__version__, "sympy": sympy.__version__,
            "platform": platform.platform(), "machine": platform.machine()}


def _get(d, path):
    cur = d
    for part in path.split("."):
        if "[" in part:
            key, idx = part[:-1].split("[")
            if key:
                cur = cur[key]
            cur = cur[int(idx)]
        else:
            cur = cur[part]
    return cur


def run(item):
    spec = ITEMS[item]
    os.makedirs(OUT, exist_ok=True)
    for name in spec["inputs"]:
        dst = os.path.join(OUT, name)
        if not os.path.exists(dst):
            shutil.copy(os.path.join(DATA, name), dst)
    if "argv" in spec:
        sys.argv = list(spec["argv"])
    os.environ.setdefault("MPLBACKEND", "Agg")
    path = os.path.join(HERE, spec["script"])
    t0 = time.time()
    mspec = importlib.util.spec_from_file_location(f"r0b_{item}", path)
    mod = importlib.util.module_from_spec(mspec)
    mspec.loader.exec_module(mod)
    mod.DATA = OUT
    if hasattr(mod, "PLOTS"):
        mod.PLOTS = OUT
    getattr(mod, spec["call"])(**spec.get("kwargs", {}))
    dt = time.time() - t0
    meta = {"item": item, "script": spec["script"], "call": spec["call"],
            "kwargs": spec.get("kwargs", {}), "argv": spec.get("argv"),
            "runtime_s": round(dt, 1), "versions": _versions(),
            "command": f"{PY} m5_32_r0_b_baseline.py run {item}",
            "reproduced_json": os.path.join("m5_32_r0_b", spec["ref"])}
    with open(os.path.join(OUT, f"{item}.meta.json"), "w") as f:
        json.dump(meta, f, indent=1)
    print(json.dumps({"item": item, "runtime_s": meta["runtime_s"]}))


def _leaves(obj, prefix=""):
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k in SKIP_KEYS:
                continue
            yield from _leaves(v, f"{prefix}.{k}" if prefix else k)
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from _leaves(v, f"{prefix}[{i}]")
    else:
        yield prefix, obj


def _cmp_leaf(a, b):
    """returns (rel_diff or None, ok)."""
    if isinstance(a, bool) or isinstance(b, bool) or isinstance(a, str) \
            or isinstance(b, str) or a is None or b is None:
        return None, a == b
    try:
        fa, fb = float(a), float(b)
    except (TypeError, ValueError):
        return None, a == b
    if fa != fa and fb != fb:      # both nan
        return 0.0, True
    den = max(abs(fa), abs(fb))
    if den < ABS_FLOOR:
        return 0.0, True
    rel = abs(fa - fb) / den
    if abs(fa - fb) < ABS_FLOOR:
        return rel, True
    return rel, rel <= REL_TOL


def compare_item(item):
    spec = ITEMS[item]
    ref_path = os.path.join(DATA, spec["ref"])
    rep_path = os.path.join(OUT, spec["ref"])
    meta_path = os.path.join(OUT, f"{item}.meta.json")
    row = {"item": item, "script": spec["script"], "reference_json": spec["ref"],
           "command": f"{PY} m5_32_r0_b_baseline.py run {item}"}
    with open(ref_path) as f:
        ref = json.load(f)
    row["reference_runtime_s"] = ref.get("runtime_s", ref.get("wall_s"))
    if not os.path.exists(meta_path) or not os.path.exists(rep_path):
        partial = os.path.exists(rep_path)
        row.update({"verdict": "TIMEOUT" if partial else "NOT_RUN",
                    "partial_output": partial})
        if not partial:
            return row
        row["runtime_s"] = None
    else:
        with open(meta_path) as f:
            meta = json.load(f)
        row["runtime_s"] = meta["runtime_s"]
        row["versions"] = meta["versions"]
    with open(rep_path) as f:
        rep = json.load(f)
    ref_leaves = dict(_leaves(ref))
    rep_leaves = dict(_leaves(rep))
    n_ok = n_fail = n_missing = 0
    worst = (0.0, None)
    fails = []
    for k, va in ref_leaves.items():
        if k not in rep_leaves:
            n_missing += 1
            continue
        rel, ok = _cmp_leaf(va, rep_leaves[k])
        if rel is not None and rel > worst[0]:
            worst = (rel, k)
        if ok:
            n_ok += 1
        else:
            n_fail += 1
            if len(fails) < 40:
                fails.append({"key": k, "reference": va,
                              "reproduced": rep_leaves[k], "rel_diff": rel})
    heads = []
    for p in HEADLINES.get(item, []):
        try:
            a = _get(ref, p)
        except (KeyError, IndexError, TypeError):
            continue
        try:
            b = _get(rep, p)
        except (KeyError, IndexError, TypeError):
            b = None
        rel, ok = _cmp_leaf(a, b)
        heads.append({"key": p, "reference": a, "reproduced": b,
                      "rel_diff": rel, "ok": ok})
    row.update({"n_leaves_compared": n_ok + n_fail, "n_ok": n_ok,
                "n_fail": n_fail, "n_missing_in_repro": n_missing,
                "max_rel_diff": worst[0], "max_rel_diff_key": worst[1],
                "headlines": heads, "failures": fails})
    if row.get("verdict") == "TIMEOUT":
        return row
    row["verdict"] = "PASS" if (n_fail == 0 and n_missing == 0) else "FAIL"
    return row


def compare():
    out = {"task": "M5.32 R0 (b) baseline reproduction",
           "gate": {"rel_tol": REL_TOL, "abs_floor": ABS_FLOOR,
                    "skipped_keys": sorted(SKIP_KEYS)},
           "versions": _versions(), "items": {}}
    for item in ITEMS:
        out["items"][item] = compare_item(item)
    out["all_pass"] = all(r["verdict"] == "PASS" for r in out["items"].values())
    with open(os.path.join(DATA, "m5_32_r0_baseline.json"), "w") as f:
        json.dump(out, f, indent=1)
    for item, r in out["items"].items():
        print(f"{item:5s} {r['verdict']:8s} leaves={r.get('n_leaves_compared')} "
              f"fail={r.get('n_fail')} maxrel={r.get('max_rel_diff')} "
              f"({r.get('max_rel_diff_key')}) rt={r.get('runtime_s')}s "
              f"ref_rt={r.get('reference_runtime_s')}s")
    print("all_pass", out["all_pass"])


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "list"
    if cmd == "run":
        run(sys.argv[2])
    elif cmd == "compare":
        compare()
    else:
        print("\n".join(ITEMS))
