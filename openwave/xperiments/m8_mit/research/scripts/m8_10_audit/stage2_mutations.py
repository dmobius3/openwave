"""Stage 2, tasks 2 and 3.
(a) Reproduction: solver_rerun/results.json (from its run_all.sh) versus its shipped results.json.
(b) Mutation tests on copies of its code (fresh folder per mutation under solver_mut/), never touching solver_work/.
    M1 item 9 residual: scale two levels of N inside its engine -> does its Casimir residual notice?
    M2 identify.py: corrupt one 110-digit engine number -> does its pipeline fail?
    M3 quadcheck: break DN in the float quadrature -> does anything fail?
    M4 item5 controls: point a 'control' at a zero level -> does anything fail?
Runs are sequential, each through ./py in the copy."""
import json
import os
import re
import shutil
import subprocess

import lib

R = lib.ROOM
SRC = R / "solver_work"
out = {}


def fresh(name):
    d = R / "solver_mut" / name
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    for f in SRC.iterdir():
        if f.suffix in (".py", ".json", ".sh"):
            shutil.copy(f, d / f.name)
    shutil.copy(R / "py", d / "py")
    return d


def run(d, args, env=None):
    e = dict(os.environ)
    if env:
        e.update(env)
    p = subprocess.run(["./py"] + args, cwd=d, capture_output=True, text=True, env=e)
    return p.returncode, p.stdout, p.stderr


def patch(path, old, new, count=1):
    s = path.read_text()
    assert s.count(old) == count, (path.name, old[:50], s.count(old))
    path.write_text(s.replace(old, new))


# ---------------- (a) reproduction
orig = json.load(open(R / "solver_rerun" / "results_original.json"))
new = json.load(open(R / "solver_rerun" / "results.json"))


def diff(a, b, path=""):
    d = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in set(a) | set(b):
            if k not in a or k not in b:
                d.append(path + "/" + str(k) + " (missing on one side)")
            else:
                d += diff(a[k], b[k], path + "/" + str(k))
    elif isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        for i, (x, y) in enumerate(zip(a, b)):
            d += diff(x, y, path + f"[{i}]")
    elif a != b:
        d.append(f"{path}: {a!r} -> {b!r}")
    return d


dd = diff(orig, new)
out["reproduction"] = {"identical": len(dd) == 0, "n_differences": len(dd), "differences": dd[:50]}
print("reproduction: identical" if not dd else f"reproduction: {len(dd)} differences, first: {dd[:5]}")

# ---------------- M1: item 9 residual blind to errors in N
d = fresh("M1")
eng = d / "engine.py"
patch(eng, '        N = {J: product(Sf, Phi, J, "sv") for J in range(0, 10)}\n',
      '        N = {J: product(Sf, Phi, J, "sv") for J in range(0, 10)}\n'
      '        if os.environ.get("MUT") == "1":\n'
      '            N[5] = N[5] * mp.mpf("1.37"); N[7] = N[7] * mp.mpf("0.5"); N[9] = N[9] * 0\n')
patch(eng, "import sys\n", "import sys\nimport os\n")
res = {}
for tag, env in (("base", {"MUT": "0"}), ("mutant", {"MUT": "1"})):
    rc, so, se = run(d, ["engine.py", "30"], env)
    lines = [l for l in so.splitlines() if "casimir resid=" in l]
    resid = [float(re.search(r"casimir resid=(\S+)", l).group(1)) for l in lines]
    lam = [float(re.search(r"DN along=\((\S+?)[ ,]", l).group(1)) if "DN along=(" in l
           else float(re.search(r"DN along=(\S+)", l).group(1)) for l in lines]
    res[tag] = {"exit": rc, "max_casimir_resid": max(resid) if resid else None, "lambda4_first_rays": lam[:10],
                "stderr_tail": se[-300:]}
    print(f"M1 {tag}: exit {rc}, max Casimir residual {max(resid) if resid else None}, lambda4 {lam[:3]}")
out["M1_item9_residual_with_corrupted_N"] = res

# ---------------- M2: identify.py with one corrupted engine value
d = fresh("M2")
e110 = json.load(open(d / "out_engine_110.json"))
key = "3:R1"
x = e110["rays"][key]["xi_level_norm2_over_g2"]["10"]
import mpmath
mpmath.mp.dps = 110
e110["rays"][key]["xi_level_norm2_over_g2"]["10"] = mpmath.nstr(mpmath.mpf(x) * (1 + mpmath.mpf(10) ** -40), 110)
json.dump(e110, open(d / "out_engine_110.json", "w"), indent=1)
rc, so, se = run(d, ["identify.py"])
r2 = json.load(open(d / "results.json"))
out["M2_identify_with_corrupted_value"] = {"exit": rc, "failures": r2["checks"]["identification_failures"],
                                         "item4_value_written": r2["item4"][key]["norm2_over_g2_by_level"].get("10")}
print("M2:", out["M2_identify_with_corrupted_value"])

# ---------------- M3: quadrature DN broken, pipeline still passes?
d = fresh("M3")
patch(d / "quadcheck.py", "DN = s[..., None] * xi + (dot + dot.conj())[..., None] * Phi",
      "DN = s[..., None] * xi + (dot)[..., None] * Phi")
rc1, so1, se1 = run(d, ["quadcheck.py"])
rc2, so2, se2 = run(d, ["identify.py"])
r3 = json.load(open(d / "results.json"))
out["M3_quadrature_DN_broken"] = {"exit_quadcheck": rc1, "exit_identify": rc2,
                                  "max_quadrature_vs_engine_abs_dev": r3["checks"]["max_quadrature_vs_engine_abs_dev"],
                                  "failures": r3["checks"]["identification_failures"]}
print("M3:", out["M3_quadrature_DN_broken"])

# ---------------- M4: item5 'control' pointed at a zero level
d = fresh("M4")
patch(d / "item5_zeros.py", 'for name, v, J in (("R2", R2, 5), ("R3", R3, 7), ("R4", R4, 9)):',
      'for name, v, J in (("R2", R2, 4), ("R3", R3, 7), ("R4", R4, 9)):')
rc, so, se = run(d, ["item5_zeros.py"])
ctl = [l for l in so.splitlines() if "control" in l]
out["M4_item5_control_at_zero_level"] = {"exit": rc, "control_lines": ctl}
print("M4:", out["M4_item5_control_at_zero_level"])

json.dump(out, open(R / "stage2_mutations.json", "w"), indent=1, default=str)
print("written stage2_mutations.json")
