"""Stage 2, task 3: mutation tests of the other worker's PASS checks, run in copies.

Each mutation: copy solver_rerun/ (code + regenerated results/) to mut/<name>/, apply one textual
edit (the old text must occur exactly once) to one file, run the listed scripts with the same
interpreter and flags as ./py (sys.executable -S, as the solver's run_all.py does), and record
the PASS/FAIL lines that mention the target.  Nothing in solver_work/ or solver_rerun/ is edited.
Usage: ./py st2_mutations.py [name ...]"""
import json, os, re, shutil, subprocess, sys, time

SRC = "solver_rerun"
MUT = "mut"


def edit_file(path, old, new):
    t = open(path).read()
    assert t.count(old) == 1, (path, old, t.count(old))
    open(path, "w").write(t.replace(old, new))


def edit_json(path, fn):
    d = json.load(open(path))
    fn(d)
    json.dump(d, open(path, "w"), indent=1)


def fresh(name):
    dst = os.path.join(MUT, name)
    if os.path.exists(dst):
        shutil.rmtree(dst)
    shutil.copytree(SRC, dst, ignore=shutil.ignore_patterns("__pycache__"))
    return dst


def run(dst, script):
    p = subprocess.run([sys.executable, "-S", script], cwd=dst, capture_output=True, text=True)
    lines = [l for l in p.stdout.splitlines() if l.startswith("PASS") or l.startswith("FAIL")]
    return p.returncode, lines, p.stdout, p.stderr


MUTATIONS = {}


def mutation(f):
    MUTATIONS[f.__name__] = f
    return f


# ---- M1: s09 attaches 'ANTIUNITARY' without reading any antiunitary verification
@mutation
def m1_s09_antiunitary_unverified(dst):
    def f(d):
        for k in ("U5_sector3", "U5_sector4", "U6_sector3", "U6_sector4"):
            d[k]["antiunitary_residual"] = "1.0"      # as if A' did NOT fix Phi, xi, E1
    edit_json(os.path.join(dst, "results", "numeric.json"), f)
    return ["s09_zeros.py"], r"ANTIUNITARY|no unresolved"


# ---- M2: s09 accepts spurious zeros at places where the named reason cannot apply
@mutation
def m2_s09_spurious_zero(dst):
    def f(d):
        d["U1_sector3"]["form_cross"] = "0"             # U1 has no tau_x, tau_y
        d["U1_sector3"]["DNxi_perp_comp_ie_t"] = "0"    # U1 has no e_t
        d["U6_sector3"]["form_ie_t"] = "0"              # U6 has no e_t
    edit_json(os.path.join(dst, "results", "numeric.json"), f)
    return ["s09_zeros.py"], r"^(PASS|FAIL) U1 sector 3: (form cross|DNxi_perp_comp_ie_t)|^(PASS|FAIL) U6 sector 3: form on i e_t|no unresolved"


# ---- M3: s10 'H0' check with a level below 6 inserted into the 4-dim sector
@mutation
def m3_s10_level_below_6(dst):
    def f(d):
        d["dims"]["4"]["2"] = 1
        d["dims"]["4"]["4"] = 1
    edit_json(os.path.join(dst, "results", "item1.json"), f)
    return ["s10_existence.py"], r"48 is an eigenvalue|H1|H2"


# ---- M4: s03 'every element of the generated group fixes [u]' with a wrong character
@mutation
def m4_s03_wrong_character(dst):
    edit_file(os.path.join(dst, "common.py"),
              "'U4': [(rot((0, 0, 1), pi / 3), mpc(-1)), (rot((1, 0, 0), pi), mpc(-1))],",
              "'U4': [(rot((0, 0, 1), pi / 3), mpc(1)), (rot((1, 0, 0), pi), mpc(-1))],")
    return ["s03_stabilizers.py"], r"U4"


# ---- M5: s04 antiunitary check with the minus sign dropped from A' (should FAIL: positive control)
@mutation
def m5_s04_antiunitary_sign(dst):
    edit_file(os.path.join(dst, "s04_pipeline.py"),
              "return {J: [vscale(mat_vec(rho_left(J, g), theta(u)), -1) for u in us] for J, us in sec_vecs.items()}",
              "return {J: [vscale(mat_vec(rho_left(J, g), theta(u)), 1) for u in us] for J, us in sec_vecs.items()}")
    return ["s04_pipeline.py", "s09_zeros.py"], r"antiunitary|ANTIUNITARY|all exact values|no unresolved"


# ---- M6: s04 kappa solve corrupted (pseudo-inverse eigenvalue doubled): is the with-kappa check independent?
@mutation
def m6_s04_kappa_solve(dst):
    edit_file(os.path.join(dst, "s04_pipeline.py"),
              "kap = [mp.fsum(V[i, a] * mp.fsum(V[k, a] * b[k] for k in range(n)) / ev[a]",
              "kap = [mp.fsum(V[i, a] * mp.fsum(V[k, a] * b[k] for k in range(n)) / (2 * ev[a])")
    return ["s04_pipeline.py"], r"(U5|U6) s3 (order-a\^5|with kappa)"


# ---- M7: item-12 Casimir factor 4 -> 3.9: does the 'residual equation' line see it?
@mutation
def m7_s12_casimir(dst):
    edit_file(os.path.join(dst, "s12_quadrature.py"), "    return vscale(out, 4)\n", "    return vscale(out, mpf('3.9'))\n")
    # restrict the run to U1 in both sectors to keep it short (same code path)
    edit_file(os.path.join(dst, "s12_quadrature.py"),
              "        for key in ('U1', 'U2', 'U3', 'U4', 'U5', 'U6'):\n            Phi = sec.block(U[key])",
              "        for key in ('U1',):\n            Phi = sec.block(U[key])")
    return ["s12_quadrature.py"], r"U1 sector"


# ---- M8: control for M7 (unmutated code, same restriction): all U1 lines PASS
@mutation
def m8_s12_control(dst):
    edit_file(os.path.join(dst, "s12_quadrature.py"),
              "        for key in ('U1', 'U2', 'U3', 'U4', 'U5', 'U6'):\n            Phi = sec.block(U[key])",
              "        for key in ('U1',):\n            Phi = sec.block(U[key])")
    return ["s12_quadrature.py"], r"U1 sector"


# ---- M9: s01 'exact dims vanish at every odd n' can only fail through chi(-h) != chi(h)
@mutation
def m9_s01_odd_levels(dst):
    # replace the exact dims by garbage at even n only: the odd-n check must still pass
    edit_file(os.path.join(dst, "s01_levels.py"),
              "    dims[name] = row\n",
              "    row = {k: (v + 7 if k % 2 == 0 else v) for k, v in row.items()}\n    dims[name] = row\n")
    return ["s01_levels.py"], r"odd n|numerical intertwiner counts"


# ---- M10: results.json field 'xi_transforms_by_same_character' after the covariance check fails
@mutation
def m10_s99_xi_character_literal(dst):
    edit_file(os.path.join(dst, "common.py"),
              "'U4': [(rot((0, 0, 1), pi / 3), mpc(-1)), (rot((1, 0, 0), pi), mpc(-1))],",
              "'U4': [(rot((0, 0, 1), pi / 3), mpc(1)), (rot((1, 0, 0), pi), mpc(-1))],")
    return ["s04_pipeline.py", "s99_assemble.py"], r"U4 s(3|4) xi"


# ---- M11: results.json field 'sign_of_lambda4' after lambda4 is made positive in numeric.json
@mutation
def m11_s99_sign_literal(dst):
    def f(d):
        for k, r in d.items():
            if "lambda4_over_g2_with_kappa" in r:
                r["lambda4_over_g2_with_kappa"] = r["lambda4_over_g2_with_kappa"].lstrip("-")
                r["lambda4_over_g2_without_kappa"] = r["lambda4_over_g2_without_kappa"].lstrip("-")
    edit_json(os.path.join(dst, "results", "numeric.json"), f)
    return ["s99_assemble.py"], r"."


if __name__ == "__main__":
    names = sys.argv[1:] or list(MUTATIONS)
    os.makedirs(MUT, exist_ok=True)
    os.makedirs("out2", exist_ok=True)
    logp = "out2/mutations.json"
    allres = json.load(open(logp)) if os.path.exists(logp) else {}
    for name in names:
        t0 = time.time()
        dst = fresh(name)
        scripts, pat = MUTATIONS[name](dst)
        rec = {"scripts": scripts, "lines": [], "returncodes": []}
        for sc in scripts:
            rc, lines, so, se = run(dst, sc)
            rec["returncodes"].append(rc)
            rec["lines"] += [l for l in lines if re.search(pat, l[5:] if pat.startswith("^") is False else l)]
            if rc != 0:
                rec["stderr_tail"] = se[-1500:]
        rec["n_PASS"] = sum(l.startswith("PASS") for l in rec["lines"])
        rec["n_FAIL"] = sum(l.startswith("FAIL") for l in rec["lines"])
        rec["seconds"] = round(time.time() - t0)
        allres[name] = rec
        print("==", name, "PASS", rec["n_PASS"], "FAIL", rec["n_FAIL"], "rc", rec["returncodes"], "%ds" % rec["seconds"])
        for l in rec["lines"][:40]:
            print("   ", l[:170])
        json.dump(allres, open(logp, "w"), indent=1)
