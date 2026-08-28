"""M5.32 R8 AUDIT MERGE: assembles ../data/m5_32_r8_audit.json from the
checkpoints written by m5_32_r8_audit.py. One entry per claim C1..C10.
"""
from __future__ import annotations

import json
import os
import time

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
CKPT = os.path.join(HERE, "..", "checkpoints", "m5_32_r8")


def ck(tag):
    with open(os.path.join(CKPT, f"{tag}.json")) as f:
        return json.load(f)


def rel(a, b):
    d = max(abs(a), abs(b), 1e-300)
    return abs(a - b) / d


def main():
    c1 = ck("c1_degree")
    lad = ck("ladder")
    ladm = ck("ladder_mutation")
    lada = ck("ladder_avgjets")
    sh = ck("shell")
    st = ck("string")
    hr = ck("h_refine")
    ov = ck("overlap")
    th = ck("theorem")
    tl = ck("tail")
    cf = ck("coeff")

    T = lad["terms"]
    C5T = ("Q_I1sq", "Q_I4sq", "Q_Fpair")
    out = {}

    # ---------------- C1 ----------------
    bi = c1["terms"]["Q_BI"]
    out["C1"] = {
        "verdict": "QUALIFIED",
        "own_number": bi["deg6_rel_to_C4"],
        "producer_number": 1e-12,
        "rel_dev": rel(bi["C4"], -5.7115152607423176e-05),
        "note": ("6 terms exactly quartic, odd part 0 at 1e-15 over 9 omegas and "
                 "a degree-6 fit; Born-Infeld is not a polynomial, its degree-6 "
                 "content is 0.16% of C4 and its extracted C4 moves 2.1%.")}

    # ---------------- C2 ----------------
    fp_slope = sh["shells"]["c2_Q_Fpair"]["slope"]
    out["C2"] = {
        "verdict": "QUALIFIED",
        "own_number": fp_slope,
        "producer_number": 0.09601494989272187,
        "rel_dev": rel(T["Q_Fpair"]["C2_exponent_in_L"], 0.09601494989272187),
        "note": ("A's L-exponents reproduce, but A itself diverges as h^-5.00 "
                 "so its L-convergence is empty; the omega^2 integrals are not "
                 "IR-finite either, 80 to 86% sits on the z-axis string and "
                 "Q_Fpair's shell integrand decays only as r^-0.64. Negligible "
                 "vs kin regardless (0.33%).")}

    # ---------------- C3 ----------------
    out["C3"] = {
        "verdict": "CONFIRMED",
        "own_number": lad["I1_control"]["kin_L48"],
        "producer_number": 426.5070121483972,
        "rel_dev": lad["I1_control"]["rel_dev"],
        "note": ("independent I1 reproduces the certified inertia to 2.7e-16 and "
                 "C4(I1) = 9.5e-14; the eta -> identity mutation breaks it "
                 "(rel dev 2.0), so the control gate is falsifiable.")}

    # ---------------- C4 ----------------
    out["C4"] = {
        "verdict": "QUALIFIED",
        "own_number": hr["exponents_in_h"]["Q_I1sq_C4"],
        "producer_number": 0.6093753326161384,
        "rel_dev": rel(T["Q_I1sq"]["C4_exponent_in_L"], 0.6093753326161384),
        "note": ("every C4 and exponent reproduces to 1e-4, but 78 to 98% of C4 "
                 "sits in the 256 cells with rho < h and it diverges as h^-1.65 "
                 "at fixed L, so it is a UV string artifact, not an IR effect; "
                 "the exponent also halves under the other branch read.")}

    # ---------------- C5 ----------------
    out["C5"] = {
        "verdict": "CONFIRMED",
        "own_number": T["Q_C6a"]["C4_exponent_in_L"],
        "producer_number": 3.0000000000000466,
        "rel_dev": rel(T["Q_C6a"]["C4_exponent_in_L"], 3.0000000000000466),
        "note": ("exponent 3.0 to 1e-13, C4 3583.18 -> 28665.4 (ratio 8.000); "
                 "strengthened: h-independent to 4e-14, shell slope +1.990, tube "
                 "fraction equals the cell fraction, so genuinely uniform. C6a "
                 "and C6b share one identical C4, not two.")}

    # ---------------- C6 ----------------
    sw = sh["tube_sweep"]["exclude_rho_lt_2.01h"]
    out["C6"] = {
        "verdict": "REFUTED",
        "own_number": sw["c2_squared"]["slope"],
        "producer_number": -0.292,
        "rel_dev": rel(sw["c2_squared"]["slope"], -0.292),
        "note": ("sign-definiteness holds (c2 <= 0 in 100% of cells). The "
                 "mechanism does not: excluding rho < 2h restores slope -1.88, "
                 "the r^-2 the claim says a naive r^-4 density would give. The "
                 "slow decay is the z-axis string (73.6% of c2^2).")}

    # ---------------- C7 ----------------
    g = th["generators"]
    out["C7"] = {
        "verdict": "QUALIFIED",
        "own_number": g["rot_12"]["tangent_max_abs"],
        "producer_number": 1.3,
        "rel_dev": rel(g["rot_12"]["tangent_max_abs"], 1.3),
        "note": ("conclusion and algebra survive (min 0.3), but the quoted "
                 "magnitudes come from X M - M X^T, antisymmetric for diagonal M "
                 "and not a tangent of symmetric M; true values 0.7/1.0/0.3 and "
                 "33.0/32.3/32.0. The escape kills the clock: kin = 0 exactly.")}

    # ---------------- C8 ----------------
    out["C8"] = {
        "verdict": "CONFIRMED",
        "own_number": tl["jet_exponent"],
        "producer_number": -0.9692154755203344,
        "rel_dev": rel(tl["jet_exponent"], -0.9692154755203344),
        "note": ("deviation flat at 0.960902353693305 (exponent 3e-18), a0 flat "
                 "at 0.4242641, M eta spectrum equals the vacuum spectrum to "
                 "1.3e-15 at every shell; a forced-decay mutation gives exponent "
                 "-1.81, so the detector can see decay.")}

    # ---------------- C9 ----------------
    hex_c5 = None
    hs = [r["h"] for r in hr["rows"]]
    reqs = [4.0 * abs(r["I1_C2"]) / abs(r["Q_I1sq_C2"]) for r in hr["rows"]]
    import numpy as np
    hex_c5 = float(np.polyfit(np.log(hs), np.log(reqs), 1)[0])
    out["C9"] = {
        "verdict": "REFUTED",
        "own_number": hex_c5,
        "producer_number": 1.0280978827865082,
        "rel_dev": rel(cf["per_term"]["Q_I1sq"]["exponent_in_L"],
                       1.0280978827865082),
        "note": ("the ladder reproduces exactly at delta 0.3 h 1.5, but C2 of "
                 "every C5 term diverges as h^-3.0, so the required c5 falls as "
                 "h^+2.99 (312.6 -> 39.8 -> 11.7) and vanishes in the continuum; "
                 "at delta 1 its L-exponent drops to 0.135.")}

    # ---------------- C10 ----------------
    fj = cf["fixed_J_omega_star"]["rows"]["Q_I1sq"]
    r0 = [x for x in fj if x["c5"] == 0.0][0]
    r300 = [x for x in fj if x["c5"] == 300.0][0]
    out["C10"] = {
        "verdict": "CONFIRMED",
        "own_number": r300["omega_star_L_spread"],
        "producer_number": 13.85095520344533,
        "rel_dev": rel(
            cf["per_term"]["Q_I1sq"]["omega_where_quartic_matches_at_c5_1"][0],
            13.85095520344533),
        "note": ("solved the fixed-J cubic instead of comparing coefficients: "
                 "omega* L spread is 1.04 at c5 = 0 and still 1.31 at c5 = 300, "
                 "so omega* stays ~1/L. Negative c5 and Born-Infeld both give a "
                 "NEGATIVE quartic and an unbounded H.")}

    for k in out:
        out[k]["own_number"] = float(out[k]["own_number"])
        out[k]["producer_number"] = float(out[k]["producer_number"])
        out[k]["rel_dev"] = float(out[k]["rel_dev"])

    vs = [v["verdict"] for v in out.values()]
    doc = {
        "task": "M5.32 R8 INDEPENDENT ADVERSARIAL AUDIT (classes C5 and C6)",
        "auditor_note": ("the producer's three R8 scripts were not read and not "
                         "imported; all six quartic densities plus I1 and "
                         "Born-Infeld were reimplemented from the written "
                         "definitions and reproduce the producer's A, C2 and C4 "
                         "to 1e-4 on every polynomial term, which validates the "
                         "instrument before the attack"),
        "claims": out,
        "n_confirmed": vs.count("CONFIRMED"),
        "n_qualified": vs.count("QUALIFIED"),
        "n_refuted": vs.count("REFUTED"),
        "summary": (
            "The arithmetic is right and reproduces bit-close on every "
            "polynomial term, and the two clean results (C3 the certified "
            "control, C5 the volume-extensive C6 class, C8 the orbit tail, C10 "
            "the 1/L clock) survive. What does not survive is the physical "
            "reading of the C5 class. The hedgehog ansatz M = Qh d4 Qh^T is "
            "discontinuous along the whole z axis (a phi-ring of M has spread "
            "0.15 at rho = 1e-4 at every z), and 73 to 98 percent of every C5 "
            "omega^4 coefficient, and 80 to 86 percent of every C5 omega^2 "
            "coefficient, sits in the 256 lattice columns adjacent to that "
            "line, which is 0.098 percent of the cells. Remove a two-cell tube "
            "and the shell integrand returns to r^-1.88, the convergent r^-2 "
            "the claim says does not happen. At fixed L every C5 coefficient "
            "diverges under refinement (C4 as h^-1.65 to h^-2.03, C2 as "
            "h^-3.0) while the certified I1 coefficients are h-stable at "
            "h^-0.05, so the C5 quartics have no continuum limit on this "
            "ansatz and their box growth is not an IR statement. The one "
            "coefficient the campaign would act on, the c5 that flips the "
            "energy's omega^2 sign, is therefore not a physical number: it "
            "scales as h^+2.99. The C6 class is the opposite case and is "
            "clean: its omega^4 inertia is exactly h-independent, exactly "
            "uniform in space and exactly volume extensive, and the two C6 "
            "terms share one and the same coefficient. C7's conclusion holds "
            "but its numbers come from a map that is not the tangent, and the "
            "degenerate-vacuum escape is closed by a theorem the producer did "
            "not state: on M = Q d4 Q^T with Q in SO(1,3), a generator that "
            "annihilates the vacuum annihilates the whole hedgehog, so the "
            "clock inertia is exactly zero (measured at delta = 0)."),
        "new_findings": [
            "the hedgehog ansatz is discontinuous along the entire z axis: the "
            "spread of M over a phi ring is 0.15 at rho = 1e-4 for z = 5 to 40, "
            "so the Euler-angle frame carries a coordinate string to the wall",
            "73.6 to 97.7 percent of every C5 omega^4 coefficient sits in the "
            "256 cells with rho < h (0.098 percent of the box); excluding a "
            "two-cell tube restores the shell slope to -1.88, i.e. convergent",
            "the C5 STATIC coefficient diverges as h^-5.001 for all three "
            "terms (1.71 -> 54.7 -> 415.4 at L = 48), so C2's IR-finiteness of "
            "A is a statement about a quantity with no continuum limit",
            "no continuum limit: at fixed L = 48 the C5 omega^4 coefficients "
            "diverge as h^-1.65 (Q_I1sq), h^-2.03 (Q_I4sq), h^-1.79 (Q_Fpair) "
            "while the certified kin is h-stable at h^-0.05",
            "the C5 omega^2 coefficients diverge faster still, as h^-3.0, so the "
            "c5 required to flip the energy's omega^2 sign scales as h^+2.99 "
            "(312.6 at h = 1.5, 39.8 at h = 0.75, 11.7 at h = 0.5)",
            "Q_C6a and Q_C6b have the same omega^4 coefficient at every box, agreeing "
            "to 4e-16 (3583.180704, 12093.234878, 28665.445636): both reduce to "
            "h^3 sum tr(a0 eta a0 eta)^2, so C5 is one result, not two",
            "the quartic C4 is branch-read dependent: squaring the averaged "
            "density instead of averaging the squares moves Q_I1sq's exponent "
            "from 0.609 to 0.189 and Q_I4sq's C4 by a factor 15",
            "the producer's generator magnitudes come from X M - M X^T, which "
            "is antisymmetric for diagonal M and cannot be a tangent of the "
            "symmetric field; the true tangents are 0.7/1.0/0.3 and 33.0/32.3/32.0",
            "escape theorem: for M = Q d4 Q^T with Q in SO(1,3) the co-moving "
            "flow is Q (X d4 + d4 X^T) Q^T, so a vacuum-annihilating generator "
            "annihilates the hedgehog; measured at delta = 0, a0 = 0 exactly and "
            "kin = 0, while the unconjugated generator stays flat at 0.708",
            "at delta = 1 the required c5 is nearly box-stable (exponent 0.135, "
            "22.0 to 24.1 over a factor 2 in L) and the quartic match omega "
            "falls to 0.813, i.e. 1.03 times the radiation window, so C9's and "
            "C10's margins are specific to delta = 0.3",
            "negative c5 opens nothing and costs boundedness: C2_q < 0 for every "
            "C5 term so H2 only grows, while H4 = 3 c5 C4_q turns negative and "
            "H(omega = 30) = -2.3e6; Born-Infeld is in that class, C4 = -5.6e-5",
            "the interior-cell overlap test passes exactly (max abs difference "
            "0.0 between the n = 32 and n = 64 boxes), so the box growth is "
            "added volume and NOT a boundary or stencil artifact",
            "the certified control is reproduced to 2.7e-16 and the C6 class "
            "shell integrand is +1.990, both with mutations that fire, so the "
            "instrument is calibrated in both directions"],
        "runtime_s": 255.0}
    with open(os.path.join(DATA, "m5_32_r8_audit.json"), "w") as f:
        json.dump(doc, f, indent=1)
    print(json.dumps({k: {"verdict": v["verdict"], "own": v["own_number"],
                          "prod": v["producer_number"]}
                      for k, v in out.items()}, indent=1))
    return doc


if __name__ == "__main__":
    main()
