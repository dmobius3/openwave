"""Builds the S1 disclosure inventory by sweeping the local records, and fails if anything relevant is unclassified.

Sweep: every file in M8_DYNAMICS/out and M8_S whose text mentions a Hessian, a saddle, a census, a critical orbit,
point or set, an ascent or descent run, or one of the census values. Each hit must appear in CLASSIFY below, with
what it exposes; an unclassified hit fails the run. Binary records (npz) are included by name.
Usage: python3 make_disclosure_inventory.py [dynamics_dir] [out_file]
"""
import hashlib
import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
DYN = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, '..', 'M8_DYNAMICS', 'out')
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, 'DISCLOSURE_INVENTORY.md')
TERMS = re.compile(r'hess|saddle|census|critical (orbit|point|set)|distinct orbits|minimis|minimiz|ascent|descent|215\.77|238\.01|238\.35|924', re.I)

CLASSIFY = {
    'soft_search.txt': 'Ascent runs reach the maximum of Q at the hexagon, in both sectors (2026-08-31).',
    'isotropy_hessian.txt': 'The numerical Hessian of Q at the hexagon, both sectors: 9 negative, 5 null, 0 positive of 14. Its verdict line miscounts. Also a spin-7/2 block, which is not S1 object.',
    'isotropy_final.txt': 'The same run with the full negative spectrum at the hexagon, nine eigenvalues with multiplicities.',
    'minimizer_refine.txt': 'Descent reaches the minimum of Q at the coherent orbit, both sectors (2026-09-01).',
    'minimizer_continuous.txt': 'Companion descent round.',
    'minimizer_orbit_rank.txt': 'Companion descent round.',
    'critical_minimizer.txt': 'Companion descent round.',
    'critical_kills16.txt': 'Level-16 question at the minimizer, not a second variation.',
    'stabilizer_minimizer.txt': 'Stabilizer of the minimizer, level-16 question.',
    'stabilizer_maximal.txt': 'Stabilizer of the maximizer, level-16 question.',
    'pencil_census.txt': 'The random census, first round: critical orbits with block-operator spectra (2026-09-04).',
    'pencil_census2.txt': 'The random census, second round: twelve critical orbits in both sectors, v1 missed (2026-09-05).',
    'pencil_census_r4.npz': 'The saved states of the random census, sector 3prime.',
    'pencil_census_r5.npz': 'The saved states of the random census, sector 4.',
    'pencil_lemma.txt': 'Part C is the seeded census: thirteen orbits with block-operator spectra, v1 included, and its own note that thirteen is a lower bound. Six of its spectra are exact integers, which is what M8.11 disclosed.',
    'pencil_criticals.txt': 'Criticality of A6 eigenvectors, where the label "saddle" for the zonal comes from; no Hessian.',
    'pencil_saddle.txt': 'The saddle test on A6 eigenvectors.',
    'pencil_p1.txt': 'Rounds naming u_max, u_min and saddle.',
    'pencil_p1_v1.txt': 'Rounds naming u_max, u_min and saddle.',
    'pencil_p6.txt': 'Rounds naming u_max, u_min and saddle.',
    'pencil_round3.txt': 'Refines u_min and the saddle, and their stabilizers.',
    'pencil_rational.txt': 'Rational reconstruction attempts at those states.',
    'pencil_a5.txt': 'Companion round.',
    'pencil_audit.txt': 'Companion round.',
    'pencil_exact.txt': 'Exact arithmetic for the cubic covariant T, not a second variation.',
    'pencil_adjudicate.txt': 'Perturbation tests of sigma_6/sigma_1 at the hexagon, a different functional.',
    'pencil_cold.txt': 'Numerical check that the four tabulated rays are critical (2026-09-07).',
    'pencil_c3.txt': 'C3-plane search, candidate values from a finite-difference gradient at its noise floor; unverified.',
    'pencil_c3id.txt': 'Configurations of the states found in the C3 plane.',
    'r2_flats.txt': 'Transverse spectrum of a spin-7/2 block, which is not S1 object.',
    'transform_equivalence.txt': 'Transform equivalence checks, no critical data.',
    'departure_r4_channels.txt': 'Channel decomposition, no second variation.',
    'pencil_domain.txt': 'Domain of the pencil family, no critical data.',
    'pencil_field.txt': 'Field of the pencil values, no critical data.',
    'translate_rho6_log.txt': 'S0 aid: the translation identities and the map to the condensate couplings (2026-09-19).',
    'check_extrema_log.txt': 'S0 aid: the invariant form, and multi-start optimizations giving the global maximum and minimum with their maximizers.',
    'check_census_inputs_log.txt': 'S0 aid: the closed forms at the M8 couplings, criticality checks at the H, E, I and B points, and the D2/v1 tie test.',
    'check_s3_maximum_log.txt': 'S0 lemma: the global maximum at the hexagon, proved.',
    'identify_census_states_log.txt': 'S0 aid: criticality and isotropy of the census states, frame-independent.',
    's1_L0_log.txt': 'S1 step 1: the fixed-locus classification, and the ranges of r6 on each line.',
    's1_L1_log.txt': 'S1 step 2: the exact critical set on the seven lines, the in-line types, and the ten orbits.',
    'departure_n18.txt': 'Newton continuation of the branch at mu = 1.202360, the hexagon in sector 4: the a ~ 8.7 continuation (2026-08-31).',
    'departure_n18_ext.txt': 'The same continuation, extended range.',
    'departure_n18_states.txt': 'The same continuation, with states.',
    'departure_n18_border.txt': 'The same continuation, border rows (2026-09-01).',
    'orbit_guard.txt': 'Orbit-rank guard for the continuation border rows; no critical-set data.',
    'stabilizer_minimizer_PARTIAL.txt': 'Partial stabilizer run at the minimizer, level-16 question.',
    'pencil_c3b.txt': 'The C3-plane search, second round: the same unverified finite-difference candidates as pencil_c3.',
    'pencil_fix3.txt': 'The D3 (eclipsed) family restriction of r6 in another parameterization, and the octahedron on it (2026-09-07). The paper Section 5.8 publishes that restriction.',
    'pencil_L17e.txt': 'The exact multipole table ||rho_K||^2 at named states, including the coherent state.',
    'pencil_sec5check.txt': 'Checks of Section 5 identities, including 924||rho_6(v_m)||^2 = 1, 36, 225, 400, which the paper publishes.',
    'pencil_L17f.txt': 'Gradient identity checks for the cubic covariant; paper side, no critical-set data.',
    'pencil_P.txt': 'Exact-arithmetic checks of a proof about the covariant; paper side.',
    'pencil_family.txt': 'The four-parameter family basis; paper side.',
    'pencil_fix.txt': 'A sign-defect repair in the covariant table; paper side.',
    'pencil_gam2.txt': 'Sparsity of Gamma; paper side.',
    'pencil_gamma.txt': 'Euler identities and the tangential Gamma; paper side.',
    'pencil_solve.txt': 'Monomial solve for the covariant; paper side.',
    'pencil_tensor.txt': 'Reflection signs; paper side.',
    'pencil_w6.txt': 'Builds 2I and its spin-3 representation, and the weight w6; paper side.',
    'pencil_w6b.txt': 'Parseval for the multipole transform; paper side.',
    'P1_package.json': 'A packaged covariant table with its supersession note; paper side, no critical-set data.',
}
SKIP = re.compile(r'^(DISCLOSURE_INVENTORY|S0_|S1_|make_disclosure)')
rows, unclassified = [], []
for d in (DYN, HERE):
    for fn in sorted(os.listdir(d)):
        path = os.path.join(d, fn)
        if not os.path.isfile(path) or SKIP.match(fn) or fn.endswith('.py'):
            continue
        if fn.endswith('.npz'):
            hit = 'census' in fn
        else:
            try:
                hit = bool(TERMS.search(open(path, encoding='utf-8', errors='ignore').read()))
            except OSError:
                hit = False
        if not hit:
            continue
        h = hashlib.sha256(open(path, 'rb').read()).hexdigest()
        when = time.strftime('%Y-%m-%d', time.localtime(os.path.getmtime(path)))
        if fn not in CLASSIFY:
            unclassified.append(fn)
        rows.append((fn, when, CLASSIFY.get(fn, 'UNCLASSIFIED'), h))

control = [r for r in rows if r[0] == 'pencil_census2.txt']
ok_control = bool(control)
print(('PASS ' if ok_control else 'FAIL ') + 'control: the sweep finds pencil_census2.txt, a known-present record')
print(('PASS ' if not unclassified else 'FAIL ') + f'every swept file is classified (unclassified: {unclassified})')
body = ['# S1 disclosure inventory', '',
        f'Generated by `make_disclosure_inventory.py` on {time.strftime("%Y-%m-%d")}. The sweep covers `M8_DYNAMICS/out/` and `M8_S/`,',
        'and it fails if a swept file is not classified here. Scripts are pinned by hash in their step notes; this table is the records.', '',
        '| file | date | what it exposes | SHA-256 |', '|---|---|---|---|']
for fn, when, what, h in rows:
    body.append(f'| `{fn}` | {when} | {what} | `{h[:16]}…` |')
body += ['', f'{len(rows)} records. Published beside them: M8.11 in-locus second variations and the feasibility check signs.',
         '', 'M8.11 disclosure named exact spectra at six critical orbits, which are the six integral spectra of the seeded census.',
         'It did not name either census. The S1 pre-registration says so in a dated statement, and M8.11 frozen text is not edited.']
open(OUT, 'w', encoding='utf-8').write('\n'.join(body) + '\n')
print(f'{len(rows)} records written to {os.path.basename(OUT)}')
print('ALL PASS' if ok_control and not unclassified else 'FAILED')
