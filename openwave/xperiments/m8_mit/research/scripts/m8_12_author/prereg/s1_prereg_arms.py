"""The mutation suite for the S1 pre-registration gate: every arm plants one defect and requires the gate it was
built to trip to fail, not merely some gate.

Each arm runs against a verified-green parent, in a scratch copy of the program folder, so nothing here can touch
the real files. An arm that plants a defect the gate does not catch is itself reported as a failure.
Usage: python3 s1_prereg_arms.py
"""
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
M8S = os.path.dirname(HERE)
SUPPORT = ['s1_L0_log.txt', 's1_L1_log.txt', 's1_H_log.txt', 's1_L0.py', 's1_L1.py', 's1_H.py', 's1_control.py',
           's1_control_log.txt', 'DISCLOSURE_INVENTORY.md', 'S1_STEP1_L0.md', 'S1_STEP2_L1.md', 'S1_STEP3_H.md',
           'S0_S3_MAXIMUM.md']
# (file, find, replace, what the arm plants)
ARMS = [
    ('task', '`v₃ + v₋₃` | 463 |', '`v₃ + v₋₃` | 464 |', 'a census value'),
    ('task', '463 | 9 | (9, 0, 0)', '463 | 9 | (9, 1, 0)', 'a census signature'),
    ('task', '288 | 9 | (6, 0, 3) | 6 | 3 |', '288 | 9 | (6, 0, 3) | 6 | 4 |', 'an index for g < 0'),
    ('task', '`v₂ + 2·v₋₁`', '`v₂ + 3·v₋₁`', 'an orbit representative'),
    ('task', '| zonal `v₀` | `v₀` | 400 | 10 | (8, 0, 2) | 8 | 2 |', '', 'a census row deleted'),
    ('task', '671187', '671189', 'a characteristic-polynomial coefficient'),
    ('task', '−(125/132)s² + (10/11)s + 3/77', '−(125/132)s² + (10/11)s + 4/77', 'a quadratic, in the claim row'),
    ('task', '| C₅ `{v₃, v₋₂}` | `−(125/132)s²', '| C₅ `{v₃, v₋₂}` | `−(125/133)s²', 'a quadratic, in the frozen table'),
    ('task', 's* = 12/25', 's* = 13/25', 'an in-line vertex'),
    ('task', '148y² + 463', '148y² + 461', 'a dihedral chart formula'),
    ('task', '`z = ±√230/10`', '`z = ±√231/10`', 'a chart critical point'),
    ('task', '`z = ±i√30/5`; `z = ∞`', '`z = ±i√30/5`', 'the point at infinity dropped'),
    ('task', '−(8/33)s² + 75/308', '−(8/33)s² + (1/100)s + 75/308', 'the exactly-zero coefficient made nonzero'),
    ('task', '57 records', '58 records', 'the inventory size'),
    ('task', 'Negative controls', 'Removed controls', 'a claim group dropped'),
    ('task', '## THE FIREWALL', '## THE WALL', 'the firewall section renamed'),
    ('task', '## DISCLOSURE', '## NOTES', 'the disclosure section renamed'),
    ('task', 'Never "stable"', 'The branch is unstable. Never "stable"', 'a bare short form'),
    ('task', '−104/55', '−104/56', 'an M8.11 parent value, first copy'),
    ("task", "M8.11's `−104/55`, and the null orbit", "M8.11's `−104/56`, and the null orbit", 'the same value, second copy'),
    ('task', '920/473', '920/474', 'the prism parent value'),
    ('task', '−56/165', '−56/166', 'the sector-bridge value'),
    ('task', '`s1_L0.py` | 128 checks', '`s1_L0.py` | 129 checks', 'a step check count'),
    ('task', 'It now runs 16 checks', 'It now runs 9 checks', "the dated correction's count"),
    ('task', '| 7 checks, on top of re-running step 3 |', '| 8 checks, on top of re-running step 3 |', 'the filed row rewritten, correction intact'),
    ('task', '**Correction, 2026-09-20, after filing.**', '**Note.**', 'the correction note undated'),
    ('task', 'TrN̄² ≤ 171/2', 'TrN̄² ≤ 171/3', 'a G2 bound'),
    ('task', '−5/231 − \\|f\\|²/22', '−5/231 − \\|f\\|²/23', 'the G2 invariant form'),
    ('task', '`C₁` (7), `C₂` (4 and 3), `C₃` (3)', '`C₁` (7), `C₂` (4 and 3)', 'a higher-dimensional space dropped'),
    ('task', 'Six points: `v₃`, `v₂`, `v₁`, `v₀`,', 'Six points: `v₃`, `v₂`, `v₁`,', 'a point class dropped'),
    ('task', '1, 36, 225 and 400 at the weight states', '1 and 400 at the weight states', 'a parent value dropped from P1'),
    ('task', 'the weight states from paper § 7.2', 'the weight states from paper § 5.8', 'P1 miscited'),
    ('task', '| The fourth bedrock paper, § 7.2 |', '| The fourth bedrock paper, § 7.9 |', 'the § 7.2 source row'),
    ('task', 'which also grades the arguments L0, L1, L2', 'which also grades the arguments L0, L2', 'L1 out of the audit set'),
    ('task', 'The audit grades L0, L1, L2, G1 and G2', 'The audit grades L0, L2, G1 and G2', 'L1 out of the definition of done'),
    ('task', 'may break this ordering', 'cannot fill the gaps', 'the retired monotonicity inference'),
    ('task', 'The run is therefore offline', 'The run may be online', 'the network posture'),
    ('task', 'items 5 to 7, are new', 'all items are new', 'the exposed items misstated'),
    ('task', 'after monic normalization', 'exactly', "H2's monic rule"),
    ('task', 'so `−224/165` and `−42/55` reproduce', 'so `−224/165` and `−21/55` reproduce', 'the bridge arithmetic'),
    ('task', '`924·r̂₆ = 204`', '`924·r̂₆ = 205`', 'the control value'),
    ('task', '`2√1002/297`', '`2√1003/297`', 'the control gradient'),
    ('task', '`2√(312√15 + 4170)/297`', '`2√(312√15 + 4171)/297`', 'the control residual, altered but still nonzero'),
    ('task', 'all six reported quantities reproduce, of which five are independent, and all four residuals are nonzero',
     'all three exact quantities reproduce, and the residual is nonzero', "N1's pass condition, scoring one residual by sign"),
    ('task', 'any of the six differs. The control itself ceases to discriminate only if', 'either vanishes, which means', "N1's fail condition weakened"),
    ('task', '`2\u221a(4170 \u2212 312\u221a15)/297` at `\u2212i\u00b7J_y u`', '`2\u221a(4171 \u2212 312\u221a15)/297` at `\u2212i\u00b7J_y u`', 'one frozen residual, still nonzero'),
    ('task', 'the four squares summing to `17368/29403`', 'the four squares summing to `17369/29403`', 'the sum rule'),
    ('task', 'M_u = Hess N(u) \u2212 4N(u)\u00b7I', 'H_u', 'the matrix renamed back to the second variation'),
    ('task', "carried onto it by the room's own stated change of chart", "carried onto it by the room's own stated change of chart, or otherwise", 'P2 chart device loosened'),
    ('task', 'a dihedral row either in the paper\'s chart', 'a dihedral row in any chart', 'L1 chart device loosened'),
    ('work', '`\u2016M_u d\u2016` for each of the four generators', 'whether the directions of `O_u` are annihilated by the matrix', 'item 11 back to a yes or no'),
    ('task', '| N1 | At the non-critical point the worklist supplies, `u = (v₃ + v₁ − v₋₂)/√3`',
     '| N1 | At `s = 1/2` on the C₃ line', 'N1 drifting from the point the worklist supplies'),
    ('task', '`(100\\|z\\|⁴ − 20x²', '`(100|z|⁴ − 20x²', 'a pipe left unescaped in a table'),
    ('work', '**11.** At `u = (v₃ + v₁ − v₋₂)/√3`', '**11.** At `u = (v₃ + v₁ − v₀)/√3`', 'the worklist point drifting from N1'),
    ('work', '**8.** Over the whole', '**8.** The extremes are 1/924 and 463/924. Over the whole', 'frozen values leaked into the handout'),
    ('work', '28/39', '28/40', 'the D7 input'),
    ('work', 'reported as missing, never as zero', 'reported as zero', 'the zeros rule'),
    ('work', '**4.** Unite the critical', '**4.** There are ten of them. Unite the critical', 'a count of the answer leaked'),
    ('work', 'u = (v₃ + v₁ − v₋₂)/√3', 'u = v₂ + v₋₂', 'a census vector used as the control point'),
    ('work', '**5b.**', '**5c.**', 'item 5b renamed'),
    ('work', 'every diagonal entry and every off-diagonal entry', 'every diagonal entry', 'the cross term dropped'),
    ('work', 'multiplied by each `w₆`', 'in each sector', 'the sector-scaled numbers dropped'),
    ('work', 'Condon-Shortley', 'some', 'the Clebsch-Gordan convention unpinned'),
    ('work', 'Wherever two of them have the same stabilizer, give that separation exactly', 'Show that two of them differ', 'the separation clause weakened'),
    ('work', 'no rotation carries one onto another', 'no rotation carries the two C₄ lines onto each other', 'the hard pair named in the handout'),
    ('work', 'first say whether its projection onto `N_u` is zero', 'project them onto `N_u`', '5b normalizing every projection again'),
    ('work', 'Where it is zero, say which direction of `O_u` it is', 'Where it is zero, skip it', 'the zero direction left unidentified'),
    ('work', 'as an orbit-null control', 'as a check', 'the orbit-null control unnamed'),
    ('work', 'Normalize the rest in `Re⟨·,·⟩`', 'Normalize them all in `Re⟨·,·⟩`', 'the zero projection normalized'),
    ('log', 'PASS   arm: at a critical point the gradient vanishes', 'FAIL   arm: at a critical point the gradient vanishes', "the control's own arm failing"),
    ('s0', 'TrN̄² ≤ 171/2', 'TrN̄² ≤ 171/2 ', 'the frozen S0 proof edited'),
    ('inventory', '<DELETE FILE>', '', 'the gate inventory deleted'),
    ('inventory', 'the exposure table states what was visible per orbit\n', '', 'a line cut from the pin, to match a deleted gate'),
]
EXPECT = {
    'the gate inventory deleted': 'gate inventory is present; re-pin deliberately',
    'a line cut from the pin, to match a deleted gate': 'gate inventory matches the pin',
    'a census value': 'hexagon        value, dim, signature and both indices match the log row exactl',
    'a census signature': 'hexagon        value, dim, signature and both indices match the log row exactl',
    'an index for g < 0': 'octahedron     value, dim, signature and both indices match the log row exactl',
    'an orbit representative': "C3 ray         representative 'v₂ + 3·v₋₁' is the vector the step-3 script use",
    'a census row deleted': 'the census table in the doc has ten rows (found 9)',
    'a characteristic-polynomial coefficient': "charpoly of prism         : coefficients ['473', '920', '473', '696', '1419',",
    'a quadratic, in the claim row': 'the doc restates 9 restrictions, and each is one the step-2 log printed',
    'a quadratic, in the frozen table': 'C5{v3,v-2}     restriction matches the step-2 log exactly',
    'an in-line vertex': 'C5{v3,v-2}     vertex matches, interior or endpoint as logged',
    'a dihedral chart formula': 'the doc restates 9 restrictions, and each is one the step-2 log printed',
    'a chart critical point': 'D3{v3+v-3,v0}    chart critical set matches the logged one',
    'the point at infinity dropped': 'the infinite chart point is named on both dihedral rows',
    'the exactly-zero coefficient made nonzero': 'C4{v3,v-1}     restriction matches the step-2 log exactly',
    'the inventory size': 'disclosure carries the inventory, the correction and the deferred orbits',
    'a claim group dropped': 'all claim groups are present',
    'the firewall section renamed': "every required section heading is present (missing: ['## THE FIREWALL'])",
    'the disclosure section renamed': "every required section heading is present (missing: ['## DISCLOSURE'])",
    'a bare short form': 'the short-form rule is stated, and no bare stable or unstable is used',
    'an M8.11 parent value, first copy': 'the doc carries the parent value −104/55 exactly 2 time(s)',
    'the same value, second copy': 'the doc carries the parent value −104/55 exactly 2 time(s)',
    'the prism parent value': 'the doc carries the parent value 920/473 exactly 2 time(s)',
    'the sector-bridge value': 'the doc carries the parent value −56/165 exactly 1 time(s)',
    'a step check count': 'step L0: the doc reports 128 checks, and the log has 128 passing and 0 failing',
    "the dated correction's count": "the control script's check count is recorded truthfully",
    'the filed row rewritten, correction intact': "the control script's check count is recorded truthfully",
    'the correction note undated': "the control script's check count is recorded truthfully",
    'a G2 bound': "G2 constant 'TrN̄² ≤ 171/2' is in the frozen proof",
    'the G2 invariant form': "G2 constant '−5/231 − |f|²/22 + (7/11)|a₀₀|² + TrN̄²/198' is in the frozen pro",
    'a higher-dimensional space dropped': 'the spaces of projective dimension at least 2 are the logged four',
    'a point class dropped': 'the six point classes are the logged ones',
    'a parent value dropped from P1': 'P1 carries all eight parent values',
    'P1 miscited': 'P1 sources the weight-state values to the paper Section 7.2',
    'the § 7.2 source row': 'the sources table names Section 7.2 and what it records',
    'L1 out of the audit set': 'the audit set includes L1: an adversarial audit follows, which al...',
    'L1 out of the definition of done': 'the audit set includes L1: The audit grades L0, L1, L2, G1 and G2...',
    'the retired monotonicity inference': 'D2 states the deferred indices are unknown',
    'the network posture': 'the firewall states the network posture',
    'the exposed items misstated': 'it names the items a search could reach',
    "H2's monic rule": 'H2 scores after monic normalization',
    'the bridge arithmetic': 'and the doc states 4 L_T as −42/55',
    'the control value': 'N1 states the control value, and the log measured it',
    'the control gradient': 'the gradient norm N1 states is the one s1_control.py measured',
    'the control residual, altered but still nonzero': 'the residual N1 freezes at -i Jx u  is the one s1_control.py measured',
    "N1's pass condition, scoring one residual by sign": 'N1 scores all six quantities and separates a mismatch from a control that cann',
    "N1's fail condition weakened": 'N1 scores all six quantities and separates a mismatch from a control that cann',
    'one frozen residual, still nonzero': 'the residual N1 freezes at -i Jy u  is the one s1_control.py measured',
    'the sum rule': 'N1 and the control log carry the same sum rule',
    'the matrix renamed back to the second variation': 'N1 names the matrix rather than the second variation, and says why they agree',
    'P2 chart device loosened': 'P2 carries the agreed pass condition exactly, with nothing added or dropped',
    'L1 chart device loosened': 'L1 carries the agreed pass condition exactly, with nothing added or dropped',
    'item 11 back to a yes or no': 'worklist item 11 requests all six scored quantities, un-normalized, naming the',
    'N1 drifting from the point the worklist supplies': 'N1 and worklist item 11 name the same control point',
    'a pipe left unescaped in a table': 'P2 carries the agreed pass condition exactly, with nothing added or dropped',
    'the worklist point drifting from N1': 'N1 and worklist item 11 name the same control point',
    'frozen values leaked into the handout': 'no frozen value, count, name or representative leaks into the worklist (hits:',
    'the D7 input': 'the worklist gives D7 as an input, with both weights',
    'the zeros rule': 'the worklist carries the zeros rule',
    'a count of the answer leaked': "every numeral in the worklist is a reviewed one (new: {'ten': 1}, gone: {})",
    'a census vector used as the control point': 'N1 and worklist item 11 name the same control point',
    'item 5b renamed': 'the worklist asks for the in-locus second variations, with the off-diagonal en',
    'the cross term dropped': 'the worklist asks for the in-locus second variations, with the off-diagonal en',
    'the sector-scaled numbers dropped': 'the worklist asks for the sector-scaled numbers',
    'the Clebsch-Gordan convention unpinned': 'the worklist pins the Clebsch-Gordan phase convention',
    'the separation clause weakened': 'the worklist forces the separation of same-stabilizer classes, without naming',
    'the hard pair named in the handout': 'the worklist forces the separation of same-stabilizer classes, without naming',
    '5b normalizing every projection again': 'the worklist asks for the in-locus second variations, with the off-diagonal en',
    'the zero direction left unidentified': 'the worklist handles the tangent whose projection vanishes, rather than normal',
    'the orbit-null control unnamed': 'the worklist asks for the in-locus second variations, with the off-diagonal en',
    'the zero projection normalized': 'the worklist asks for the in-locus second variations, with the off-diagonal en',
    "the control's own arm failing": "the control script's check count is recorded truthfully",
    'the frozen S0 proof edited': 'the frozen S0 maximum proof is unchanged',
}

PATHS = {'task': 'prereg/s1_task_details.md', 'work': 'prereg/worklist.md',
         'log': 's1_control_log.txt', 's0': 'S0_S3_MAXIMUM.md', 'inventory': 'prereg/gate_inventory.txt'}


def build(root):
    return subprocess.run([sys.executable, 's1_prereg_build.py'], cwd=os.path.join(root, 'prereg'),
                          capture_output=True, text=True)


def stage(root):
    os.makedirs(os.path.join(root, 'prereg'), exist_ok=True)
    for f in SUPPORT:
        shutil.copy(os.path.join(M8S, f), os.path.join(root, f))
    # every file the gate reads, including its own pinned reference data, or the staged parent is not the real one
    for f in ('s1_task_details.md', 'worklist.md', 's1_prereg_build.py', 'gate_inventory.txt', 's1_prereg_arms.py'):
        shutil.copy(os.path.join(HERE, f), os.path.join(root, 'prereg', f))


root = tempfile.mkdtemp(prefix='s1_arms_')
caught = missed = 0
try:
    stage(root)
    r = build(root)
    print('parent: ' + ('GREEN' if r.returncode == 0 else 'RED, so no arm below proves anything'), flush=True)
    if r.returncode != 0:
        print(r.stdout[-2000:])
        sys.exit(1)
    for i, (which, a, b, what) in enumerate(ARMS, 1):
        stage(root)  # every arm starts from the green parent
        p = os.path.join(root, PATHS[which])
        if a == '<DELETE FILE>':
            os.remove(p)
            out = build(root)
            hit = [ln for ln in out.stdout.split('\n') if ln.startswith('FAIL')]
            want = EXPECT.get(what)
            ok = hit and (not want or any(want in ln for ln in hit))
            caught, missed = (caught + 1, missed) if ok else (caught, missed + 1)
            print(f'{i:3d} {"caught    " if ok else "NOT CAUGHT"} {what}' + (f'  ->  {hit[0][5:].strip()[:70]}' if hit else ''))
            continue
        s = open(p, encoding='utf-8').read()
        if a not in s:
            print(f'{i:3d} UNAPPLIED  {what}: the text to mutate was not found')
            missed += 1
            continue
        open(p, 'w', encoding='utf-8').write(s.replace(a, b, 1))
        out = build(root)
        hit = [ln for ln in out.stdout.split('\n') if ln.startswith('FAIL')]
        want = EXPECT.get(what)
        if not hit:
            missed += 1
            # a mutation that breaks the build produces no gate failure at all, which is a broken arm, not a catch
            print(f'{i:3d} {"BROKE THE BUILD" if out.returncode else "NOT CAUGHT"} {what}')
        elif want and not any(want in ln for ln in hit):
            missed += 1
            print(f'{i:3d} WRONG GATE {what}: expected {want[:44]!r}, tripped {hit[0][5:].strip()[:44]!r}')
        else:
            caught += 1
            print(f'{i:3d} caught     {what}  ->  {hit[0][5:].strip()[:78]}')
    stage(root)
    print('restored: ' + ('GREEN' if build(root).returncode == 0 else 'RED'), flush=True)
finally:
    shutil.rmtree(root, ignore_errors=True)
print(f'{caught} of {len(ARMS)} arms caught' + ('' if not missed else f', {missed} NOT CAUGHT'))
sys.exit(1 if missed else 0)
