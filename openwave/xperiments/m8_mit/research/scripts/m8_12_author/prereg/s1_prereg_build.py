"""Gates the S1 pre-registration: the task doc against the three step logs, and the worklist against leakage.

Every frozen value in the task doc must appear in the log of the step that produced it; the worklist must contain no
value, count or representative that the blind agents are meant to derive; the house rules (no em-dash, the short
forms, the firewall list, the disclosure) must hold; and the budget must be inside the Section 12.2 cap.
Each gate that can be armed is armed with a planted violation.
Usage: python3 s1_prereg_build.py
"""
import hashlib
import os
import re
import sys
from collections import Counter

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
M8S = os.path.dirname(HERE)
doc = open(os.path.join(HERE, 's1_task_details.md'), encoding='utf-8').read()
work = open(os.path.join(HERE, 'worklist.md'), encoding='utf-8').read()
logs = {n: open(os.path.join(M8S, f's1_{n}_log.txt'), encoding='utf-8').read() for n in ('L0', 'L1', 'H')}
cells = lambda row: re.split(r'(?<!\\)\|', row)  # a table cell can carry an escaped pipe, as the chart rows do
prose = lambda s: re.sub(r'`[^`]*`', ' ', s)  # the prose gates read prose, not the code spans, where I is the identity
unquoted = lambda s: s.replace('"stable"', ' ').replace('"unstable"', ' ')  # the rule may name what it forbids
fails = []


names = []


def gate(name, ok):
    names.append(name)
    print(('PASS ' if ok else 'FAIL ') + name)
    if not ok:
        fails.append(name)


# Block edits of this script have twice deleted gates silently, found only by the mutation suite ten minutes later.
# The run's gate names are pinned in gate_inventory.txt, normalized so that a changed count or value is not a change
# of gate. Any gate added or removed fails here, and re-pinning is deliberate.
INVENTORY = os.path.join(HERE, 'gate_inventory.txt')
norm = lambda n: re.sub(r'\([^)]*\)', '()', n).strip()  # parenthetical detail only: erasing digits collapsed 24
# gates onto shared names, which let a substitution inside a group pass. Digits in a gate name are identity, an
# orbit or a claim row or a step's check count, and a pin should report them when they move (xrodz, #571).


print('(A) the task doc')
no_em = lambda s: '—' not in s
gate('  no em-dash, and a planted one is refused', no_em(doc) and not no_em(doc.replace('## TASK PLANNING', '## TASK PLANNING —', 1)))
words = len(doc.split())
gate(f'  budget: inside the Section 12.2 cap of 8000 ({words} words)', words <= 8000)
first = lambda t: re.search(r'(?<![A-Za-z])(I|my|My)(?![A-Za-z])', prose(t))
gate('  no first-person voice outside the code spans', not first(doc))
gate('  arm: a planted first person is caught', bool(first(doc + ' I checked this by hand.')))
bare = lambda s: re.search(r'\b(un)?stable\b', unquoted(s))
gate('  the short-form rule is stated, and no bare stable or unstable is used',
     'Never "stable" or "unstable" alone' in doc and not bare(doc))
gate('  arm: a planted bare unstable is caught', bool(bare(doc + ' the branch is unstable')))

# the census rows against the step-3 log, field by field against its summary table, never by substring
CEN = {}
for ln in logs['H'].split('\n'):
    m = re.match(r'^  (coherent v3|v2|v1|zonal v0|octahedron|hexagon|pyramid|prism|C3 ray|D2 ray) +([\d/]+) +(\d+) +\((\d+), (\d+), (\d+)\) +(\d+) +(\d+)$', ln)
    if m:
        CEN[m.group(1)] = (m.group(2), m.group(3), f'({m.group(4)}, {m.group(5)}, {m.group(6)})', m.group(7), m.group(8))
gate(f'  the step-3 log carries a census summary of ten orbits (found {len(CEN)})', len(CEN) == 10)
NAME = {'coherent `v₃`': 'coherent v3', '`v₂`': 'v2', '`v₁`': 'v1', 'zonal `v₀`': 'zonal v0', 'octahedron': 'octahedron',
        'hexagon': 'hexagon', 'pyramid': 'pyramid', 'prism': 'prism', 'C3 ray': 'C3 ray', 'D2 ray': 'D2 ray'}
rows = re.findall(r'^\| (coherent `v₃`|`v₂`|prism|`v₁`|D2 ray|C3 ray|pyramid|octahedron|zonal `v₀`|hexagon) \| .*?\| `?([\d/]+)`? \| (\d+) \| \((\d+), \*?\*?(\d+)\*?\*?, (\d+)\) \| (\d+) \| (\d+) \|$', doc, re.M)
gate(f'  the census table in the doc has ten rows (found {len(rows)})', len(rows) == 10)
seen = set()
for name, val, k, a_, b_, c_, ip, im in rows:
    key = NAME[name]
    seen.add(key)
    got = (val, k, f'({a_}, {b_}, {c_})', ip, im)
    gate(f'    {key:14s} value, dim, signature and both indices match the log row exactly', CEN.get(key) == got)
gate('  every logged orbit appears in the doc', seen == set(CEN))
gate('  arm: one altered field is caught', CEN['hexagon'] != ('464', '9', '(9, 0, 0)', '9', '0'))
# the representative column against the ORB dictionary the step-3 script actually ran on
src = open(os.path.join(M8S, 's1_H.py'), encoding='utf-8').read()
ns = {'sp': sp}
exec(src[src.index('ORB = {'):src.index('}\n', src.index("'D2 ray'")) + 1], ns)
SUB = {'\u2083': 3, '\u2082': 2, '\u2081': 1, '\u2080': 0, '\u208b\u2081': -1, '\u208b\u2082': -2, '\u208b\u2083': -3}


def parse_rep(text):
    """The doc's representative, as the {m: coefficient} dictionary the script uses."""
    out = {}
    for term in text.split(' + '):
        co, _, sub = term.partition('v')
        co = co.rstrip('\u00b7').strip()
        i = sp.I if co.startswith('i') else 1
        co = co[1:] if co.startswith('i') else co
        if co.startswith('\u221a'):
            body = co[1:].strip('()')
            c = sp.sqrt(sp.Rational(body) if '/' in body else sp.Integer(body))
        else:
            c = sp.Integer(co) if co else sp.Integer(1)
        out[SUB[sub]] = sp.simplify(i * c)
    return out


census_block = doc[doc.index('### The census'):doc.index('###', doc.index('### The census') + 5)]
reps = re.findall(r'^\| (coherent `v\u2083`|`v\u2082`|prism|`v\u2081`|D2 ray|C3 ray|pyramid|octahedron|zonal `v\u2080`|hexagon) \| `([^`]+)` \|', census_block, re.M)
gate(f'  the census table gives a representative for every orbit ({len(reps)})', len(reps) == 10)
for name, rep in reps:
    want = ns['ORB'][NAME[name]]
    got = parse_rep(rep)
    gate(f'    {NAME[name]:14s} representative {rep!r} is the vector the step-3 script used',
         set(want) == set(got) and all(sp.simplify(want[m] - got[m]) == 0 for m in want))
gate('  arm: an altered representative is caught', parse_rep('v\u2082 + 3\u00b7v\u208b\u2081') != ns['ORB']['C3 ray'])

# and every signature written anywhere in the doc, claim tables included, is one the log carries
sigs = set(re.findall(r'\((\d+), \*?\*?(\d+)\*?\*?, (\d+)\)', doc))
logsigs = {(x[2][1:-1].split(', ')[0], x[2][1:-1].split(', ')[1], x[2][1:-1].split(', ')[2]) for x in CEN.values()}
gate(f'  every signature triple in the doc is a logged one ({len(sigs)} distinct)', sigs <= logsigs)
gate('  arm: a signature triple the log does not carry is caught', not {('9', '1', '0')} <= logsigs)

# the seven lines: the doc's own rows parsed and compared with the step-2 log, not fragments looked up in it
sx, sy, ss = sp.symbols('x y s', real=True)
SUBS = str.maketrans('\u2080\u2081\u2082\u2083\u2084\u2085\u2086', '0123456')


def norm_line_label(t):
    return t.replace('\u2083', '3').replace('\u2082', '2').replace('\u2084', '4').replace('\u2085', '5').replace('\u2086', '6') \
            .replace('`', '').replace(' ', '').replace('v\u208b', 'v-').replace('\u2081', '1').replace('\u2080', '0').replace('v\u2083', 'v3')


def to_expr(t):
    """The doc's typeset restriction, as the expression the step-2 script printed."""
    t = t.replace('\\|', '|').replace('\u2212', '-').replace('\u00b7', '*')  # the table cells escape their pipes
    t = t.replace('|z|\u2074', '(x**2 + y**2)**2').replace('|z|\u00b2', '(x**2 + y**2)')
    t = t.replace('\u2074', '**4').replace('\u00b3', '**3').replace('\u00b2', '**2')
    t = re.sub(r'(\d)\s*([xys(])', r'\1*\2', t)
    t = re.sub(r'\)\s*([xys(])', r')*\1', t)
    return sp.sympify(t, locals={'x': sx, 'y': sy, 's': ss})


LOGQ, LOGV = {}, {}
lab = None
for ln in logs['L1'].split('\n'):
    m = re.match(r'^PASS   (C\d \{[^}]+\}):', ln)
    if m:
        lab = m.group(1).replace(' ', '')
    m = re.match(r'^    q\(s\) = (.+)$', ln)
    if m and lab:
        LOGQ[lab] = sp.sympify(m.group(1), locals={'s': ss})
    m = re.match(r'^    vertex s\* = (\S+)', ln)
    if m and lab:
        LOGV[lab] = sp.Rational(m.group(1))
gate(f'  the step-2 log carries five cyclic quadratics with their vertices ({len(LOGQ)})', len(LOGQ) == 5 == len(LOGV))

rowsL = re.findall(r'^\| (C[\u2083\u2084\u2085\u2086] `\{[^`]+\}`) \| `([^`]+)` \| ([^|]+) \|$', doc, re.M)
gate(f'  the doc freezes five cyclic lines ({len(rowsL)})', len(rowsL) == 5)
for label, poly, crit in rowsL:
    key = norm_line_label(label)
    ok = key in LOGQ and sp.expand(to_expr(poly) - LOGQ[key]) == 0
    gate(f'    {key:14s} restriction matches the step-2 log exactly', ok)
    m = re.search(r's\* = ([\d/]+)', crit)
    want = LOGV.get(key)
    gate(f'    {key:14s} vertex matches, interior or endpoint as logged',
         (sp.Rational(m.group(1)) == want and 0 < want < 1) if m else (want in (0, 1) and 'endpoint' in crit))
gate('  arm: an altered restriction is caught', sp.expand(to_expr('\u2212(125/132)s\u00b2 + (10/11)s + 4/77') - LOGQ['C5{v3,v-2}']) != 0)
gate('  arm: an altered vertex is caught', sp.Rational('13/25') != LOGV['C5{v3,v-2}'])

LOGF, LOGZ = {}, {}
lab = None
for ln in logs['L1'].split('\n'):
    m = re.match(r'^    f\(x, y\) = (.+)$', ln)
    if m:
        pend = sp.sympify(m.group(1), locals={'x': sx, 'y': sy})
    m = re.match(r'^PASS   (D\d \{[^}]+\}): Poincare', ln)
    if m:
        lab = m.group(1).replace(' ', '')
        LOGF[lab] = pend
        LOGZ[lab] = zs
        zs = []
    m = re.match(r'^    z = (\S+) \+ i\*\((\S+)\):', ln)
    if m:
        zs = (zs if 'zs' in dir() else []) + [sp.sympify(m.group(1)) + sp.I * sp.sympify(m.group(2))]
gate(f'  the step-2 log carries both dihedral charts with their critical sets ({len(LOGF)})', len(LOGF) == 2 and all(len(v) == 5 for v in LOGZ.values()))

rowsD = re.findall(r'^\| (D[\u2082\u2083] `\{[^`]+\}`) \| `([^`]+)` \| ([^|]+) \|$', doc, re.M)
gate(f'  the doc freezes two dihedral lines ({len(rowsD)})', len(rowsD) == 2)
for label, f_, crit in rowsD:
    key = norm_line_label(label).replace('+', '+')
    ok = key in LOGF and sp.simplify(to_expr(f_) - LOGF[key]) == 0
    gate(f'    {key:16s} chart formula matches the step-2 log exactly', ok)
    pts = set()
    for tok in crit.split(';'):
        tok = tok.strip().strip('`').replace('z = ', '').replace('\u00b1', '')
        if tok in ('\u221e', '`\u221e`', ''):
            continue
        t = tok.replace('`', '').replace('\u221a', 'sqrt')
        t = re.sub(r'sqrt(\d+)', r'sqrt(\1)', t)
        im = t.startswith('i')
        v = sp.sympify(t[1:] if im else t)
        pts |= {sp.I * v, -sp.I * v} if im else ({v, -v} if v != 0 else {sp.Integer(0)})
    gate(f'    {key:16s} chart critical set matches the logged one', pts == set(LOGZ.get(key, [])))
gate('  arm: an altered chart point is caught', {sp.sqrt(231) / 10, -sp.sqrt(231) / 10} != set(LOGZ['D3{v3+v-3,v0}']))
gate('  the infinite chart point is named on both dihedral rows', all('\u221e' in r[2] for r in rowsD))

# every restatement of a restriction anywhere in the doc, claim rows included, not only the frozen table
spans = [t for t in re.findall(r'`([^`]+)`', doc.replace('\\|', '|')) if re.search(r's\u00b2|\|z\|\u2074', t) and 'r\u0302' not in t]
known = list(LOGQ.values()) + list(LOGF.values())
gate(f'  the doc restates {len(spans)} restrictions, and each is one the step-2 log printed',
     len(spans) >= 9 and all(any(sp.simplify(to_expr(t) - k) == 0 for k in known) for t in spans))
gate('  arm: an altered restatement is caught',
     not any(sp.simplify(to_expr('(100|z|\u2074 \u2212 20x\u00b2 + 148y\u00b2 + 461)/(231(|z|\u00b2 + 2)\u00b2)') - k) == 0 for k in known))

# the filed document's own structure
HEADS = ['## TASK PLANNING', '## THE LAW, PINNED', '## SETTING AND DERIVATION', '## THE FIREWALL', '## DISCLOSURE',
         '## CANDIDATE PRE-REGISTERED CLAIMS', '## FROZEN VALUES', "## FEASIBILITY, AND THE AUTHOR'S DERIVATION",
         '## TO BE FIXED AT GO', '## DEFINITION OF DONE']
missing = [h for h in HEADS if h not in doc]
gate(f'  every required section heading is present (missing: {missing})', not missing)
gate('  arm: a renamed heading is caught', '## THE WALL' not in HEADS)

# the characteristic polynomials: every coefficient of three digits or more must appear in the step-3 log
for line in doc.split('\n'):
    m = re.match(r'^\| (coherent `v₃`|`v₂`|`v₁`|zonal `v₀`|octahedron|hexagon|pyramid|prism|C3 ray|D2 ray) \| `\((.+)\)?`? \|$', line)
    if m and 'λ' in line:
        big = [c for c in re.findall(r'\d{3,}', line)]
        gate(f'    charpoly of {m.group(1).strip("`"):14s}: coefficients {big} are in the step-3 log', all(c in logs['H'] for c in big))

# the classification against the step-1 log
for cls in ('C3 {v2,v-1}', 'C4 {v3,v-1}', 'C4 {v2,v-2}', 'C5 {v3,v-2}', 'C6 {v3,v-3}', 'D3 {v3+v-3,v0}', 'D2 {v2+v-2,v0}'):
    gate(f'    line class {cls!r} is in the step-1 log', cls in logs['L0'])
gate('  the four spaces of dimension at least 3 are the logged ones',
     "[('C1', 7), ('C2', 3), ('C2', 4), ('C3', 3)]" in logs['L0'] and 'four such, from three groups' in doc)

# the classification section against the step-1 log
gate('  the six point classes are the logged ones',
     "exactly the six canonical points (['hexagon', 'octahedron', 'v0', 'v1', 'v2', 'v3'])" in logs['L0']
     and all(t in doc[doc.index('### The classification'):] for t in
             ('`v\u2083`', '`v\u2082`', '`v\u2081`', '`v\u2080`', '`v\u2082 + v\u208b\u2082`', '`v\u2083 + v\u208b\u2083`')))
gate('  the seven line classes are the logged ones', 'exactly the seven canonical lines (7)' in logs['L0'])
gate('  the spaces of projective dimension at least 2 are the logged four',
     "[('C1', 7), ('C2', 3), ('C2', 4), ('C3', 3)]" in logs['L0']
     and '`C\u2081` (7), `C\u2082` (4 and 3), `C\u2083` (3)' in doc)
gate('  arm: a fifth such space is not in the log', "('C4', 2)]" not in logs['L0'])

# P4 and P5, the M8.11 parents, against the step-3 log that reproduced them
for val, where in (('-104/55', 'pyramid: H along the unit C5-line tangent'), ('920/473', 'prism: H along the chart tangent v0'),
                   ('8/11', 'prism: H along the chart tangent i v0'), ('-56/165 and -21/110', 'the sector bridge gives M8.11 L_T')):
    gate(f'    M8.11 parent {val} reproduces in the step-3 log', any(where in ln and val in ln for ln in logs['H'].split('\n')))
PARENT_N = {'−104/55': 2, '920/473': 2, '8/11': 2, '−56/165': 1, '−21/110': 1}  # every copy is counted, so altering one of two is caught
for val, n in PARENT_N.items():
    gate(f'    the doc carries the parent value {val} exactly {n} time(s)', doc.count(val) == n)
gate('  arm: altering one copy of a parent value is caught', doc.replace('\u2212104/55', '\u2212104/56', 1).count('\u2212104/55') != PARENT_N['\u2212104/55'])

# G1 and G2, whose constants come from the frozen S0 proof rather than from this task's steps
s3 = open(os.path.join(M8S, 'S0_S3_MAXIMUM.md'), encoding='utf-8').read()
gate('  the frozen S0 maximum proof is unchanged',
     hashlib.sha256(s3.encode('utf-8')).hexdigest() == '7c634a33fdb0ee05d5a934395a414f6346ba26e049b07f2619961792ab725498')
for t in ('\u22125/231 \u2212 |f|\u00b2/22 + (7/11)|a\u2080\u2080|\u00b2 + TrN\u0304\u00b2/198', '|a\u2080\u2080|\u00b2 \u2264 1/7', 'TrN\u0304\u00b2 \u2264 171/2'):
    gate(f'    G2 constant {t!r} is in the frozen proof', t in s3 and t.replace('|', '\\|') in doc or t in s3 and t in doc)
gate('  the frozen proof closes at the hexagon value', '\u22125/231 + 1/11 + 171/396 = 463/924' in s3)
gate('  arm: an altered bound is not in the frozen proof', 'TrN\u0304\u00b2 \u2264 171/3' not in s3)

for need, label in ((('S0 memo', 'lemma note', 'M8_S/', 'M8_DYNAMICS/', "M8.11's records"), 'firewall names every excluded source'),
                    (('DISCLOSURE_INVENTORY.md', '57 records', 'A dated correction', '215.775792502', '238.016528926', '238.356854743'), 'disclosure carries the inventory, the correction and the deferred orbits'),
                    (('Group P', 'Group L', 'Group H', 'Group G', 'Negative controls', 'Diagnostics'), 'all claim groups are present')):
    gate('  ' + label, all(s in doc for s in need))
claims_block = doc[doc.index('## CANDIDATE PRE-REGISTERED CLAIMS'):doc.index('## FROZEN VALUES')]
CHART = {
    'P2': "reproduced exactly in the paper's chart, `b\u2081 = v\u2083 + v\u208b\u2083` and `b\u2082 = v\u2080` unnormalized, "
          "or an equivalent restriction in the chart the room states, carried onto it by the room's own stated change of chart",
    'L1': "every row reproduces exactly, a dihedral row either in the paper's chart or as an equivalent restriction the room "
          "carries onto it by its own stated change of chart, and completeness rests on an elimination argument rather than on a solver",
}
for cid, clause in CHART.items():
    row = [ln for ln in doc.split('\n') if ln.startswith(f'| {cid} |')][0]
    gate(f'  {cid} carries the agreed pass condition exactly, with nothing added or dropped', cells(row)[4].strip() == clause)
gate('  arm: a loosened chart clause is caught', CHART['P2'] + ', or otherwise' != CHART['P2'])
gate('  arm: a dropped chart clause is caught', 'reproduced exactly' != CHART['P2'])
gate('  the exposure table states what was visible per orbit', doc.count('| orbit | already visible | new here |') == 1)

claim_rows = [cells(ln) for ln in claims_block.split('\n') if re.match(r'^\| (P\d|L\d|H\d|G\d|N\d) \|', ln)]
diag_rows = [cells(ln) for ln in claims_block.split('\n') if re.match(r'^\| D\d \|', ln)]
gate(f'  each diagnostic is recorded, with no pass or fail column ({len(diag_rows)} rows)',
     len(diag_rows) == 2 and all(len(r) == 4 and r[2].strip() for r in diag_rows))
gate(f'  every claim row has a standing, a pass and a fail condition ({len(claim_rows)} rows)',
     len(claim_rows) >= 15 and all(len(r) >= 6 and all(x.strip() for x in r[1:6]) for r in claim_rows))
# Provenance, verified in the repo on 2026-09-19: surviving-ray.md Section 7.2 records that every weight state is
# critical and that the seven realize 1, 36, 225 and 400, so P1 claims all eight values, and no row claims the
# weight-state values as new. M8.1.2 D1 carries 1, 400, 288, 463; D3 the pyramid; D4 the prism.
p1 = [ln for ln in doc.split('\n') if ln.startswith('| P1 |')]
gate('  the doc has a P1 row', len(p1) == 1)
gate('  P1 carries all eight parent values', all(re.search(r'(?<![\d/])' + v + r'(?![\d/])', p1[0]) for v in ('1', '36', '225', '400', '288', '463', '1188/5', '8800/43')))
gate('  P1 sources the weight-state values to the paper Section 7.2', '\u00a7 7.2' in p1[0])
gate('  P1 keeps the other parents on D1, D3 and D4', all(t in p1[0] for t in ('D1', 'D3', 'D4', '\u00a7 5.8')))
gate('  the sources table names Section 7.2 and what it records',
     any(ln.startswith('| The fourth bedrock paper, \u00a7 7.2 |') and '1, 36, 225 and 400' in ln for ln in doc.split('\n')))
gate('  no row claims a weight-state value as new', not [ln for ln in doc.split('\n') if 'No parent records' in ln])
gate('  arm: a row claiming one as new is caught', 'No parent records' in doc + ' No parent records either value')

# the audit grades L1, whose completeness is an argument, in both places that name the audit set
for where in ('an adversarial audit follows, which also grades the arguments L0, L1, L2, G1 and G2',
              'The audit grades L0, L1, L2, G1 and G2 as arguments'):
    gate(f'  the audit set includes L1: {where[:38]}...', where in doc)
gate('  arm: an audit set without L1 is caught', 'grades the arguments L0, L2' not in doc)

# D2 records the index sequence over the classified set and infers nothing about the deferred orbits
d2 = [ln for ln in doc.split('\n') if ln.startswith('| D2 |')]
gate('  D2 confines the sequence to the classified orbits', len(d2) == 1 and 'Among the ten classified orbits' in d2[0])
gate('  D2 states the deferred indices are unknown', "deferred orbits' indices are unknown and may break this ordering" in d2[0])
gate('  no monotonicity inference about the deferred orbits survives anywhere',
     'do not fill' not in doc and 'cannot fill' not in doc and 'monotone in the value' not in doc)
gate('  arm: the retired inference is caught', 'do not fill' in doc + ' so they do not fill the gaps')

# the network posture, which decides how items 1 to 4 and 8 are scored
posture = lambda t: 'The run is therefore offline' in t and 'reproduced or located' in t
gate('  the firewall states the network posture', posture(doc))
gate('  it names the items a search could reach', 'items 1 to 4' in doc and 'item 8' in doc and 'items 5 to 7, are new' in doc)
gate('  arm: a firewall with the posture cut is caught', not posture(doc.replace('The run is therefore offline', 'The run may be online', 1)))

# H2 is scored monic, because the frozen table prints each polynomial up to a positive constant
h2 = [ln for ln in doc.split('\n') if ln.startswith('| H2 |')]
gate('  H2 scores after monic normalization', len(h2) == 1 and 'monic normalization' in h2[0] and 'up to a positive constant' in h2[0])
gate('  the spectra table says the same', 'characteristic polynomial, up to a positive constant' in doc)

# P5, recomputed here rather than carried: what the rooms report is w6 times H, which is 4 L_T
p5 = [ln for ln in doc.split('\n') if ln.startswith('| P5 |')][0]
H_pyr = sp.Rational(-104, 55)
for w, lt, four in ((sp.Rational(28, 39), sp.Rational(-56, 165), '\u2212224/165'), (sp.Rational(21, 52), sp.Rational(-21, 110), '\u221242/55')):
    gate(f'    the bridge closes in the sector with w6 = {w}: (w6/4)H = {lt}', sp.simplify(w / 4 * H_pyr - lt) == 0)
    gate(f'    and the doc states 4 L_T as {four}', four in p5 and sp.simplify(sp.Rational(four.replace('\u2212', '-')) - 4 * lt) == 0)
gate('  arm: the old value 21/55 does not satisfy the bridge', sp.simplify(sp.Rational(-21, 55) - 4 * sp.Rational(-21, 110)) != 0)

# the worklist asks for the quantities P4 and P5 are scored on, and pins the Clebsch-Gordan phase
gate('  the worklist asks for the in-locus second variations, with the off-diagonal entries',
     '**5b.**' in work and 'every diagonal entry and every off-diagonal entry' in work and 'Gram matrix' in work
     and 'whether its projection onto `N_u` is zero' in work and 'as an orbit-null control' in work and 'Normalize the rest in `Re\u27e8\u00b7,\u00b7\u27e9`' in work)
gate('  the worklist handles the tangent whose projection vanishes, rather than normalizing zero',
     'Where it is zero, say which direction of `O_u` it is' in work and 'Normalize the rest' in work)
gate('  arm: the old wording, which normalized every projection, is caught',
     'Normalize the rest' not in work.replace('Normalize the rest', 'normalize them all', 1))
gate('  the worklist asks for the sector-scaled numbers', 'multiplied by each `w\u2086`' in work and 'items 4, 5, 5b and 6' in work)
gate('  the worklist pins the Clebsch-Gordan phase convention',
     'Condon-Shortley' in work and '\u27e83 3; 3 3 | 6 6\u27e9 = +1' in work)
sep = lambda t: 'Wherever two of them have the same stabilizer, give that separation exactly' in t and 'C\u2084' not in t
gate('  the worklist forces the separation of same-stabilizer classes, without naming the pair', sep(work))
gate('  arm: the weaker two-classes wording is caught', not sep(work.replace('Wherever two of them have the same stabilizer, give that separation exactly', 'Show that two of them differ', 1)))
gate('  arm: naming the hard pair in the worklist is caught', not sep(work + ' the two C\u2084 lines'))

# N1: the frozen control, the point the worklist hands the rooms, and the measured log must be one and the same
ctl = open(os.path.join(M8S, 's1_control_log.txt'), encoding='utf-8').read()
n1 = [ln for ln in doc.split('\n') if ln.startswith('| N1 |')][0]
w11 = [ln for ln in work.split('\n') if ln.startswith('**11.**')][0]
gate('  worklist item 11 requests all six scored quantities, un-normalized, naming the matrix',
     all(t in w11 for t in ('924\u00b7r\u0302\u2086', 'norm of the tangential gradient', '\u2016M_u d\u2016',
                            'four generators', 'without normalizing', 'M_u = Hess N(u) \u2212 4N(u)\u00b7I')))
gate('  arm: the earlier yes-or-no wording is caught',
     'whether the directions of `O_u` are annihilated' not in w11)
v_doc = re.search(r'`u = \(([^`]+)\)/\u221a3`', n1)
v_work = re.search(r'`u = \(([^`]+)\)/\u221a3`', w11)
gate('  N1 and worklist item 11 name the same control point',
     bool(v_doc) and bool(v_work) and v_doc.group(1) == v_work.group(1))
gate('  arm: a different control point in either file is caught', v_doc.group(1) != 'v\u2082 + v\u208b\u2081')
gate('  the control log measured that same point',
     'control u = (v3 + v1 - v-2)/sqrt(3)' in ctl and v_doc.group(1) == 'v\u2083 + v\u2081 \u2212 v\u208b\u2082')
gate('  the control log is green, arms included',
     ctl.rstrip().endswith('ALL PASS') and 'at a critical point the gradient vanishes' in ctl
     and 'at a critical point every orbit residual vanishes exactly' in ctl)
src_c = open(os.path.join(M8S, 's1_control.py'), encoding='utf-8').read()
nsc = {'sp': sp}
exec(src_c[src_c.index('FROZEN = {'):src_c.index('}\n', src_c.index("'-i Jz u'")) + 1], nsc)


def surd(t):
    """A residual as N1 typesets it, as the expression the control script compared against."""
    t = t.replace('\u2212', '-').replace('\u221a(', 'sqrt(')
    t = re.sub(r'\u221a(\d+)', r'sqrt(\1)', t)
    t = re.sub(r'(\d)\s*sqrt', r'\1*sqrt', t)
    return sp.sympify(t)


res = re.findall(r'`([^`]*\u221a[^`]*)` at `(i\u00b7u|\u2212i\u00b7J_[xyz] u)`', n1)
gate(f'  N1 freezes a residual at each of the four generators ({len(res)})', len(res) == 4)
KEY = {'i\u00b7u': 'i u', '\u2212i\u00b7J_x u': '-i Jx u', '\u2212i\u00b7J_y u': '-i Jy u', '\u2212i\u00b7J_z u': '-i Jz u'}
for expr, where in res:
    k = KEY[where]
    gate(f'    the residual N1 freezes at {k:8s} is the one s1_control.py measured',
         sp.simplify(surd(expr) - nsc['FROZEN'][k]) == 0)
gate('  arm: a neighbouring surd would not match', sp.simplify(surd('2\u221a3685/297') - nsc['FROZEN']['-i Jz u']) != 0)
gate('  N1 and the control log carry the same sum rule',
     '17368/29403' in n1 and 'the four squares sum to 17368/29403' in ctl)
gate('  N1 names the matrix rather than the second variation, and says why they agree on the orbit',
     'M_u = Hess N(u) \u2212 4N(u)\u00b7I' in n1 and 'orthogonal to the orbit' in n1
     and 'grad N is orthogonal to every orbit direction' in ctl)
gn = re.search(r'has norm `([^`]+)`', n1)
gate('  N1 states the tangential gradient norm', bool(gn))
gate('    the gradient norm N1 states is the one s1_control.py measured',
     bool(gn) and sp.simplify(surd(gn.group(1)) - nsc['FROZEN']['i u']) == 0 and '2*sqrt(1002)/297' in ctl)
gate('  arm: a mutated gradient norm is caught, even though the same surd appears again as a residual',
     sp.simplify(surd('2\u221a1003/297') - nsc['FROZEN']['i u']) != 0)
val_doc = re.search(r'`924\u00b7r\u0302\u2086 = (\d+)`', n1)
gate('  N1 states the control value, and the log measured it',
     bool(val_doc) and f'924*r6 = {val_doc.group(1)} ' in ctl)
gate('  arm: a frozen control value the log does not carry is caught', '2*sqrt(1003)/297' not in ctl)
def scores_all_six(row):
    """N1 is only a control if a wrong but nonzero residual fails it, and a claim mismatch is not the same event
    as the control ceasing to discriminate, which needs every residual and the gradient to vanish together."""
    c = cells(row)
    return ('all six reported quantities reproduce' in c[4] and 'five are independent' in c[4]
            and 'all four residuals are nonzero' in c[4]
            and 'any of the six differs' in c[5] and 'ceases to discriminate only if' in c[5])


gate('  N1 scores all six quantities and separates a mismatch from a control that cannot fail', scores_all_six(n1))
gate('  arm: the earlier pass condition, which scored one residual by its sign, is caught',
     not scores_all_six(n1.replace('all six reported quantities reproduce, of which five are independent, and all four residuals are nonzero',
                                   'all three exact quantities reproduce, and the residual is nonzero', 1)))
# the feasibility table's per-step counts, against the logs those steps wrote
feas = doc[doc.index('## FEASIBILITY'):]
for step in ('L0', 'L1', 'H'):
    n = len([ln for ln in logs[step].split('\n') if ln.startswith('PASS')])
    f = len([ln for ln in logs[step].split('\n') if ln.startswith('FAIL')])
    row = [ln for ln in feas.split('\n') if ln.startswith(f'| {step} |')]
    gate(f'    step {step}: the doc reports {n} checks, and the log has {n} passing and {f} failing',
         len(row) == 1 and f'{n} checks' in row[0] and f == 0)
gate('  arm: a miscounted step is caught', '999 checks' not in doc)

own = ctl[ctl.index('=== s1_control.py ==='):]  # the log also carries step 3's own checks, since it re-runs them
n_own = len([ln for ln in own.split('\n') if ln.startswith('PASS')])
stated = re.search(r'`s1_control\.py` \| (\d+) checks', doc)
note = re.search(r'\*\*Correction, (\d{4}-\d{2}-\d{2}), after filing\.\*\*[^\n]*?It now runs (\d+) checks', doc)
# The filed row stays as filed. The document is true when the row matches the script, or when a dated note records the
# divergence and that note's count matches it: a correction, never a silent rewrite (ERP section 8).
FILED_CONTROL_CHECKS = '7'  # what the merged record carries in the N1 row; it is never rewritten, only annotated
ok = (bool(stated) and stated.group(1) == FILED_CONTROL_CHECKS
      and (FILED_CONTROL_CHECKS == str(n_own) or (bool(note) and note.group(2) == str(n_own))))
print(f'    script runs {n_own}; filed row says {stated.group(1) if stated else "none"}'
      + (f'; dated correction {note.group(1)} records {note.group(2)}' if note else '; no correction note'))
gate("  the control script's check count is recorded truthfully", ok)
gate('  arm: a correction note carrying the wrong count is caught', not (note and note.group(2) == str(n_own + 1)))
gate('  arm: rewriting the filed row is caught, even with the correction intact',
     re.search(r'`s1_control\.py` \| (\d+) checks', doc.replace(f'| {FILED_CONTROL_CHECKS} checks, on top of', '| 8 checks, on top of', 1)).group(1)
     != FILED_CONTROL_CHECKS)
gate('  arm: a correction note that records the filed count rather than the live one is caught',
     not (note and str(n_own) == FILED_CONTROL_CHECKS))

inv = open(os.path.join(M8S, 'DISCLOSURE_INVENTORY.md'), encoding='utf-8').read()
rec = len([ln for ln in inv.split('\n') if ln.startswith('| `')])
gate(f'  the disclosure names the inventory size the inventory actually has ({rec})', f'{rec} records' in doc)
gate('  arm: a wrong inventory size is caught', f'{rec + 1} records' not in doc)

print('(B) the worklist, which the rooms see')
LEAK = ['1188/5', '8800/43', '75/308', '200/903', '9/35', '1/924', '924', '463', '288', '400', '225',
        'seven lines', 'ten orbits', 'six points', 'v₂ + v₋₂', 'v₃ + v₋₃', '23/10', '6/5', '12/25',
        '(9, 0, 0)', '(0, 0, 10)', 'hexagon', 'octahedron', 'pyramid', 'prism', 'coherent', 'zonal',
        'minimum of `r̂₆` is', 'maximum is']
probe = work.replace('924\u00b7r\u0302\u2086', 'the normalization')  # the agreed item 11 asks for it by name
hits = [s for s in LEAK if s in probe]
gate(f'  no frozen value, count, name or representative leaks into the worklist (hits: {hits})', not hits)
gate('  arm: a planted leak is caught', any(s in probe + ' the ten orbits' for s in LEAK))
gate('  arm: the allowance does not admit a value over that normalization',
     any(t in probe + ' the maximum is 463/924' for t in LEAK))
gate('  the worklist gives D7 as an input, with both weights', '28/39' in work and '21/52' in work and 'not asked to derive' in work)
gate('  the worklist carries the zeros rule', 'reported as missing, never as zero' in work)
gate('  the worklist asks for arguments and for what would fail if omitted',
     'what would fail if that step were omitted' in work and work.count('argument') >= 6)
gate('  the worklist defines the quartic, the tangent space and the transverse space by formula',
     all(s in work for s in ('ρ₆(u)_Q', 'r̂₆(u) = ', 'T_u = ', 'O_u = ', 'N_u = ')))
def numerals(t):
    """Every numeral in the worklist that is not an item number, a section reference or a given dimension."""
    t = prose(t)
    t = re.sub(r'\*\*\d+\.\*\*', ' ', t)
    t = re.sub(r'#+ \d(\.\d)? ', ' ', t)
    t = re.sub(r'items? \d+( and \d+)?', ' ', t)
    t = re.sub(r'\u00a7 ?[\d.]+', ' ', t)
    t = re.sub(r'(complex )?dimension \d', ' ', t)
    t = re.sub(r'\d+(,| and)(?= \d)', ' ', t)
    return Counter(m.group(0).lower() for m in re.finditer(r'\b(ten|nine|eight|seven|six|five|four|three|two|one|\d+)\b', t, re.I))


# Pinned after reading every one of them: the ambient dimensions, the two sectors, and instructions that quantify
# nothing about the answer. Any numeral added later changes this counter, and the gate fails until it is read again.
PINNED = {'0': 2, '1': 2, '13': 1, '2': 2, '3': 1, '6': 1, '7': 1, '9': 1, 'four': 1, 'one': 14, 'seven': 1, 'two': 6}
cur = numerals(work)
gate(f'  every numeral in the worklist is a reviewed one (new: {dict(cur - Counter(PINNED))}, gone: {dict(Counter(PINNED) - cur)})', cur == Counter(PINNED))
gate('  arm: a count of the answer is caught', numerals(work + ' There are ten of them.') != cur)
gate('  the worklist asks for the census without naming its size', 'At each orbit of item 4' in work)
gate('  no em-dash in the worklist, and a planted one is refused', no_em(work) and not no_em(work + '—'))
ex = re.findall(r'`u = ([^`]+)`', work)
gate(f'  the worklist example vectors are generic ({len(ex)} of them, none a census representative)',
     len(ex) >= 3 and not any(v in ''.join(ex) for v in ('v₂ + v₋₂', 'v₃ + v₋₃', '√(23/10)')))

print('(C) how the two files render on GitHub')
import subprocess
from html.parser import HTMLParser


class Tables(HTMLParser):
    """Column counts per row, per table, from the rendered GFM."""

    def __init__(self):
        super().__init__()
        self.tables, self.row = [], None

    def handle_starttag(self, tag, attrs):
        if tag == 'table':
            self.tables.append([])
        elif tag == 'tr':
            self.row = 0
        elif tag in ('td', 'th') and self.row is not None:
            self.row += 1

    def handle_endtag(self, tag):
        if tag == 'tr' and self.tables:
            self.tables[-1].append(self.row)
            self.row = None


def columns(md):
    r = subprocess.run(['pandoc', '-f', 'gfm', '-t', 'html'], input=md, capture_output=True, text=True)
    t = Tables()
    t.feed(r.stdout)
    return t.tables


def rendered_code(md):
    r = subprocess.run(['pandoc', '-f', 'gfm', '-t', 'html'], input=md, capture_output=True, text=True)
    return r.stdout, re.findall(r'<code>(.*?)</code>', r.stdout, re.S)


ENT = {'&amp;': '&', '&lt;': '<', '&gt;': '>', '&quot;': '"'}


def unent(t):
    for k, v in ENT.items():
        t = t.replace(k, v)
    return t


for label, md in (('task doc', doc), ('worklist', work)):
    html, codes = rendered_code(md)
    t = Tables()
    t.feed(html)
    ragged = [i for i, tab in enumerate(t.tables) if len(set(tab)) != 1]
    gate(f'  {label}: every table renders with a uniform column count ({len(t.tables)} tables, ragged: {ragged})', not ragged)
    want = [c.replace('\\|', '|') for ln in md.split('\n') if ln.startswith('|') for c in re.findall(r'`([^`]+)`', ln)]
    got = {unent(c) for c in codes}
    lost = [c for c in want if c not in got]
    gate(f'  {label}: every code span in a table survives the render ({len(want)} spans, lost: {lost[:2]})', not lost)

bad_md = '| a | b |\n|---|---|\n| `x\u00b2` | `(1 + |z|\u00b2)` |\n'
_, bad_codes = rendered_code(bad_md)
gate('  arm: an unescaped pipe inside a table code span is caught', '(1 + |z|\u00b2)' not in {unent(c) for c in bad_codes})
def raw_pipes(md):
    """Table rows whose code spans, which alternate with the text between backticks, still carry a bare pipe."""
    out = []
    for ln in md.split('\n'):
        if ln.startswith('|'):
            parts = ln.split('`')
            if any(re.search(r'(?<!\\)\|', parts[i]) for i in range(1, len(parts), 2)):
                out.append(ln)
    return out


gate('  no unescaped pipe survives inside a table code span', not raw_pipes(doc) and not raw_pipes(work))
gate('  arm: the pipe check sees a bare pipe in a code span', bool(raw_pipes('| a | b |\n| `(1 + |z|)` | c |')))

print('(D) hashes')
for f in ('prereg/s1_task_details.md', 'prereg/worklist.md', 'prereg/s1_prereg_build.py', 's1_L0.py', 's1_L0_log.txt',
          's1_L1.py', 's1_L1_log.txt', 's1_H.py', 's1_H_log.txt', 's1_control.py', 's1_control_log.txt',
          'prereg/s1_prereg_arms.py', 'DISCLOSURE_INVENTORY.md',
          'S1_STEP1_L0.md', 'S1_STEP2_L1.md', 'S1_STEP3_H.md'):
    p = os.path.join(M8S, f)
    if os.path.exists(p):
        print(f'  {hashlib.sha256(open(p, "rb").read()).hexdigest()}  {f}')
    else:
        gate(f'  pinned file present: {f}', False)
cur = sorted(norm(n) for n in names)
if os.path.exists(INVENTORY):
    want = sorted(l.rstrip('\n') for l in open(INVENTORY, encoding='utf-8') if l.strip())
    # a multiset, so that a substitution inside a group, one gate dropped and another run twice, is a deficit and a
    # surplus rather than a silent swap
    gone = Counter(want) - Counter(cur)
    added = Counter(cur) - Counter(want)
    print(f'    {len(cur)} gates ran, {len(set(cur))} distinct after normalization')
    print(f'    gate_inventory.txt SHA-256 {hashlib.sha256(open(INVENTORY, "rb").read()).hexdigest()}')
    gate('  gate inventory matches the pin'
         + (f'; REMOVED: {dict(list(gone.items())[:3])}' if gone else '')
         + (f'; ADDED: {dict(list(added.items())[:3])}' if added else ''), not gone and not added)
elif os.environ.get('PIN') == '1':
    open(INVENTORY, 'w', encoding='utf-8').write('\n'.join(cur) + '\n')
    print(f'  gate inventory written with {len(cur)} gates, because PIN=1 was set')
else:
    # fail closed: without this, deleting the one file the guard compares against would silently re-pin and go green,
    # which is the same fail-open shape as the missing-file traceback (xrodz, #571)
    gate('  gate inventory is present; re-pin deliberately with PIN=1', False)
print(f'{len(fails)} FAILED: {fails}' if fails else 'ALL PASS')
sys.exit(1 if fails else 0)
