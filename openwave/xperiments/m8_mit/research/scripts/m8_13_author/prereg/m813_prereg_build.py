"""Gates the M8.13 pre-registration: the task doc against the proof and the author-side check, the stage-1 worklist
against leakage, and the stage-2 instructions and N2's variant against the landed proof.

Gate names are stable identifiers; measured detail is printed beside them. Every gate that can be armed is armed with
a planted violation. The inventory of gate names is pinned in gate_inventory.txt and compared as a raw multiset; a
missing inventory fails, and re-pinning is an explicit opt-in with PIN=1.
Usage: python3 m813_prereg_build.py
"""
import hashlib
import os
import re
import subprocess
import sys
from collections import Counter

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
H = lambda s: hashlib.sha256(s.encode('utf-8')).hexdigest()
HF = lambda f: hashlib.sha256(open(os.path.join(HERE, f), 'rb').read()).hexdigest()
cells = lambda row: re.split(r'(?<!\\)\|', row)
prose = lambda s: re.sub(r'`[^`]*`', ' ', s)
fails, names = [], []


def gate(name, ok):
    names.append(name)
    print(('PASS ' if ok else 'FAIL ') + name, flush=True)
    if not ok:
        fails.append(name)


# Every input is checked for presence before any is read, so a missing file fails here rather than raising
REQUIRED = ['m8_13_task_details.md', 'worklist.md', 'stage2a_grading.md', 'stage2b_grading.md', 'step3_text.md',
            'S0_S3_MAXIMUM.md', 'check_s3_maximum.py', 'check_s3_maximum_log.txt', 'm813_equality.py', 'm813_equality_log.txt',
            'm8_12_worklist.md']
missing = [f for f in REQUIRED if not os.path.exists(os.path.join(HERE, f))]
print(f'    required inputs missing: {missing}')
gate('  every required input is present', not missing)
if missing:
    print(f'{len(fails)} FAILED: {fails}')
    sys.exit(1)
rd = lambda f: open(os.path.join(HERE, f), encoding='utf-8').read()
doc, work, st2a, st2b, var = (rd('m8_13_task_details.md'), rd('worklist.md'), rd('stage2a_grading.md'), rd('stage2b_grading.md'),
                              rd('step3_text.md'))


def sec(t, a, b):
    i = t.find(a)
    j = t.find(b, i + 1) if i >= 0 else -1
    return t[i:j] if i >= 0 and j >= 0 else ''

proof, eqlog, wl12 = rd('S0_S3_MAXIMUM.md'), rd('m813_equality_log.txt'), rd('m8_12_worklist.md')


print('(A) the task doc')
no_em = lambda s: '—' not in s
for label, text in (('task doc', doc), ('worklist', work), ('2a instructions', st2a), ('2b instructions', st2b), ('N2 variant', var)):
    gate(f'  no em-dash in the {label}', no_em(text))
gate('  arm: a planted em-dash is refused', not no_em(doc + '—'))
words = len(doc.split()) + len(work.split())
print(f'    task doc and worklist together: {words} words')
gate('  budget: inside the Section 12.2 cap of 8000', words <= 8000)
HEADS = ['## TASK PLANNING', '## SETTING', '## THE FIREWALL', '## DISCLOSURE', '## CANDIDATE PRE-REGISTERED CLAIMS',
         '## FROZEN VALUES', "## FEASIBILITY, AND THE AUTHOR'S DERIVATION", '## TO BE FIXED AT GO',
         '## AMENDMENT (before go, 2026-09-22)', '## DEFINITION OF DONE']
gate('  every required section heading is present', all(h in doc for h in HEADS))
gate('  arm: a renamed heading is caught', not all(h in doc.replace('## THE FIREWALL', '## THE WALL') for h in HEADS))
first = lambda t: re.search(r'(?<![A-Za-z])(I|my|My)(?![A-Za-z])', prose(t))
gate('  no first-person voice outside the code spans', not first(doc))
gate('  arm: a planted first person is caught', bool(first(doc + ' I checked it.')))

block = doc[doc.index('## CANDIDATE PRE-REGISTERED CLAIMS'):doc.index('## FROZEN VALUES')]
rows = {cells(ln)[1].strip(): cells(ln) for ln in block.split('\n') if re.match(r'^\| (P|U|N|D)\d \|', ln)}
for cid in ('P1', 'P2', 'U1', 'U2', 'U3', 'N1', 'N2'):
    r = rows.get(cid)
    gate(f'  claim {cid} has a claim, a standing, a pass and a fail condition', bool(r) and len(r) == 7 and all(x.strip() for x in r[2:6]))
gate('  diagnostic D1 is recorded, with no pass or fail column', 'D1' in rows and len(rows['D1']) == 4)
gate('  arm: a claim row with an empty pass condition is caught',
     not all(x.strip() for x in cells('| U9 | a claim | new |  | fails |')[2:6]))

print('(B) the frozen values, against the landed proof and the author-side check')
gate('  the landed proof is the pinned one', H(proof) == '7c634a33fdb0ee05d5a934395a414f6346ba26e049b07f2619961792ab725498')
gate('  arm: an edited proof fails the pin', H(proof + ' ') != '7c634a33fdb0ee05d5a934395a414f6346ba26e049b07f2619961792ab725498')
# Each frozen cell is pinned as written and its load-bearing quantities traced to the proof or the check: a value that
# appears more than once in the doc would let a presence test pass while one copy moved.
ftab = doc[doc.index('## FROZEN VALUES'):doc.index('## FEASIBILITY')]
frozen = {cells(ln)[1].strip(): cells(ln)[2].strip() for ln in ftab.split('\n')
          if ln.startswith('| ') and not ln.startswith('| quantity') and not ln.startswith('| ---')}
FROZEN = {
    'the maximum': ('`463/924`', proof, ('463/924',)),
    'the nematic bound': ('`TrN̄² ≤ 171/2`, equivalently `λ_max ≤ 15/√6` over unit traceless `e`', proof, ('TrN̄² ≤ 171/2', '15/√6')),
    'the constraint ellipse': ('`(2/3)a² + 8b² = 1`, with `e₃ = 2a/3` and `e₁ − e₂ = 4b`', proof, ('(2/3)a² + 8b² = 1', 'e₃ = 2a/3', 'e₁ − e₂ = 4b')),
    'the equality points on it': ('`(√6/2, 0)`, `(−√6/4, √6/8)` and `(−√6/4, −√6/8)`', eqlog, ('union over the blocks',)),
    'which blocks reach each': ('the first by M₁ and M₂; the second by M₁ and M₃; the third by M₂ and M₃', eqlog, ('M3: reaches',)),
    'the same points as `√6·e`': ('the three permutations of `(2, −1, −1)`', eqlog, ('sqrt6 * e =',)),
    'the top eigenspace at each': ('two-dimensional; at `√6·e = (−1, −1, 2)` it is `span{v₃, v₋₃}`', eqlog, ('span{v3, v-3}',)),
    'the maximizing set': ('the orbit of `(v₃ + v₋₃)/√2` under rotations and phase', proof, ('(v₃ + v₋₃)/√2',)),
    'the minimum and its set': ('`1/924`, on coherent states', proof, ('1/924',)),
}
print(f'    frozen table rows: {len(frozen)}')
gate('  the frozen table holds exactly the pinned quantities', set(frozen) == set(FROZEN))
for q, (cell, src, traces) in FROZEN.items():
    gate(f'    frozen: {q} reads as pinned, and its quantities are in its source', frozen.get(q) == cell and all(t in src for t in traces))
gate('  arm: an altered frozen cell is caught', frozen.get('the maximum') != '`464/924`')
own = [ln for ln in eqlog.split('\n') if ln.startswith(('PASS', 'FAIL'))]
print(f'    author-side check: {sum(l.startswith("PASS") for l in own)} pass, {sum(l.startswith("FAIL") for l in own)} fail')
gate('  the author-side check is green, with its arms', eqlog.rstrip().endswith('ALL PASS') and not any(l.startswith('FAIL') for l in own)
     and 'arm: the comparison would catch a missing point' in eqlog)
rerun = subprocess.run([sys.executable, 'm813_equality.py'], cwd=HERE, capture_output=True, text=True,
                       env={**os.environ, 'PYTHONHASHSEED': '0'})
print(f'    rerun of the author-side check: exit {rerun.returncode}, output {"identical to" if rerun.stdout + rerun.stderr == eqlog else "DIFFERENT FROM"} the log')
gate('  the author-side check reruns, exits 0, and prints exactly its log', rerun.returncode == 0 and rerun.stdout + rerun.stderr == eqlog)
n_eq = sum(l.startswith('PASS') for l in own)
gate('  the disclosure states the check count the log carries', f'{n_eq} checks, 0 failures' in doc)

# F1-5: the stage-2 inputs are frozen inputs, and every hash the doc states must be the hash of a file it names. The
# amendment adds the five instrument sources; the two logs are outputs of runs that read the doc, so they stay outside it.
STAGE2_IN = ('S0_S3_MAXIMUM.md', 'check_s3_maximum.py', 'check_s3_maximum_log.txt')
SOURCES = ('m813_equality.py', 'm813_equality_log.txt', 'm813_prereg_build.py', 'm813_prereg_arms.py', 'gate_inventory.txt')
named = {HF(f): f for f in STAGE2_IN + SOURCES if os.path.exists(os.path.join(HERE, f))}
stated = set(re.findall(r'`([0-9a-f]{64})`', doc))
print(f'    the doc states {len(stated)} SHA-256 values; unmatched: {[h[:16] for h in stated - set(named)]}')
gate('  every SHA-256 the task doc states is the hash of a file it names', bool(stated) and stated <= set(named))
gate('  the doc pins the proof, its checker and the checker\'s log', all(HF(f) in stated for f in STAGE2_IN))
gate('  arm: a hash that matches no file is caught', ('0' * 64) not in named)
# The chain from the review at #583, which pinned all seven instruments: the two the amendment left alone must still be
# byte-identical to the review's pins, and the doc must state those same values.
REVIEW_583 = {'m813_equality.py': '3577872c112cefe13413eab3a706e834314e635129f7a914f5488d43f3cc090d',
              'm813_equality_log.txt': 'a30b8198d329d2a038ea3253fce2253b007954db0581fa4b71aeac5d12de0354'}
gate('  the two instruments #583 left unchanged are byte-identical to its pins, and the doc states those pins',
     all(os.path.exists(os.path.join(HERE, f)) and HF(f) == h and h in stated for f, h in REVIEW_583.items()))
gate('  the doc pins the five instrument sources', all(os.path.exists(os.path.join(HERE, f)) and HF(f) in stated for f in SOURCES))


def surd(t):
    t = t.replace('−', '-').replace('√', 'sqrt')
    t = re.sub(r'sqrt(\d+)', r'sqrt(\1)', t)
    return sp.nsimplify(sp.sympify(t))


frozen_pts = re.findall(r'`\((−?√6/\d|\d), (−?√6/\d|0)\)`', frozen.get('the equality points on it', ''))
doc_pts = {(surd(x), surd(y)) for x, y in frozen_pts}
log_union = re.search(r'union over the blocks: \[(.+)\]', eqlog).group(1)
log_pts = {(sp.nsimplify(sp.sympify(p)), sp.nsimplify(sp.sympify(q))) for p, q in re.findall(r'\(([^,()]+(?:\([^()]*\))?[^,()]*), ([^()]+?(?:\([^()]*\))?[^()]*)\)', log_union)}
print(f'    doc freezes {len(doc_pts)} equality points; the log found {len(log_pts)}')
gate('  the equality points the doc freezes are exactly the ones the check found', len(doc_pts) == 3 and doc_pts == log_pts)
gate('  arm: a moved point would not match', {(surd('√6/3'), sp.Integer(0))} != {p for p in log_pts if p[1] == 0})
blocks = {nm: re.search(nm + r': reaches 15/sqrt\(6\) exactly at \[(.+)\]', eqlog).group(1) for nm in ('M1', 'M2', 'M3')}
gate('  which blocks reach which point, as the doc freezes it, is what the check logged',
     'sqrt(6)/2, 0' in blocks['M1'] and 'sqrt(6)/2, 0' in blocks['M2'] and 'sqrt(6)/2' not in blocks['M3']
     and 'the first by M₁ and M₂; the second by M₁ and M₃; the third by M₂ and M₃' in doc)
perms = re.findall(r'sqrt6 \* e = \[([^\]]+)\]', eqlog)
gate('  every equality point is a permutation of (2, -1, -1), with a two-dimensional top eigenspace',
     len(perms) == 3 and all(sorted(int(x) for x in p.split(',')) == [-1, -1, 2] for p in perms)
     and eqlog.count('15/sqrt6 has multiplicity 2') == 3 and 'the top eigenspace is span{v3, v-3}' in eqlog)

print('(C) the claims that carry the design')
u2 = rows.get('U2', [''] * 7)[2]
parts = ('the block decomposition', "each block's bound", 'the equality set of `λ_max = 15/√6` on the ellipse',
         'the passage from `‖Q‖ = 15/√6`', 'that eigenspace as `span{')
gate('  U2 grades step 3 in its five parts', all(p in u2 for p in parts))
gate('  U2 credits the route with any auditor-supplied completion named, not the text alone',
     'with any completion the auditor supplies named explicitly' in u2)
gate('  arm: a U2 missing a part is caught', not all(p in u2.replace(parts[2], 'the result') for p in parts))
gate('  U3 records the audited argument only if U1 and U2 both pass', 'U1 and U2 both pass' in rows.get('U3', [''] * 7)[4])
n2 = rows.get('N2', [''] * 7)
gate('  N2 scores noticing, names ESTABLISHED, SUPPLIED as the only scale-correct grade, and records the other two as misapplied',
     "explicitly notes that the tracing sentence omits the third block's equality case" in n2[4]
     and 'The scale-correct grade is ESTABLISHED, SUPPLIED' in n2[4]
     and 'each is recorded as a scale misapplication' in n2[4]
     and 'a derivation the grader writes is a supplied part' in n2[4]
     and "accepts the complete equality set without addressing the third block's equality case" in n2[5])
gate("  N2 is scored on the third block's case, not on the label",
     'N2 is scored on that content, not on the label' in n2[4] and 'whatever label it gives' in n2[5])
gate('  the pre-registered expectation says both texts should earn the same label',
     'the two texts are expected to earn the same label, and a matching label is not a failure of the control' in doc)
n1 = rows.get('N1', [''] * 7)
gate('  N1 is a known-answer check, and says it does not test that the maximum method transfers',
     'where the answer is known' in n1[2] and 'reaches the maximum but not the minimum' in n1[5]
     and 'with a completeness argument' in n1[4] and 'no completeness argument' in n1[5])
gate('  the doc says where uniqueness lives, and why N2 exists',
     '**Where uniqueness lives.**' in doc and 'Each of its three points is reached by two of the three blocks at once' in doc)
gate('  the doc states the two-block property for blocks, not for every traced case',
     'no single block is load-bearing for the set' in doc and 'removing any single case' not in doc
     and 'so omitting it would lose that point' in doc)
gate('  the firewall names the author-side check among what stage 1 cannot see',
     "the author's check described below" in doc and 'It receives the worklist only' in doc)
fw, run = sec(doc, '## THE FIREWALL', '## DISCLOSURE'), sec(doc, '### Ownership and run format', '### Sources of record')
gate('  the run format splits stage 2 into 2a and 2b, each return committed before the next part is handed over',
     "Stage 2 has two parts, and each return is committed before the next part's files are handed over" in run
     and "at 2a the auditor receives the grading instructions for 2a and N2's variant" in run
     and "at 2b it receives the author's argument with its checker and log" in run)
gate('  the firewall hands over the variant alone at 2a, withholding the argument and its checker until 2b',
     "At 2a it additionally receives `stage2a_grading.md` and N2's variant, `step3_text.md`, and nothing else" in fw
     and "its checker names the third block's equality point among its exact checks" in fw
     and 'At 2b it additionally receives `stage2b_grading.md`' in fw)
gate('  the firewall drops the S0 memo, which no file in the repository resolves', bool(fw) and 'S0 memo' not in fw)
amd = sec(doc, '## AMENDMENT (before go, 2026-09-22)', '## DEFINITION OF DONE')
gate('  the amendment records each change and the instrument chain',
     all(s in amd for s in ('**N2 is sequenced, so it cannot be passed by comparing texts.**', '**N2 is scored on content, not on the label.**',
                            '**The firewall names only what the repository resolves.**', '**The instruments.**',
                            "supersede the review's for those three files only", "and supersede the review's for those two files only")))
gate('  the network posture is load-bearing: a stage-1 room that reaches the public argument cannot score U1',
     '**Network posture: offline, and this is load-bearing for U1.**' in doc and 'public in this repository since #581' in doc
     and 'If network access is enabled and that argument is reached at stage 1 from any source, U1 cannot be scored as an independent reproduction, and stage 1 must be rerun offline.' in doc
     and 'exclud' not in doc
     and 'recorded as located' in doc and "U1 still passes only on the auditor's own completeness argument" in doc)
gate("  the doc pre-registers that the real note's part 3c is expected to grade ESTABLISHED, SUPPLIED",
     '**The author expects part 3c to grade ESTABLISHED, SUPPLIED, not ESTABLISHED as written.**' in doc)
gate("  a pass that needed a supplied part carries it into G2's status", "G2's status records what was supplied" in doc)
gate('  the controls heading covers a positive and a negative control', '### Controls' in doc and '### Negative controls' not in doc)
gate('  the disclosure names the author-side check and says it is withheld',
     '`m813_equality.py`' in doc and 'withheld from the auditor at both stages' in doc)
gate('  arm: a disclosure that drops the withholding is caught', 'withheld from the auditor at both stages' not in doc.replace('withheld', 'shown'))

print('(D) the stage-1 worklist, which the room sees')
setup12 = wl12[wl12.index('### 1.1 The space'):wl12.index('### 1.4 Fixed spaces')].rstrip()
gate("  the setup is M8.12's, byte for byte, so the conventions are ones a room has already run on",
     H(wl12) == 'c07f9bc64d39ffc99317188e28a91378be18c0fc8902490bd44bed76bba4efdc' and setup12 in work)
LEAK = ['463', '924', '171', '15/√6', '√6', 'hexagon', 'NOON', 'coherent', 'nematic', 'TrN', 'a₀₀', 'magnetiz', 'singlet',
        'span{', 'v₃ + v₋₃', '(2, −1, −1)', 'permutation', 'octahedron', 'prism', 'pyramid', 'ellipse', 'Kawaguchi', 'Romero']
hits = [t for t in LEAK if t in work]
print(f'    leak hits: {hits}')
gate('  no value, shape, route or name that the room is meant to derive appears in the worklist', not hits)
gate('  arm: a planted leak is caught', any(t in work + ' the hexagon' for t in LEAK))


def numerals(t):
    t = prose(t)
    t = re.sub(r'\*\*\d+\.\*\*', ' ', t)
    t = re.sub(r'#+ \d(\.\d)? ', ' ', t)
    t = re.sub(r'items? \d+( to \d+)?( and \d+)?', ' ', t)
    return Counter(m.group(0).lower() for m in re.finditer(r'\b(ten|nine|eight|seven|six|five|four|three|two|one|\d+)\b', t, re.I))


# Pinned after reading each in context: the ambient dimension, the section numbers, 'degree 0', and 'one markdown file'.
# None counts anything in the answer. A numeral added later changes the counter and fails until it is read again.
PINNED = {'one': 3, 'seven': 1, '1': 1, '0': 1, '2': 1, '3': 1}
print(f'    worklist numerals: {dict(numerals(work))}')
gate('  every numeral in the worklist is a reviewed one', numerals(work) == Counter(PINNED))
gate('  arm: a count of the answer is caught', numerals(work + ' There are three of them.') != numerals(work))
gate('  the worklist asks for completeness and what fails if a case is omitted',
     'give an argument that the set you report is complete' in work and 'what would fail if it were omitted' in work)
gate('  the worklist asks for the minimum, which is N1, without prescribing the method',
     '**3.** Answer items 1 and 2 for the minimum' in work and 'by the same method' not in work)
gate('  the worklist carries the zeros rule', 'reported as missing, never as zero' in work)
gate('  the worklist says a search establishes nothing on its own', 'on its own it establishes neither a maximum' in work)

print('(E) stage 2, in its two parts, and the N2 variant')
scales = (sec(st2a, '## How to grade', '## What to grade'), sec(st2b, '## How to grade', '## What to grade'))
gate('  the stage-2 instructions grade step 3 in parts 3a to 3e', all(f'**3{c}.**' in sec(st2b, '## What to grade', 'Then reconcile') for c in 'abcde'))
gate('  the stage-2 instructions grade the argument, not its conclusion', all('Grade the argument, not its conclusion.' in s for s in scales))
gate('  the stage-2 scale applies GAP strictly, so a completable step is never failed', all('GAP is strict: a step you can complete is never a GAP.' in s for s in scales))
gate('  the stage-2 scale treats a written-out derivation as supplied, so a completion cannot pass as ESTABLISHED',
     all('a derivation you write out is a supplied part, even when the text contains its ingredients' in s for s in scales))
gate('  the stage-2 scale defines ESTABLISHED, SUPPLIED, so a supplied case is named rather than silent',
     all(all(f'- **{v}**:' in s for v in ('ESTABLISHED', 'ESTABLISHED, SUPPLIED', 'GAP', 'DEFECT')) and 'Name exactly what you supplied.' in s for s in scales))
gate('  the 2a and 2b instructions carry the same scale, word for word', bool(scales[0]) and scales[0] == scales[1])
intro_a, intro_b = st2a[:st2a.find('## How to grade')], st2b[:st2b.find('## How to grade')]
gate('  the 2a instructions hand over exactly one file, the variant, and say the return is committed first',
     'Besides these instructions, you now receive one file: `step3_text.md`' in intro_a
     and 'Your stage-1 return is committed and will not change.' in intro_a
     and 'Your return on it is committed before any further file is handed over.' in intro_a)
HINTS = ['S0_S3_MAXIMUM', 'check_s3_maximum', 'author', 'variant', 'N2', 'third block', 'M₃', '√6/4', 'second text',
         'another text', 'control', 'planted', 'altered', 'modified', 'omitted', 'omission', 'differs']
print(f'    2a hint hits: {[t for t in HINTS if t in st2a]}')
gate('  the 2a instructions name neither the argument nor its checker, and hint at no omission or control', not any(t in st2a for t in HINTS))
gate('  the 2b instructions say both earlier returns are committed, and hand over exactly the argument, its checker and its log',
     'Your stage-1 and 2a returns are committed and will not change.' in intro_b
     and 'Besides these instructions, you now receive three files:' in intro_b
     and all(f in intro_b for f in ('`S0_S3_MAXIMUM.md`', '`check_s3_maximum.py`', '`check_s3_maximum_log.txt`')) and 'step3_text' not in intro_b)
n2para = st2a[st2a.find('## What to grade'):] if '## What to grade' in st2a else ''
gate('  the N2 paragraph asks for a grade as an argument, and names no expected verdict',
     'as an argument' in n2para and not re.search(r'\b(GAP|ESTABLISHED|DEFECT)\b', n2para))
gate('  arm: an N2 paragraph that leaks the expected verdict is caught', bool(re.search(r'\b(GAP|ESTABLISHED|DEFECT)\b', n2para + ' It is a GAP.')))
step3 = proof[proof.index('**Step 3: the nematic bound'):proof.index('**Step 4: the equality case')].rstrip()
a3 = 'Tracing the equality cases (`b = 0` with `a > 0`, `a = −2b` with `b > 0`, and `a = −√6/4` in `M₃`) shows'
b3 = 'Tracing the equality cases (`b = 0` with `a > 0` and `a = −2b` with `b > 0`) shows'
gate("  the N2 variant is the landed step 3 with exactly the third block's equality case removed",
     step3.count(a3) == 1 and var.rstrip('\n') == step3.replace(a3, b3, 1))
gate('  arm: a variant that also drops the third block\'s bound is caught',
     var.rstrip('\n') != step3.replace(a3, b3, 1).replace('- **`M₃`:**', '- **`M`:**'))
gate('  the renamed variant is byte-identical to the file #583 merged', HF('step3_text.md') == '05600d85936f8db7dadf0f131ff6cef83baf320eaab032970b53172c77e3fddf')

print('(F) how the files render on GitHub')


def rendered(md):
    return subprocess.run(['pandoc', '-f', 'gfm', '-t', 'html'], input=md, capture_output=True, text=True).stdout


ENT = {'&amp;': '&', '&lt;': '<', '&gt;': '>', '&quot;': '"'}
unent = lambda t: [t := t.replace(k, v) for k, v in ENT.items()][-1]
for label, md in (('task doc', doc), ('worklist', work), ('2a instructions', st2a), ('2b instructions', st2b)):
    html = rendered(md)
    got = {unent(c) for c in re.findall(r'<code>(.*?)</code>', html, re.S)}
    want = [c.replace('\\|', '|') for ln in md.split('\n') if ln.startswith('|') for c in re.findall(r'`([^`]+)`', ln)]
    lost = [c for c in want if c not in got]
    print(f'    {label}: {len(want)} table code spans, lost {lost[:2]}')
    gate(f'  {label}: every code span in a table survives the render', not lost)


def raw_pipes(md):
    out = []
    for ln in md.split('\n'):
        if ln.startswith('|'):
            parts = ln.split('`')
            if any(re.search(r'(?<!\\)\|', parts[i]) for i in range(1, len(parts), 2)):
                out.append(ln)
    return out


gate('  no unescaped pipe inside a table code span', not raw_pipes(doc) and not raw_pipes(work) and not raw_pipes(st2a) and not raw_pipes(st2b))
gate('  arm: a bare pipe in a table code span is caught', bool(raw_pipes('| a | `(1 + |z|)` |')))

print('(G) hashes')
for f in ('m8_13_task_details.md', 'worklist.md', 'stage2a_grading.md', 'stage2b_grading.md', 'step3_text.md', 'S0_S3_MAXIMUM.md',
          'check_s3_maximum.py', 'check_s3_maximum_log.txt', 'm813_equality.py', 'm813_equality_log.txt',
          'm8_12_worklist.md', 'm813_prereg_build.py', 'm813_prereg_arms.py'):
    p = os.path.join(HERE, f)
    if os.path.exists(p):
        print(f'  {hashlib.sha256(open(p, "rb").read()).hexdigest()}  {f}')
    else:
        gate(f'  pinned file present: {f}', False)

INVENTORY = os.path.join(HERE, 'gate_inventory.txt')
# Raw names: measured detail is printed beside a gate, never inside its name, so any renaming is a change of gate
cur = sorted(n.strip() for n in names)
if os.path.exists(INVENTORY):
    want = sorted(l.rstrip('\n') for l in open(INVENTORY, encoding='utf-8') if l.strip())
    gone, added = Counter(want) - Counter(cur), Counter(cur) - Counter(want)
    print(f'    {len(cur)} gates ran; gate_inventory.txt SHA-256 {hashlib.sha256(open(INVENTORY, "rb").read()).hexdigest()}')
    gate('  gate inventory matches the pin' + (f'; REMOVED: {dict(gone)}' if gone else '') + (f'; ADDED: {dict(added)}' if added else ''),
         not gone and not added)
elif os.environ.get('PIN') == '1':
    open(INVENTORY, 'w', encoding='utf-8').write('\n'.join(cur) + '\n')
    print(f'  gate inventory written with {len(cur)} gates, because PIN=1 was set')
else:
    gate('  gate inventory is present; re-pin deliberately with PIN=1', False)
print(f'{len(fails)} FAILED: {fails}' if fails else 'ALL PASS')
sys.exit(1 if fails else 0)
