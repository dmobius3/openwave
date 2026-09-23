"""The mutation suite for the M8.13 pre-registration gate: each arm plants one defect and requires the gate it was built
to trip to fail, not merely some gate.

Each arm runs against a fresh copy of the verified-green parent, staged in a scratch directory and removed afterwards.
An arm whose target text is missing, whose planted defect trips a different gate, or whose defect breaks the build
instead of failing a gate, is itself reported as a failure.
Usage: python3 m813_prereg_arms.py
"""
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
FILES = ['m8_13_task_details.md', 'worklist.md', 'stage2_grading.md', 'n2_variant_step3.md', 'S0_S3_MAXIMUM.md',
         'check_s3_maximum.py', 'check_s3_maximum_log.txt', 'm813_equality.py', 'm813_equality_log.txt',
         'm8_12_worklist.md', 'm813_prereg_build.py', 'm813_prereg_arms.py', 'gate_inventory.txt']
DOC, WL, ST2, VAR, INV = 'm8_13_task_details.md', 'worklist.md', 'stage2_grading.md', 'n2_variant_step3.md', 'gate_inventory.txt'

# (file, find, replace, what the arm plants, the gate it must trip)
ARMS = [
    (DOC, '## TASK PLANNING', '## TASK PLANNING —', 'an em-dash in the task doc', 'no em-dash in the task doc'),
    (DOC, '## THE FIREWALL', '## THE WALL', 'a required heading renamed', 'every required section heading is present'),
    (DOC, '## SETTING\n', '## SETTING\n\nI checked this.\n', 'a first person in the prose', 'no first-person voice'),
    (DOC, '| every step and every part of step 3 ESTABLISHED, as written or with the supplied part named |',
     '|  |', "U2's pass condition emptied", 'claim U2 has a claim'),
    (DOC, '| the maximum | `463/924` |', '| the maximum | `464/924` |', 'the frozen maximum altered', 'frozen: the maximum reads as pinned'),
    (DOC, '| the equality points on it | `(√6/2, 0)`', '| the equality points on it | `(√6/3, 0)`', 'an equality point moved',
     'the equality points the doc freezes are exactly'),
    (DOC, 'the first by M₁ and M₂; the second by M₁ and M₃', 'the first by M₁ and M₃; the second by M₁ and M₂',
     'the block membership swapped', 'which blocks reach which point'),
    (DOC, 'the equality set of `λ_max = 15/√6` on the ellipse', 'the result', "U2's equality-set part dropped",
     'U2 grades step 3 in its five parts'),
    (DOC, 'U1 and U2 both pass', 'U1 or U2 passes', "U3's rule weakened", 'U3 records the audited argument only if'),
    (DOC, "accepts the complete equality set without addressing the third block's equality case", 'grades 3c GAP',
     "N2's fail condition reverted to a verdict", 'N2 scores noticing, names ESTABLISHED, SUPPLIED as the only scale-correct grade'),
    (DOC, 'but each is recorded as a scale misapplication', 'and both are fine', "N2's scale-misapplication record dropped",
     'N2 scores noticing, names ESTABLISHED, SUPPLIED as the only scale-correct grade'),
    (DOC, 'The scale-correct grade is ESTABLISHED, SUPPLIED', 'An acceptable grade is ESTABLISHED', "N2 blessing plain ESTABLISHED again",
     'N2 scores noticing, names ESTABLISHED, SUPPLIED as the only scale-correct grade'),
    (ST2, 'a derivation you write out is a supplied part, even when the text contains its ingredients',
     'a derivation you write out may count as the text', 'written-out derivations allowed to pass as ESTABLISHED',
     'the stage-2 scale treats a written-out derivation as supplied'),
    (ST2, 'You now receive four files:', 'You now receive three files:', 'the file count wrong again',
     'the stage-2 instructions count the files they hand over correctly'),
    (DOC, 'with any completion the auditor supplies named explicitly, establishes', 'establishes',
     "U2 crediting the text alone again", 'U2 credits the route with any auditor-supplied completion named'),
    (DOC, "G2's status records what was supplied", "G2's status reads audited argument", "the supplied part dropped from G2's status",
     "a pass that needed a supplied part carries it into G2's status"),
    (DOC, '**The author expects part 3c to grade ESTABLISHED, SUPPLIED, not ESTABLISHED as written.**', '**A note on part 3c.**',
     'the pre-registered expectation dropped', "the doc pre-registers that the real note's part 3c"),
    (DOC, '### Controls', '### Negative controls', 'the controls heading reverted', 'the controls heading covers a positive'),
    (ST2, 'GAP is strict: a step you can complete is never a GAP. ', '', 'strict GAP dropped from the scale',
     'the stage-2 scale applies GAP strictly'),
    (DOC, "U1 still passes only on the auditor's own completeness argument", 'U1 may pass on located material', 'the online fallback weakened',
     'the network posture is load-bearing'),
    (DOC, 'the known value and set, with a completeness argument | a different set, or no completeness argument',
     'the known value and set | a different set', "N1's completeness requirement dropped", 'N1 is a known-answer check'),
    (DOC, 'no single block is load-bearing for the set', 'removing any single case leaves the conclusion unchanged',
     'the false any-single-case sentence restored', 'the doc states the two-block property for blocks'),
    (DOC, '`34054633f9602e864c12d017b035e6517ebbc4d8b0ba6c043dc1ef4adc764d3a`',
     '`34054633f9602e86d7db7c9fd0c35a1fc5f1b15ab4ea13e9f7aa4fde25c9e16d`',
     "a hash in the doc that matches no file, the author's own earlier slip", 'every SHA-256 the task doc states is the hash of a file'),
    (DOC, '25 checks, 0 failures', '24 checks, 0 failures', 'a disclosure count the log does not carry',
     'the disclosure states the check count the log carries'),
    ('m813_equality.py', "print('ALL PASS' if not fails else f'{len(fails)} FAILED: {fails}')",
     "print('ALL PASS' if not fails else f'{len(fails)} FAILED: {fails}')\nprint('BROKEN CHECKER')",
     'the check changed while its green log stays stale', 'the author-side check reruns, exits 0, and prints exactly its log'),
    ('m813_equality_log.txt', 'union over the blocks:', 'union over the blocks, as logged:',
     'the log changed while the check stays intact', 'the author-side check reruns, exits 0, and prints exactly its log'),
    (ST2, '<DELETE FILE>', '', 'an ordinary room input deleted', 'every required input is present'),
    (WL, '**3.** Answer items 1 and 2 for the minimum of `r̂₆` over the unit sphere.',
     '**3.** Answer items 1 and 2 for the minimum of `r̂₆` over the unit sphere, by the same method.',
     'the method-transfer instruction restored', 'the worklist asks for the minimum, which is N1, without prescribing'),
    (ST2, '- **ESTABLISHED, SUPPLIED**:', '- **ESTABLISHED (with help)**:', 'the supplied label dropped from the scale',
     'the stage-2 scale defines ESTABLISHED, SUPPLIED'),
    (DOC, 'offline, and this is load-bearing for U1', 'offline where convenient', 'the posture made optional', 'the network posture is load-bearing'),
    (DOC, 'public in this repository since #581, so', 'not public, so', 'the public exposure of the argument denied', 'the network posture is load-bearing'),
    (DOC, 'and stage 1 must be rerun offline.', 'and the room may continue.',
     'the rerun requirement dropped', 'the network posture is load-bearing'),
    (DOC, 'and stage 1 must be rerun offline.', 'and stage 1 must be rerun offline or with that repository excluded.',
     'the repository-exclusion fallback restored', 'the network posture is load-bearing'),
    (DOC, 'that argument is reached at stage 1 from any source,', 'this repository is reached at stage 1,',
     'the trigger narrowed to one repository', 'the network posture is load-bearing'),
    (DOC, "The worklist's standing request for anything looked up is the discriminator.",
     "The worklist's standing request for anything looked up is the discriminator. A room may instead run with the author's fork excluded.",
     'an exclusion fallback added in another sentence', 'the network posture is load-bearing'),
    (DOC, 'withheld from the auditor at both stages', 'shown to the auditor at stage 2', 'the withholding dropped',
     'the disclosure names the author-side check and says it is withheld'),
    (DOC, '| `(2/3)a² + 8b² = 1`, with', '| `(2/3)a² + 8b² = 2`, with', 'the frozen ellipse altered',
     'frozen: the constraint ellipse reads as pinned'),
    (WL, '**1.** Find the maximum', '**1.** The maximum is attained at the hexagon. Find the maximum', 'a shape leaked into the worklist',
     'no value, shape, route or name'),
    (WL, '**2.** Find the complete set', '**2.** There are three cases. Find the complete set', 'a count leaked into the worklist',
     'every numeral in the worklist is a reviewed one'),
    (WL, 'in the Condon-Shortley convention', 'in the standard convention', 'the setup altered from M8.12',
     "the setup is M8.12's, byte for byte"),
    (WL, 'give an argument that the set you report is complete', 'say what the set is', 'the completeness request dropped',
     'the worklist asks for completeness'),
    (WL, 'reported as missing, never as zero', 'reported as zero', 'the zeros rule dropped', 'the worklist carries the zeros rule'),
    (ST2, '   - **3c.**', '   - **3x.**', 'part 3c dropped from stage 2', 'the stage-2 instructions grade step 3 in parts 3a to 3e'),
    (ST2, 'as an argument, and say why.', 'as an argument; it should be a GAP.', 'the expected N2 verdict leaked',
     'the N2 paragraph asks for a grade as an argument, and names no expected verdict'),
    (ST2, 'Grade the argument, not its conclusion. ', '', 'the argument-not-conclusion rule dropped',
     'the stage-2 instructions grade the argument, not its conclusion'),
    (VAR, '(`b = 0` with `a > 0` and `a = −2b` with `b > 0`)', '(`b = 0` with `a > 0`, `a = −2b` with `b > 0`, and `a = −√6/4` in `M₃`)',
     'the variant restored to the full tracing', 'the N2 variant is the landed step 3 with exactly'),
    (VAR, '- **`M₃`:**', '- **`M`:**', "the variant also altering the third block's bound", 'the N2 variant is the landed step 3 with exactly'),
    ('S0_S3_MAXIMUM.md', '**Step 4: the equality case.**', '**Step 4: the equality case.** ', 'the landed proof edited',
     'the landed proof is the pinned one'),
    ('m813_equality_log.txt', 'PASS every equality point is reached by exactly two blocks', 'FAIL every equality point is reached by exactly two blocks',
     "the author-side check failing", 'the author-side check is green'),
    (DOC, '| U2 | The route of `S0_S3_MAXIMUM.md` steps 1 to 4', '| U2 | The route of `(1 + |z|)` and `S0_S3_MAXIMUM.md` steps 1 to 4',
     'an unescaped pipe in a table code span', 'no unescaped pipe inside a table code span'),
    (INV, '<DELETE FILE>', '', 'the gate inventory deleted', 'gate inventory is present; re-pin deliberately'),
    (INV, 'the landed proof is the pinned one\n', '', 'a line cut from the pin, to match a deleted gate', 'gate inventory matches the pin'),
    ('m813_prereg_build.py', "gate('  every equality point is a permutation of (2, -1, -1),", "gate('  every equality point is a permutation of (1, 1, -2),",
     'a gate renamed inside its parentheses', 'gate inventory matches the pin'),
]


def build(root):
    return subprocess.run([sys.executable, 'm813_prereg_build.py'], cwd=root, capture_output=True, text=True,
                          env={**os.environ, 'PYTHONHASHSEED': '0'})


def stage(root):
    for f in FILES:
        shutil.copy(os.path.join(HERE, f), os.path.join(root, f))


root = tempfile.mkdtemp(prefix='m813_arms_')
caught = missed = 0
try:
    stage(root)
    r = build(root)
    print('parent: ' + ('GREEN' if r.returncode == 0 else 'RED, so no arm below proves anything'), flush=True)
    if r.returncode != 0:
        print(r.stdout[-1500:])
        sys.exit(1)
    for i, (f, a, b, what, want) in enumerate(ARMS, 1):
        stage(root)
        p = os.path.join(root, f)
        if a == '<DELETE FILE>':
            os.remove(p)
        else:
            s = open(p, encoding='utf-8').read()
            if a not in s:
                print(f'{i:3d} UNAPPLIED  {what}: the text to mutate was not found')
                missed += 1
                continue
            open(p, 'w', encoding='utf-8').write(s.replace(a, b, 1))
        out = build(root)
        hit = [ln for ln in out.stdout.split('\n') if ln.startswith('FAIL')]
        if not hit:
            missed += 1
            print(f'{i:3d} {"BROKE THE BUILD" if out.returncode else "NOT CAUGHT"} {what}')
        elif not any(want in ln for ln in hit):
            missed += 1
            print(f'{i:3d} WRONG GATE {what}: expected {want[:40]!r}, tripped {hit[0][5:].strip()[:40]!r}')
        else:
            caught += 1
            print(f'{i:3d} caught     {what}')
    stage(root)
    print('restored: ' + ('GREEN' if build(root).returncode == 0 else 'RED'), flush=True)
finally:
    shutil.rmtree(root, ignore_errors=True)
print(f'{caught} of {len(ARMS)} arms caught' + ('' if not missed else f', {missed} NOT CAUGHT'))
sys.exit(1 if missed else 0)
