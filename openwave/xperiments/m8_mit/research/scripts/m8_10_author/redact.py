"""M8.10 author package: the privacy redaction of two pinned files, the landing tree, and its gates.

Fourteen of the sixteen files pinned in the #546 task doc land byte-identical. Two carry machine-local
paths and, in the transcript, the author's startup-environment and account records; they land as
deterministic redactions of the pinned originals. This script is the whole transformation. It fails
closed: it refuses an existing output directory, assembles and checks the tree in a temporary
directory, and moves the tree into place only after every gate passes.

Rules, in order:
  T1  transcript only: the one system:init line and every rate_limit_event line become marker lines that
      keep the event type and carry the SHA-256 of the exact original line (without its newline); the
      init marker also keeps the model name
  P1  the author's Claude Code session-data directory under the home directory -> <claude-project>
  P2  the author's home directory -> <home>

Private inputs, supplied at run time so that this file names none of them:
  M8_10_HOME            the author's home directory, the only private input the rules use
  M8_10_PRIVATE         colon-separated identifiers that must be absent afterwards, in any case
  M8_10_SENSITIVE_FILE  regular expressions, one per line, naming fields of the removed records; each must
                        match the pinned transcript, none may match afterwards, and the file's SHA-256 goes
                        in the manifest

Gates:
  G1  each of the sixteen inputs matches its SHA-256 pinned in the task doc
  G2  zero residuals in the landed UTF-8 text of all eighteen files (the sixteen, this script and the
      manifest): the four check_no_local.py patterns, copied verbatim; JSON-escaped and percent-encoded
      home paths, macOS temp paths and email addresses; the home directory; the sensitive patterns; the
      private identifiers. A home path the rules cannot reach, such as a JSON-escaped one, therefore
      refuses the write
  G3  the written transcript against the pinned original: equal line counts and a final newline; every
      line parses as JSON; each startup or rate-limit record is replaced by its marker, whose hash is the
      SHA-256 of that record's line; and every other line, with the path rules reversed, equals its
      original exactly
  G4  the written audit script against its original: with the path rules reversed, it equals the
      original exactly; and the fourteen byte-identical files, re-read from disk, match their pins
  G5  the landed audit script, run on the landed transcript, reproduces the containment audit (55 tool
      calls, CLEAN), and the landed comparison script reproduces 103 matches with no mismatch; both
      must exit 0
G3 and G4 run twice: before G5 executes any landed code, and again on the finished tree. The manifest's
table is then read back and checked against the pins and the landed files, and after the move every
landed file is re-hashed against the tree the gates checked.

Usage:
  M8_10_HOME=... M8_10_PRIVATE=... M8_10_SENSITIVE_FILE=... python3 redact.py <pinned_dir> <task_doc> <out_dir>
"""
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REDACT = ('audit_dryrun.py', 'dryrun/transcript.jsonl')
NO_LOCAL = [  # dev_docs/utils/check_no_local.py PATTERNS, verbatim
    r"/(?:Users|home)/[A-Za-z0-9._-]+/",
    r"-(?:Users|home)-[A-Za-z0-9._-]+-",
    r"[A-Za-z]:\\\\?Users\\\\?[A-Za-z0-9._-]+",
    r"/private/(?:tmp|var)/|/tmp/[A-Za-z0-9._-]+-\d+/",
]
# JSON-escaped and percent-encoded home paths, macOS temp paths, email addresses. G2 sweeps this file
# too, so each pattern is written in a form that does not match itself.
STRUCTURAL = [r"\\/Users\\/", r"(?i)%2FUsers%2[F]", r"/var/folder[s]/", r"[\w.+-]+@[\w-]+\.[A-Za-z]{2,}"]
PLACEHOLDERS = ('<home>', '<claude-project>')
MANIFEST_ROW = r'^\| `([^`]+)` \| ([A-Z-]+) \| `([0-9a-f]{64})` \| `([0-9a-f]{64})` \| [^|]+ \|$'
sha = lambda b: hashlib.sha256(b).hexdigest()


def fail(msg):
    sys.exit(f'REDACTION REFUSED: {msg}')


def pins_from(task_doc):
    rows = re.findall(r'^\| `([^`]+)` \| `([0-9a-f]{64})` \|$', Path(task_doc).read_text(encoding='utf-8'), flags=re.M)
    if len(rows) != 16 or len(dict(rows)) != 16:
        fail(f'expected 16 distinct pinned files in the task doc, found {len(rows)} rows')
    return dict(rows)


def project_dir_pattern(home):
    return re.escape(home) + r'/\.claude/projects/[A-Za-z0-9._-]+'


def forward(text, home, counts):  # P1, then P2
    text, n = re.subn(project_dir_pattern(home), '<claude-project>', text)
    counts['P1'] = counts.get('P1', 0) + n
    counts['P2'] = counts.get('P2', 0) + text.count(home)
    return text.replace(home, '<home>')


def reverse(text, home, proj):
    return text.replace('<claude-project>', proj).replace('<home>', home)


def is_removed_record(ev):
    return ev.get('type') == 'rate_limit_event' or (ev.get('type') == 'system' and ev.get('subtype') == 'init')


def marker_for(ev, line):  # T1
    marker = {'type': ev['type']}
    if ev['type'] == 'system':
        marker.update(subtype='init', model=ev.get('model'), redacted='startup environment record removed')
    else:
        marker['redacted'] = 'record removed'
    marker['sha256_removed'] = sha(line.encode('utf-8'))
    return marker


def redact_transcript(raw, home, counts):
    lines = raw.split('\n')
    if lines[-1] != '':
        fail('the transcript does not end with a newline')
    out = []
    for line in lines[:-1]:
        ev = json.loads(line)
        if is_removed_record(ev):
            counts['T1'] = counts.get('T1', 0) + 1
            out.append(json.dumps(marker_for(ev, line), ensure_ascii=False, separators=(',', ':')))
        else:
            out.append(forward(line, home, counts))
    return '\n'.join(out) + '\n'


def check_written(tmp, src, home, pins):  # G3 and G4: the files on disk against the pinned originals
    for rel, pin in pins.items():
        landed, orig = (tmp / rel).read_bytes(), (src / rel).read_bytes()
        if sha(orig) != pin:
            fail(f'{rel}: the original no longer matches its pin')
        if rel not in REDACT:
            if sha(landed) != pin:
                fail(f'{rel}: the written bytes do not match the pin')
            continue
        lt, ot = landed.decode('utf-8'), orig.decode('utf-8')
        if any(p in ot for p in PLACEHOLDERS):
            fail(f'{rel}: the original contains a placeholder, so reversal would be ambiguous')
        dirs = set(re.findall(project_dir_pattern(home), ot))
        if len(dirs) > 1:
            fail(f'{rel}: more than one session-data directory, so reversal would be ambiguous')
        proj = dirs.pop() if dirs else '<claude-project>'
        if rel.endswith('.py'):
            if reverse(lt, home, proj) != ot:
                fail(f'{rel}: the written file does not reverse to its original')
            continue
        ll, ol = lt.split('\n'), ot.split('\n')
        if len(ll) != len(ol) or ll[-1] != '' or ol[-1] != '':
            fail(f'{rel}: the line count or the final newline changed')
        for i, (a, b) in enumerate(zip(ll[:-1], ol[:-1]), 1):
            try:
                ea = json.loads(a)
            except ValueError:
                fail(f'{rel} line {i}: not JSON')
            eb = json.loads(b)
            if is_removed_record(eb):
                if ea != marker_for(eb, b):
                    fail(f'{rel} line {i}: a startup or rate-limit record is not replaced by its marker')
            elif reverse(a, home, proj) != b:
                fail(f'{rel} line {i}: the line differs from its original beyond the path rules')


def write_manifest(tmp, rows, sensitive_sha):
    L = ['# M8.10 author package: landing manifest', '',
         'Sixteen files were pinned by SHA-256 in [`../../tasks/m8_10_task_details.md`](../../tasks/m8_10_task_details.md). '
         'Fourteen land byte-identical. Two land as deterministic privacy redactions of the pinned originals, produced by '
         '[`redact.py`](redact.py), whose docstring states the rules and the gates. With `redact.py` and this manifest, the '
         'directory holds eighteen files.', '',
         '| file | status | pinned SHA-256 | landed SHA-256 | substitutions |', '| --- | --- | --- | --- | --- |']
    for rel, status, pin, landed, counts in rows:
        cnt = ', '.join(f'{k} {v}' for k, v in sorted(counts.items())) or '-'
        L.append(f'| `{rel}` | {status} | `{pin}` | `{landed}` | {cnt} |')
    L += ['', 'T1 replaces a whole line with a marker carrying the SHA-256 of the removed line; P1 and P2 are the path rules. '
          'G3 and G4 checked the files as written against the pinned originals, before any landed code ran and again on '
          'the finished tree: every startup and rate-limit record is replaced by the marker carrying its hash, and every '
          'other line, with the path rules reversed, equals its original exactly. On the landed files, `audit_dryrun.py` '
          'reproduces the containment audit (55 tool calls, CLEAN) and `compare_dryrun.py` reproduces 103 matches with no '
          'mismatch, both exiting 0.', '',
          'G2 found zero residuals in the landed UTF-8 text of all eighteen files, this manifest and `redact.py` included, '
          'so a home path the rules could not reach, such as a JSON-escaped one, would have refused the write. The '
          'sensitive-pattern list G2 applied was supplied at run time, so this repository names none of its entries; '
          f"each entry matched the pinned transcript, so none was dead. The list's SHA-256 is `{sensitive_sha}`. "
          '`check_no_local.py` does not read `.jsonl` files, so G2 is what covers the transcript.']
    (tmp / 'MANIFEST.md').write_text('\n'.join(L) + '\n', encoding='utf-8')


def main(pinned_dir, task_doc, out_dir):
    home = os.environ.get('M8_10_HOME', '')
    private = [s for s in os.environ.get('M8_10_PRIVATE', '').split(':') if s]
    sens_file = os.environ.get('M8_10_SENSITIVE_FILE', '')
    if not home.startswith('/') or not private or not sens_file or not Path(sens_file).is_file():
        fail('set M8_10_HOME, M8_10_PRIVATE and M8_10_SENSITIVE_FILE')
    sens_bytes = Path(sens_file).read_bytes()
    sensitive = [re.compile(p) for p in sens_bytes.decode('utf-8').split('\n') if p]
    if not sensitive:
        fail('the sensitive-pattern file is empty')
    out = Path(out_dir)
    if out.exists():
        fail(f'the output directory already exists: {out}')
    pins, src = pins_from(task_doc), Path(pinned_dir)
    tmp = Path(tempfile.mkdtemp(prefix='m8_10_landing_'))
    rows = []
    for rel, pin in pins.items():
        raw = (src / rel).read_bytes()
        if sha(raw) != pin:  # G1
            fail(f'{rel} does not match its pin')
        dest = tmp / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        if rel not in REDACT:
            dest.write_bytes(raw)
            rows.append([rel, 'BYTE-IDENTICAL', pin, None, {}])
            continue
        counts, text = {}, raw.decode('utf-8')
        if rel.endswith('.jsonl'):
            dead = [rx.pattern for rx in sensitive if not rx.search(text)]
            if dead:
                fail(f'{len(dead)} sensitive pattern(s) match nothing in the pinned transcript')
            text = redact_transcript(text, home, counts)
        else:
            text = forward(text, home, counts)
        dest.write_bytes(text.encode('utf-8'))
        rows.append([rel, 'PRIVACY-REDACTED', pin, None, counts])
    check_written(tmp, src, home, pins)  # G3, G4, before any landed code runs
    a = subprocess.run([sys.executable, 'audit_dryrun.py', 'dryrun/transcript.jsonl'], cwd=tmp, capture_output=True, text=True)
    c = subprocess.run([sys.executable, 'compare_dryrun.py', 'dryrun/RETURN.md'], cwd=tmp, capture_output=True, text=True)
    if a.returncode != 0 or '55 tool calls audited; 0 flags' not in a.stdout or 'CLEAN' not in a.stdout:  # G5
        fail(f'the containment audit on the landed files (exit {a.returncode}): {a.stdout.strip()[-300:]} {a.stderr.strip()[-300:]}')
    if c.returncode != 0 or '103 values match exactly; 0 mismatches' not in c.stdout:
        fail(f'the comparison on the landed files (exit {c.returncode}): {c.stdout.strip()[-300:]} {c.stderr.strip()[-300:]}')
    for junk in list(tmp.rglob('__pycache__')):
        shutil.rmtree(junk)
    shutil.copy2(__file__, tmp / 'redact.py')
    for row in rows:
        row[3] = sha((tmp / row[0]).read_bytes())
    write_manifest(tmp, rows, sha(sens_bytes))
    check_written(tmp, src, home, pins)  # G3, G4 again, on the finished tree
    files = sorted(p for p in tmp.rglob('*') if p.is_file())
    if len(files) != 18:
        fail(f'expected eighteen landed files, found {len(files)}')
    patterns = [re.compile(p) for p in NO_LOCAL + STRUCTURAL] + sensitive
    for p in files:  # G2
        try:
            text = p.read_bytes().decode('utf-8')
        except UnicodeDecodeError:
            fail(f'{p.relative_to(tmp)} is not UTF-8')
        low = text.lower()
        hits = [(rx.pattern, m.group(0)) for rx in patterns for m in rx.finditer(text)]
        hits += [(s, s) for s in [home] + private if s.lower() in low]
        if hits:
            fail(f'{p.relative_to(tmp)}: {len(hits)} residual(s), first {hits[:2]}')
    checked = {str(p.relative_to(tmp)): sha(p.read_bytes()) for p in files}
    table = re.findall(MANIFEST_ROW, (tmp / 'MANIFEST.md').read_text(encoding='utf-8'), flags=re.M)
    if sorted((r, p, l) for r, _, p, l in table) != sorted((r, pins[r], checked[r]) for r in pins) or \
            any((s == 'BYTE-IDENTICAL') != (r not in REDACT) for r, s, _, _ in table):
        fail('the manifest table does not match the pins and the landed files')
    if checked['redact.py'] != sha(Path(__file__).read_bytes()):
        fail('the landed redact.py is not this script')
    if out.exists():
        fail(f'the output directory appeared during the run: {out}')
    shutil.move(str(tmp), str(out))
    moved = {str(p.relative_to(out)): sha(p.read_bytes()) for p in out.rglob('*') if p.is_file()}
    if moved != checked:
        fail(f'the moved tree at {out} differs from the checked tree; delete it before rerunning')
    print(f'landing tree written: {out} (18 files)')
    for rel, status, pin, landed, counts in rows:
        if status != 'BYTE-IDENTICAL':
            print(f'  {rel}: {status}, landed {landed[:16]}..., {counts}')
    print('G1 to G5 passed; the manifest table and the moved tree match the checked tree')


if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
