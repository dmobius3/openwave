"""s2b_06_checker_coverage.py -- numbers about check_s3_maximum.py and its log,
computed rather than counted by eye.  The file is read as TEXT; it is never executed.
"""
import re, numpy as np
from pathlib import Path
HERE = Path(__file__).parent
src = (HERE/"check_s3_maximum.py").read_text().splitlines()
log = (HERE/"check_s3_maximum_log.txt").read_text().splitlines()

# gate() call sites, and which print-section they fall under
secs, cur = {}, None
order = []
for i, line in enumerate(src, 1):
    m = re.match(r"\s*print\('\((L\d)\)", line)
    if m:
        cur = m.group(1); order.append(cur); secs.setdefault(cur, [])
    if re.match(r"\s*gate\(", line):
        secs.setdefault(cur, []).append(i)
print("gate() call sites in check_s3_maximum.py, by section:")
tot = 0
for k in order:
    print("   %s : %d gates  (lines %s)" % (k, len(secs[k]), secs[k]))
    tot += len(secs[k])
print("   TOTAL gates = %d   (the note claims '30 checks')" % tot)

pas = [l for l in log if l.strip().startswith("PASS")]
fail = [l for l in log if l.strip().startswith("FAIL")]
print("\nlog: %d PASS lines in total, %d FAIL lines" % (len(pas), len(fail)))
i0 = next(i for i,l in enumerate(log) if l.startswith("(L1)"))
tail = log[i0:]
print("   PASS lines from (L1) onward (i.e. attributable to this script): %d"
      % len([l for l in tail if l.strip().startswith("PASS")]))
print("   sections printed by this script:", order)
heads = [l.strip() for l in log if re.match(r"^\([A-L]\d?\)", l.strip())]
print("   section headers present in the log:", heads)
print("   -> headers (A)-(F) are NOT printed by check_s3_maximum.py; the log is a concatenation.")

# the sweep
th = np.linspace(0, 2*np.pi, 20001)
spacing = th[1]-th[0]
axes = np.array([0, 2*np.pi/3, 4*np.pi/3])
dist = np.min(np.abs(((th[:,None]-axes[None,:])+np.pi) % (2*np.pi) - np.pi), axis=1)
print("\nthe numerical sweep in the checker:")
print("   20001 points, grid spacing = %.6e rad" % spacing)
print("   points EXCLUDED from the strictness gate (dist <= 0.01 rad): %d of %d"
      % (int((dist <= 0.01).sum()), len(th)))
print("   that is %.1f grid points per axis on each side" % (0.01/spacing))
print("   so within those bands the checker asserts nothing except at the axis itself.")

# does the checker's e(t) parametrisation cover the unit traceless diagonal circle?
e = np.sqrt(2/3)*np.stack([np.cos(th), np.cos(th-2*np.pi/3), np.cos(th+2*np.pi/3)], axis=1)
print("   parametrisation e(t): max |sum e| = %.2e, max ||e||-1 = %.2e  -> full unit traceless circle"
      % (np.abs(e.sum(axis=1)).max(), np.abs(np.linalg.norm(e,axis=1)-1).max()))
print("   e(0)*sqrt6 =", np.round(e[0]*np.sqrt(6), 12))

# dead line
print("\ndead code: line 61 is", src[60].strip()[:78])
print("            line 62 overwrites Bm:", src[61].strip()[:78])
print("   '.subs({e1 - e2: 4*b})' cannot fire (e1-e2 is not an atom), so line 61 is inert;")
print("   harmless, because line 62 rebuilds Bm correctly.")

# imports that are not supplied
imps = [l.strip() for l in src if l.strip().startswith("from ")]
print("\nunsupplied imports the checker depends on:", imps)
names = re.findall(r"import ([^\n]+)", " ".join(imps))
print("   names taken on trust:", names)
