"""Compares a dry-run RETURN.md against RESULTS.md, value by value, in exact rational arithmetic.

Items compared: 1 (multiplicities), 3 (lambda_2/g, plus R5 against 1 + w_6 * 9/35), 4 (every per-level
||P_n xi||^2/g^2 and ||P_n N||^2), 7 (lambda_4/g^2) and 8 (the ratios). A planted mutation of one value runs
after the real comparison and must produce exactly one mismatch, so a clean result can fail.
"""
import re
import sys
from fractions import Fraction as F

RAY = {'R1': 'coherent v3', 'R2': 'zonal v0', 'R3': 'octahedron', 'R4': 'hexagon'}
SEC = {'3-dim': "3′", '4-dim': '4'}
LEVELS = {"3′": [10, 14, 16, 18], '4': [8, 12, 14, 16, 18]}
den = lambda n: n * (n + 2) - 48
res = open('RESULTS.md', encoding='utf-8').read()


def rtables(sec):
    body = res.split(f'## {sec}\n', 1)[1].split('\n## ', 1)[0]
    out, cur = [], []
    for line in body.split('\n') + ['']:
        if line.startswith('|'):
            cur.append([c.strip() for c in line.strip('|').split('|')])
        elif cur:
            out.append(cur[2:])
            cur = []
    return out


MINE = {}
for s in ("3′", '4'):
    main, lv = rtables(f'Sector {s}')
    for (ray, q, l4, _, x2), row in zip(main, lv):
        MINE[(s, ray)] = dict(Q=F(q), lam4=F(l4), xi={n: F(v) for n, v in zip(LEVELS[s], row[1:])})
RATIOS = {r[0]: F(r[1]) for r in rtables('Sector comparison')[0]}


def section(text, head):
    return text.split(head, 1)[1].split('\n## ', 1)[0]


def compare(ret):
    ret = ret.replace('−', '-')
    ok, bad = 0, []

    def cmp(label, a, b):
        nonlocal ok
        if a == b:
            ok += 1
        else:
            bad.append((label, a, b))
    got = {int(r[0]): (int(r[1]), int(r[2])) for r in re.findall(r'\|\s*(\d+)\s*\|\s*(\d)\s*\|\s*(\d)\s*\|', section(ret, '## Item 1'))}
    for n in range(0, 19, 2):
        cmp(f'item1 n={n}', got.get(n), (int(n in (6, 10, 14, 16, 18)), 2 if n == 18 else int(n in (6, 8, 12, 14, 16))))
    for R, a, b in re.findall(r'\|\s*(R\d)[^|]*\|\s*([\d/]+)\s*\|\s*([\d/]+)\s*\|', section(ret, '## Item 3')):
        if R in RAY:
            cmp(f'item3 {R} 3′', F(a), MINE[("3′", RAY[R])]['Q'])
            cmp(f'item3 {R} 4', F(b), MINE[('4', RAY[R])]['Q'])
        else:
            cmp('item3 R5', (F(a), F(b)), (1 + F(28, 39) * F(9, 35), 1 + F(21, 52) * F(9, 35)))
    for secname, body in re.findall(r'### (\d-dim) sector\n(.*?)(?=\n### |\Z)', section(ret, '## Item 4'), flags=re.S):
        s = SEC[secname]
        for R, tbl in re.findall(r'\*\*(R\d):\*\*\n(.*?)(?=\n\*\*R\d|\Z)', body, flags=re.S):
            if R not in RAY:
                continue
            for n, cells in re.findall(r'^\|\s*(\d+)\s*\|(.*)\|\s*$', tbl, flags=re.M):
                cells = [c.strip().strip('*').strip() for c in cells.split('|')]
                xi = MINE[(s, RAY[R])]['xi'][int(n)]
                cmp(f'item4 {s} {R} n={n} xi', F(cells[-1]), xi)
                cmp(f'item4 {s} {R} n={n} N', F(cells[0]), den(int(n)) ** 2 * xi)
    for R, a, b in re.findall(r'\|\s*(R\d)\s*\|\s*(-?[\d/]+)\s*\|\s*(-?[\d/]+)\s*\|', section(ret, '## Item 7')):
        cmp(f'item7 {R} 3′', F(a), MINE[("3′", RAY[R])]['lam4'])
        cmp(f'item7 {R} 4', F(b), MINE[('4', RAY[R])]['lam4'])
    for R, q in re.findall(r'\|\s*(R\d)\s*\|\s*([\d/]+)\s*≈', section(ret, '## Item 8')):
        cmp(f'item8 {R}', F(q), RATIOS[RAY[R]])
    return ok, bad


if __name__ == '__main__':
    text = open(sys.argv[1] if len(sys.argv) > 1 else 'dryrun/RETURN.md', encoding='utf-8').read()
    ok, bad = compare(text)
    print(f'{ok} values match exactly; {len(bad)} mismatches')
    for b in bad:
        print('  MISMATCH', b)
    target = '1871763250/229518083223'
    assert target in text.replace('−', '-'), 'mutation target missing'
    m_ok, m_bad = compare(text.replace('−', '-').replace(target, '1871763251/229518083223', 1))
    assert len(m_bad) == 1 and m_ok == ok - 1, (m_ok, m_bad)
    print('planted mutation: exactly one mismatch, as required')
