#!/usr/bin/env python3
"""Mutation testing of the audited script's PASS-style assertions.

Each mutation corrupts the OBJECT being checked (never the check itself) and
we record whether the run dies, and how the reported verdicts move.  An
assertion that survives its mutation proves nothing."""

import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(HERE, 'mut', 'base.py')

ANCHOR_SORT = "    chars.sort(key=lambda v: (F.to_fraction(v[0]), [x for x in v[1:]]))"
ANCHOR_VERIFY = "        chk = verify_table(gd)"
ANCHOR_STORE = "        G[nm] = gd\n        SD[nm] = GroupSD(gd)"

MUTATIONS = [
    ('M1_char_value_negated',
     'negate one 2I character value -> row/column orthogonality must break',
     ANCHOR_SORT,
     ANCHOR_SORT + "\n    if name == '2I':\n        chars[3][2] = F.neg(chars[3][2])"),

    ('M2_dimension_falsified',
     'set one 2O irreducible degree from 3 to 4 -> sum of squares must break',
     ANCHOR_SORT,
     ANCHOR_SORT + "\n    if name == '2O':\n        chars[5][0] = F.rat(4)"),

    ('M3_chi_def_shifted',
     'add 1 to one 2T defining-character value -> tensor positivity must break',
     ANCHOR_VERIFY,
     ANCHOR_VERIFY + "\n        if nm == '2T':\n            gd.chi_def[1] = gd.F.add(gd.chi_def[1], gd.F.ONE)"),

    ('M4_class_size_bumped',
     'BD_4 class size off by one after the table checks -> U12 must go red',
     ANCHOR_STORE,
     ANCHOR_STORE + "\n        if nm == 'BD_4':\n            gd.sizes[1] += 1"),

    ('M5_S_value_poisoned',
     'poison one S value in the U3 pool -> the residual test must go red',
     "    # pooled exact fit on a 3-row subset, then verify on the rest",
     "    data[-1][5] += Fraction(1, 7)\n"
     "    # pooled exact fit on a 3-row subset, then verify on the rest"),

    ('M6_chi_def_falsified_rational',
     'set one 2T defining-character value from 1 to 0 (stays rational, so no '
     'downstream rationality guard fires) -> does U2 notice?',
     "        G[nm] = gd\n        SD[nm] = GroupSD(gd)",
     "        G[nm] = gd\n        if nm == '2T':\n"
     "            gd.chi_def[1] = gd.F.rat(0)\n        SD[nm] = GroupSD(gd)"),

    ('M8_chi_def_permuted_2I',
     'swap two 2I defining-character values after all table checks -> which '
     'guard catches it?',
     "        G[nm] = gd\n        SD[nm] = GroupSD(gd)",
     "        G[nm] = gd\n        if nm == '2I':\n"
     "            gd.chi_def[2], gd.chi_def[5] = gd.chi_def[5], gd.chi_def[2]\n"
     "        SD[nm] = GroupSD(gd)"),

    ('M7_wrong_dynkin_diagram',
     'replace E8 by D8 (still positive definite) -> does anything assert?',
     "E8_EDGES = [(0, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (1, 3)]",
     "E8_EDGES = [(0, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (1, 2)]"),
]


def run(tag, src):
    d = os.path.join(HERE, 'mut', tag)
    os.makedirs(d, exist_ok=True)
    p = os.path.join(d, 'run.py')
    with open(p, 'w') as fh:
        fh.write(src)
    r = subprocess.run([sys.executable, 'run.py'], cwd=d, capture_output=True,
                       text=True, timeout=600)
    out = {'returncode': r.returncode,
           'died': r.returncode != 0,
           'last_error': r.stderr.strip().splitlines()[-1] if r.stderr.strip() else None,
           'trace_tail': r.stderr.strip().splitlines()[-6:] if r.stderr.strip() else []}
    j = os.path.join(d, 'm8_1_1_defect.json')
    if os.path.exists(j):
        with open(j) as fh:
            out['json'] = json.load(fh)
    return out


def main():
    base_src = open(BASE).read()
    base = run('BASE', base_src)
    assert not base['died'], 'baseline run failed'
    bj = base['json']
    results = []
    for tag, desc, old, new in MUTATIONS:
        assert old in base_src, f'{tag}: anchor not found'
        src = base_src.replace(old, new, 1)
        r = run(tag, src)
        rec = {'mutation': tag, 'what': desc, 'script_died': r['died'],
               'error': r['last_error'], 'trace_tail': r['trace_tail']}
        if 'json' in r:
            j = r['json']
            rec['U2_exact_identity_holds'] = j['U2']['exact_identity_holds']
            rec['U2_largest_numeric_deviation'] = j['U2']['largest_deviation_decimal']
            rec['U12_all_agree'] = j['U12']['all_agree']
            rec['U3_nonzero_residuals'] = j['U3']['pooled_nonzero_residuals']
            rec['U3_triple'] = j['U3']['universal_triple']
            rec['U11_roots'] = j['U11']['norm_minus2_vectors']
            rec['U11_det'] = j['U11']['det_positive_cartan']
            rec['U11_mod4_counts'] = j['U11']['mod4_value_counts']
            rec['U5_all_match'] = all(v['equals_(num_classes-1)/8'] for v in j['U5'].values())
            rec['2I_D_defining'] = [row['D'] for row in j['U1']['2I']['rows']
                                    if row['is_defining']]
            rec['detected'] = (r['died']
                               or not j['U12']['all_agree']
                               or j['U3']['pooled_nonzero_residuals'] > 0
                               or not j['U2']['exact_identity_holds']
                               or not rec['U5_all_match'])
        else:
            rec['detected'] = r['died']
        results.append(rec)
        print(f"{tag:26s} died={r['died']!s:6s} detected={rec['detected']}"
              f"  {r['last_error'] or ''}")
    base_rec = {'mutation': 'BASE', 'script_died': False,
                'U2_exact_identity_holds': bj['U2']['exact_identity_holds'],
                'U2_largest_numeric_deviation': bj['U2']['largest_deviation_decimal'],
                'U12_all_agree': bj['U12']['all_agree'],
                'U3_nonzero_residuals': bj['U3']['pooled_nonzero_residuals'],
                'U3_triple': bj['U3']['universal_triple'],
                'U11_roots': bj['U11']['norm_minus2_vectors'],
                'U11_det': bj['U11']['det_positive_cartan'],
                'U11_mod4_counts': bj['U11']['mod4_value_counts'],
                'U5_all_match': all(v['equals_(num_classes-1)/8'] for v in bj['U5'].values()),
                '2I_D_defining': [row['D'] for row in bj['U1']['2I']['rows']
                                  if row['is_defining']]}
    with open(os.path.join(HERE, 'mutation_results.json'), 'w') as fh:
        json.dump({'baseline': base_rec, 'mutations': results}, fh, indent=1)
    print('\nwrote mutation_results.json')


if __name__ == '__main__':
    main()
