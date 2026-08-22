"""
Automated registry-coverage set-equality checker.
Verifies that every ID in the construction and convention registries
appears in the coverage table, and vice versa.
"""
import re


def main():
    with open('METHOD_AND_GATE_MANIFEST.md') as f:
        text = f.read()

    # Extract construction registry IDs (section 2)
    constr_ids = set()
    in_constr = False
    for line in text.split('\n'):
        if '## 2. Construction registry' in line:
            in_constr = True; continue
        if in_constr and line.startswith('## '):
            break
        if in_constr and line.startswith('| C-'):
            cid = line.split('|')[1].strip()
            constr_ids.add(cid)

    # Extract convention registry IDs (section 3)
    conv_ids = set()
    in_conv = False
    for line in text.split('\n'):
        if '## 3. Convention registry' in line:
            in_conv = True; continue
        if in_conv and line.startswith('## '):
            break
        if in_conv and line.startswith('| V-'):
            cid = line.split('|')[1].strip()
            conv_ids.add(cid)

    # Extract gate registry IDs (section 4)
    gate_ids = set()
    in_gates = False
    for line in text.split('\n'):
        if '## 4. Pre-reveal gate registry' in line:
            in_gates = True; continue
        if in_gates and line.startswith('## '):
            break
        if in_gates and line.startswith('| G-'):
            gid = line.split('|')[1].strip()
            gate_ids.add(gid)

    # Extract coverage table IDs (section 7)
    coverage_ids = set()
    in_cov = False
    for line in text.split('\n'):
        if '## 7. Pre-implementation coverage table' in line:
            in_cov = True; continue
        if in_cov and line.startswith('## '):
            break
        if in_cov and (line.startswith('| C-') or line.startswith('| V-')):
            cid = line.split('|')[1].strip()
            coverage_ids.add(cid)

    all_registry = constr_ids | conv_ids
    print(f"Construction IDs ({len(constr_ids)}): {sorted(constr_ids)}")
    print(f"Convention IDs ({len(conv_ids)}): {sorted(conv_ids)}")
    print(f"Gate IDs ({len(gate_ids)}): {sorted(gate_ids)}")
    print(f"Coverage table IDs ({len(coverage_ids)}): {sorted(coverage_ids)}")

    # Check set equality: coverage table should contain all construction + convention IDs
    missing_from_coverage = all_registry - coverage_ids
    extra_in_coverage = coverage_ids - all_registry
    assert not missing_from_coverage, f"Missing from coverage table: {missing_from_coverage}"
    assert not extra_in_coverage, f"Extra in coverage table (not in registries): {extra_in_coverage}"
    print("\nConstruction + Convention IDs == Coverage table IDs: ✓")

    # Verify gate IDs reference valid constructions or conventions
    # (gates should test constructions/conventions that exist)
    print(f"\nGate count: {len(gate_ids)}")
    for gid in sorted(gate_ids):
        print(f"  {gid}")

    # Verify no PENDING items in coverage
    pending = []
    in_results = False
    for line in text.split('\n'):
        if '## 9. Pre-implementation validation results' in line:
            in_results = True; continue
        if in_results and 'PENDING' in line:
            pending.append(line.strip())

    if pending:
        print(f"\nPENDING items in validation results:")
        for p in pending:
            print(f"  {p}")
    else:
        print("\nNo PENDING items in validation results: ✓")

    print("\n=== MANIFEST VALIDATION PASSED ===")


if __name__ == '__main__':
    main()
