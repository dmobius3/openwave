#!/usr/bin/env python3
"""Automated manifest validator for M8.8.

Reads METHOD_AND_GATE_MANIFEST.md and verifies:
1. Exact set equality between construction+convention registry IDs
   and the coverage table's Registry ID column.
2. Every checkable row (Checkable Pre-Impl? = YES) has Result = PASS.
3. Every gate registry ID appears in the coverage table.
4. No duplicates in any registry or the coverage table.
5. No row with a prose-only check (must name a validation artifact or
   explicitly state NOT checkable).
"""

import re
import sys

def parse_table(text, header_pattern):
    """Parse a markdown table that starts with a header matching the pattern.
    Returns list of dicts keyed by header names."""
    lines = text.split('\n')
    table_lines = []
    in_table = False
    headers = None

    for line in lines:
        stripped = line.strip()
        if not in_table:
            if '|' in stripped and re.search(header_pattern, stripped):
                headers = [h.strip() for h in stripped.split('|')[1:-1]]
                in_table = True
                continue
        else:
            if stripped.startswith('|') and '---' in stripped:
                continue  # separator line
            elif stripped.startswith('|') and '|' in stripped[1:]:
                values = [v.strip() for v in stripped.split('|')[1:-1]]
                if len(values) == len(headers):
                    table_lines.append(dict(zip(headers, values)))
                elif len(values) < len(headers):
                    padded = values + [''] * (len(headers) - len(values))
                    table_lines.append(dict(zip(headers, padded)))
            else:
                in_table = False

    return table_lines

def main():
    with open('METHOD_AND_GATE_MANIFEST.md', 'r') as f:
        manifest = f.read()

    errors = []
    warnings = []

    # 1. Extract construction registry IDs
    const_rows = parse_table(manifest, r'ID.*Name.*Description')
    const_ids = set()
    for row in const_rows:
        rid = row.get('ID', '').strip()
        if rid.startswith('CONST-'):
            if rid in const_ids:
                errors.append(f"Duplicate construction ID: {rid}")
            const_ids.add(rid)

    # 2. Extract convention registry IDs
    conv_rows = parse_table(manifest, r'ID.*Name.*Source.*Value')
    conv_ids = set()
    for row in conv_rows:
        rid = row.get('ID', '').strip()
        if rid.startswith('CONV-'):
            if rid in conv_ids:
                errors.append(f"Duplicate convention ID: {rid}")
            conv_ids.add(rid)

    # 3. Extract gate registry IDs from all gate tables
    gate_ids = set()
    gate_tables = parse_table(manifest, r'ID.*Gate.*Establishes.*Mutation')
    for row in gate_tables:
        rid = row.get('ID', '').strip()
        if rid.startswith('GATE-'):
            if rid in gate_ids:
                errors.append(f"Duplicate gate ID: {rid}")
            gate_ids.add(rid)

    # 4. Extract coverage table
    coverage_rows = parse_table(manifest, r'Registry ID.*Checkable')
    coverage_ids = set()
    coverage_map = {}
    for row in coverage_rows:
        rid = row.get('Registry ID', '').strip()
        if not rid:
            continue
        if rid in coverage_ids:
            errors.append(f"Duplicate coverage table entry: {rid}")
        coverage_ids.add(rid)
        coverage_map[rid] = row

    # 5. Check set equality: (constructions ∪ conventions ∪ gates) == coverage IDs
    registry_ids = const_ids | conv_ids | gate_ids

    missing_from_coverage = registry_ids - coverage_ids
    extra_in_coverage = coverage_ids - registry_ids

    if missing_from_coverage:
        errors.append(f"Registry IDs missing from coverage table: {sorted(missing_from_coverage)}")
    if extra_in_coverage:
        errors.append(f"Coverage table IDs not in any registry: {sorted(extra_in_coverage)}")

    # 6. Check that every checkable row has PASS
    for rid, row in coverage_map.items():
        checkable = row.get('Checkable Pre-Impl?', '').strip()
        result = row.get('Result', '').strip()
        artifact = row.get('Validation Artifact', '').strip()

        if checkable == 'YES':
            if 'PASS' not in result:
                errors.append(f"{rid}: checkable=YES but result is '{result}', not PASS")
            if not artifact or artifact == '—':
                errors.append(f"{rid}: checkable=YES but no validation artifact named")
        elif checkable == 'NO':
            if result and result != '—':
                warnings.append(f"{rid}: marked not checkable but has result '{result}'")

    # 7. Check MANIFEST STATUS line
    if 'MANIFEST STATUS: FINAL' not in manifest:
        errors.append("Missing 'MANIFEST STATUS: FINAL' line")

    # Report
    print("=" * 60)
    print("MANIFEST VALIDATION REPORT")
    print("=" * 60)
    print(f"\nConstruction IDs: {sorted(const_ids)}")
    print(f"Convention IDs:   {sorted(conv_ids)}")
    print(f"Gate IDs:         {sorted(gate_ids)}")
    print(f"Coverage IDs:     {sorted(coverage_ids)}")
    print(f"\nTotal registry:  {len(registry_ids)}")
    print(f"Total coverage:  {len(coverage_ids)}")
    print(f"Set equality:    {registry_ids == coverage_ids}")

    if errors:
        print(f"\nERRORS ({len(errors)}):")
        for e in errors:
            print(f"  ✗ {e}")
    else:
        print("\nNo errors.")

    if warnings:
        print(f"\nWARNINGS ({len(warnings)}):")
        for w in warnings:
            print(f"  ⚠ {w}")

    if not errors:
        print("\n✓ MANIFEST VALIDATED: exact set equality confirmed, all checkable rows PASS.")
        return 0
    else:
        print("\n✗ MANIFEST VALIDATION FAILED")
        return 1

if __name__ == '__main__':
    sys.exit(main())
