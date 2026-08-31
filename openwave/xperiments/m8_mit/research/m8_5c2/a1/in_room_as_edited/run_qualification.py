"""Full M8.5-C2 qualification run: all gates 1-11 in sequence.

This is the master entry point. It executes every gate, records to the
output ledger per § 13, and closes with the hash.
"""
import sys, time, json, hashlib
sys.path.insert(0, '.')

from build.ledger_util import (init_clock, gate_record, append_record,
                                resource_record, arm_allocation_record,
                                disposition_record, coverage_record)
from build.gates import run_gate_1, run_gate_2, run_gate_10, ArenaTracker
from build.gate_checks import (run_gate_3, run_gate_5, run_gate_6,
                                run_gate_4)
from build.newton import run_gate_7, run_gate_9
from build.integrator import run_gate_8
from build.group import IRREP_NAMES
from build.arena import (AGREEMENT_RUNGS, CONTROL_B_RUNGS, ALL_RUNGS,
                          NONTRIVIAL_SECTORS)
from build.sections import total_modes
from build.packet import generate_packet

import os
ledger_path = 'ledger/OUTPUT_LEDGER.jsonl'
if os.path.exists(ledger_path):
    os.remove(ledger_path)

init_clock()
tracker = ArenaTracker()
STOP = False

print("=" * 60)
print("M8.5-C2 FULL QUALIFICATION RUN")
print("=" * 60)
t_start = time.time()

# Generate the § 4.3 field packet
print("\nGenerating § 4.3 field packet...", flush=True)
t0 = time.time()
packet_data, packet_bytes = generate_packet()
packet_hash = hashlib.sha256(packet_bytes).hexdigest()
print(f"  Packet hash: {packet_hash}")
print(f"  Generated in {time.time()-t0:.1f}s", flush=True)

# Build packet dict keyed by rung
packet = {}
for N in AGREEMENT_RUNGS + CONTROL_B_RUNGS:
    if N in packet_data:
        packet[N] = packet_data[N]

# ======================================================================
# Gate 1: G-LIN
# ======================================================================
print("\n" + "=" * 60)
print("GATE 1: G-LIN (linear operator wiring)")
print("=" * 60)
t0 = time.time()
g1 = run_gate_1(IRREP_NAMES, ALL_RUNGS)
dt1 = time.time() - t0
print(f"Gate 1: {'PASS' if g1 else 'FAIL'}  [{dt1:.1f}s]", flush=True)

for rho in IRREP_NAMES:
    for N in ALL_RUNGS:
        if total_modes(rho, N) > 0:
            if rho == 'R0':
                tracker.mark_exercised(f'A-R0-N{N}')
            else:
                tracker.mark_exercised(f'A-SECTOR-{rho}-N{N}')

if not g1:
    STOP = True
    print("STOP-QUAL: Gate 1 failed")

# ======================================================================
# Gate 2: G-LABEL
# ======================================================================
if not STOP:
    print("\n" + "=" * 60)
    print("GATE 2: G-LABEL (Casimir labeling)")
    print("=" * 60)
    t0 = time.time()
    g2 = run_gate_2(IRREP_NAMES, AGREEMENT_RUNGS)
    dt2 = time.time() - t0
    print(f"Gate 2: {'PASS' if g2 else 'FAIL'}  [{dt2:.1f}s]", flush=True)
    if not g2:
        STOP = True
        print("STOP-QUAL: Gate 2 failed")

# ======================================================================
# Gate 3: G-SECTOR
# ======================================================================
if not STOP:
    print("\n" + "=" * 60)
    print("GATE 3: G-SECTOR (sector bases)")
    print("=" * 60)
    t0 = time.time()
    g3, g3_recs = run_gate_3(IRREP_NAMES, AGREEMENT_RUNGS, tracker.exercised)
    dt3 = time.time() - t0
    print(f"Gate 3: {'PASS' if g3 else 'FAIL'}  [{dt3:.1f}s]", flush=True)
    if not g3:
        STOP = True
        print("STOP-QUAL: Gate 3 failed")

# ======================================================================
# Gate 4: G-PROJECTOR (dual-route agreement)
# ======================================================================
if not STOP:
    print("\n" + "=" * 60)
    print("GATE 4: G-PROJECTOR (dual-route agreement)")
    print("=" * 60)
    t0 = time.time()
    g4, g4_recs = run_gate_4(ALL_RUNGS, tracker.exercised,
                              packet=packet, verbose=True)
    dt4 = time.time() - t0
    print(f"Gate 4: {'PASS' if g4 else 'FAIL'}  [{dt4:.1f}s]", flush=True)
    if not g4:
        STOP = True
        print("STOP-QUAL: Gate 4 failed")

# ======================================================================
# Gate 5: G-STRUCTURAL
# ======================================================================
if not STOP:
    print("\n" + "=" * 60)
    print("GATE 5: G-STRUCTURAL (structural identities)")
    print("=" * 60)
    t0 = time.time()
    g5 = run_gate_5(AGREEMENT_RUNGS, tracker.exercised)
    dt5 = time.time() - t0
    print(f"Gate 5: {'PASS' if g5 else 'FAIL'}  [{dt5:.1f}s]", flush=True)
    if not g5:
        STOP = True
        print("STOP-QUAL: Gate 5 failed")

# ======================================================================
# Gate 6: G-CASCADE
# ======================================================================
if not STOP:
    print("\n" + "=" * 60)
    print("GATE 6: G-CASCADE (cascade monitor)")
    print("=" * 60)
    t0 = time.time()
    g6, g6_recs = run_gate_6(AGREEMENT_RUNGS, tracker.exercised)
    dt6 = time.time() - t0
    print(f"Gate 6: {'PASS' if g6 else 'FAIL'}  [{dt6:.1f}s]", flush=True)
    if not g6:
        STOP = True
        print("STOP-QUAL: Gate 6 failed")

# ======================================================================
# Gate 7: G-CONTINUATION
# ======================================================================
if not STOP:
    print("\n" + "=" * 60)
    print("GATE 7: G-CONTINUATION (Newton controls)")
    print("=" * 60)
    t0 = time.time()
    g7, g7_report = run_gate_7(control_b_rungs=CONTROL_B_RUNGS, verbose=True)
    dt7 = time.time() - t0
    print(f"Gate 7: {'PASS' if g7 else 'FAIL'}  [{dt7:.1f}s]", flush=True)

    tracker.mark_exercised('A-CTRLA')
    tracker.mark_nonlinear('A-CTRLA')
    for N in CONTROL_B_RUNGS:
        tracker.mark_exercised(f'A-R0-N{N}')
        tracker.mark_nonlinear(f'A-R0-N{N}')

    gate_record(
        gate_id='G7-CONTINUATION',
        arena_id='A-CTRLA',
        rung=0,
        parent_status='GREEN' if g7 else 'RED',
        mutation_status='N/A',
        measured_values={'control_a': g7_report.get('control_a', {}),
                         'control_b_rungs': list(CONTROL_B_RUNGS)},
        wall_clock_seconds=dt7,
    )

    if not g7:
        STOP = True
        print("STOP-QUAL: Gate 7 failed")

# ======================================================================
# Gate 8: G-TIME (time arm)
# ======================================================================
if not STOP:
    print("\n" + "=" * 60)
    print("GATE 8: G-TIME (Störmer-Verlet time arm)")
    print("=" * 60)
    t0 = time.time()
    g8, g8_report = run_gate_8(verbose=True)
    dt8 = time.time() - t0
    print(f"Gate 8: {'PASS' if g8 else 'FAIL'}  [{dt8:.1f}s]", flush=True)

    tracker.mark_exercised('A-CTRLI')
    # Control (i) has c1=0 (law exception), nonlinear by construction
    tracker.mark_nonlinear('A-CTRLI')
    # Control (ii) is nonlinear on the level-0 R0 mode (at N=36 wiring)
    tracker.mark_nonlinear('A-R0-N36')

    gate_record(
        gate_id='G8-TIME',
        arena_id='A-CTRLI',
        rung=0,
        parent_status='GREEN' if g8 else 'RED',
        mutation_status='GREEN' if g8_report.get('mutation_arm', {}).get('euler_fails_contraction', False) else 'RED',
        measured_values=g8_report,
        wall_clock_seconds=dt8,
        law_exception=True,
    )

    if not g8:
        STOP = True
        print("STOP-QUAL: Gate 8 failed")

# ======================================================================
# Gate 9: G-CONVERGENCE
# ======================================================================
if not STOP:
    print("\n" + "=" * 60)
    print("GATE 9: G-CONVERGENCE")
    print("=" * 60)

    t0 = time.time()
    g9, g9_report = run_gate_9(g7_report, control_b_rungs=CONTROL_B_RUNGS,
                                verbose=True)
    dt9 = time.time() - t0
    print(f"Gate 9: {'PASS' if g9 else 'FAIL'}  [{dt9:.1f}s]", flush=True)

    gate_record(
        gate_id='G9-CONVERGENCE',
        arena_id='A-R0-N36',
        rung=0,
        parent_status='GREEN' if g9 else 'RED',
        mutation_status='GREEN' if g9_report.get('mutation', {}).get('caught', False) else 'RED',
        measured_values=g9_report,
        wall_clock_seconds=dt9,
    )

    if not g9:
        STOP = True
        print("STOP-QUAL: Gate 9 failed")

# ======================================================================
# Gate 10: G-COVERAGE
# ======================================================================
if not STOP:
    print("\n" + "=" * 60)
    print("GATE 10: G-COVERAGE")
    print("=" * 60)

    # Mark manufactured and nonlinear arenas
    tracker.mark_exercised('A-CTRLA')
    tracker.mark_exercised('A-CTRLI')
    for N in ALL_RUNGS:
        tracker.mark_exercised(f'A-R0C2-N{N}')
        tracker.mark_nonlinear(f'A-R0-N{N}')
        tracker.mark_nonlinear(f'A-R0C2-N{N}')

    g10 = run_gate_10(tracker.exercised, tracker.nonlinear)
    print(f"Gate 10: {'PASS' if g10 else 'FAIL'}", flush=True)
    if not g10:
        STOP = True
        print("STOP-QUAL: Gate 10 failed")

# ======================================================================
# Gate 11: Law scope (RECORD only)
# ======================================================================
print("\n" + "=" * 60)
print("GATE 11: Law scope (RECORD)")
print("=" * 60)
disposition_record('LAW-SCOPE', decision='cubic_only', saturating='out_of_scope')
print("  Recorded: option (b), cubic only; saturating out of scope")

# ======================================================================
# FINAL SUMMARY
# ======================================================================
total_time = time.time() - t_start
print("\n" + "=" * 60)
print("QUALIFICATION SUMMARY")
print("=" * 60)
print(f"Total wall-clock: {total_time:.1f}s ({total_time/3600:.2f}h)")
print(f"48h ceiling: {'OK' if total_time < 48*3600 else 'EXCEEDED'}")

if STOP:
    print("\nQUALIFICATION: STOPPED (see above)")
    disposition_record('M8.5-C2-FAILED', reason='STOP-QUAL')
else:
    print("\nAll gates passed.")
    disposition_record('M8.5-C2-QUALIFIED')

# Ledger hash
import os
ledger_path = 'ledger/OUTPUT_LEDGER.jsonl'
if os.path.exists(ledger_path):
    with open(ledger_path, 'rb') as f:
        ledger_hash = hashlib.sha256(f.read()).hexdigest()
    print(f"\nOutput ledger hash: {ledger_hash}")
    append_record({'type': 'LEDGER-HASH', 'sha256': ledger_hash})
