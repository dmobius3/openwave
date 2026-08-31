"""Execute all qualification gates for the M8.5-C2 commission."""
import sys, time, json
sys.path.insert(0, '.')

from build.ledger_util import init_clock, gate_record, append_record
from build.gates import run_gate_1, run_gate_2, run_gate_10, ArenaTracker
from build.gate_checks import run_gate_3, run_gate_5, run_gate_6
from build.group import IRREP_NAMES
from build.arena import AGREEMENT_RUNGS, CONTROL_B_RUNGS, ALL_RUNGS, NONTRIVIAL_SECTORS

init_clock()
tracker = ArenaTracker()

print("=" * 60)
print("M8.5-C2 GATE EXECUTION")
print("=" * 60)

# Gate 1: G-LIN
t0 = time.time()
g1 = run_gate_1(IRREP_NAMES, ALL_RUNGS)
dt = time.time() - t0
print(f"Gate 1 (G-LIN):    {'PASS' if g1 else 'FAIL'}  [{dt:.1f}s]", flush=True)

# Mark all arenas exercised by Gate 1
for rho in IRREP_NAMES:
    for N in ALL_RUNGS:
        from build.sections import total_modes
        if total_modes(rho, N) > 0:
            if rho == 'R0':
                tracker.mark_exercised(f'A-R0-N{N}')
            else:
                tracker.mark_exercised(f'A-SECTOR-{rho}-N{N}')

# Gate 2: G-LABEL
t0 = time.time()
g2 = run_gate_2(IRREP_NAMES, AGREEMENT_RUNGS)
dt = time.time() - t0
print(f"Gate 2 (G-LABEL):  {'PASS' if g2 else 'FAIL'}  [{dt:.1f}s]", flush=True)

# Gate 3: G-SECTOR
t0 = time.time()
g3, g3_recs = run_gate_3(IRREP_NAMES, AGREEMENT_RUNGS, tracker.exercised)
dt = time.time() - t0
print(f"Gate 3 (G-SECTOR): {'PASS' if g3 else 'FAIL'}  [{dt:.1f}s]", flush=True)

# Gate 5: Structural identities (R0 only)
t0 = time.time()
g5 = run_gate_5(AGREEMENT_RUNGS, tracker.exercised)
dt = time.time() - t0
print(f"Gate 5 (G-STRUCT): {'PASS' if g5 else 'FAIL'}  [{dt:.1f}s]", flush=True)

# Gate 6: Cascade monitor
t0 = time.time()
g6 = run_gate_6(AGREEMENT_RUNGS, tracker.exercised)
dt = time.time() - t0
print(f"Gate 6 (G-CASCADE):{'PASS' if g6 else 'FAIL'}  [{dt:.1f}s]", flush=True)

# Gate 10: Coverage check (preliminary — will be finalized after gates 4, 7, 8, 9)
# Mark manufactured arenas
tracker.mark_exercised('A-CTRLA')
tracker.mark_exercised('A-CTRLI')
for N in AGREEMENT_RUNGS:
    tracker.mark_exercised(f'A-R0C2-N{N}')
for N in CONTROL_B_RUNGS:
    tracker.mark_exercised(f'A-R0-N{N}')
    tracker.mark_exercised(f'A-R0C2-N{N}')
    for rho in NONTRIVIAL_SECTORS:
        tracker.mark_exercised(f'A-SECTOR-{rho}-N{N}')

print()
print(f"Exercised arenas: {len(tracker.exercised)}")
print(f"Nonlinear arenas: {len(tracker.nonlinear)}")

results = {
    'G1': g1, 'G2': g2, 'G3': g3,
    'G5': g5, 'G6': g6,
}
all_pass = all(results.values())
print()
for k, v in results.items():
    print(f"  {k}: {'PASS' if v else 'FAIL'}")
print(f"\nOverall (gates 1-3, 5-6): {'ALL PASS' if all_pass else 'FAILURES DETECTED'}")
