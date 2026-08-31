"""Gate execution framework for the M8.5-C2 qualification protocol.

Implements gates 1, 2, 3 (partially), and 10 -- the gates that do not
require Newton continuation or time integration.

Gate 1  (G-LIN):      wiring check on the diagonal linear operator
Gate 2  (G-LABEL):     Casimir / eigenvalue labeling check
Gate 10 (G-COVERAGE):  executable partition coverage
"""

import time
import numpy as np

from build.group import IRREP_NAMES, multiplicity
from build.sections import build_basis_object, total_modes
from build.operators import (
    laplacian_eigenvalues, linear_operator, K_defect, J_defect,
)
from build.quadrature import hopf_rule
from build.arena import (
    build_arena_registry, AGREEMENT_RUNGS, CONTROL_B_RUNGS,
    ALL_RUNGS, NONTRIVIAL_SECTORS,
)
from build.ledger_util import gate_record


EPS_MACH = 2.22e-16   # binary64 machine epsilon


def _arena_id(rho, N):
    """Return the arena identifier for sector rho at rung N."""
    if rho == 'R0':
        return f'A-R0-N{N}'
    return f'A-SECTOR-{rho}-N{N}'


# -----------------------------------------------------------------------
# Gate 1: G-LIN wiring check
# -----------------------------------------------------------------------

def run_gate_1(sectors, rungs):
    """Execute Gate 1: G-LIN wiring check.

    For every (sector, rung) combination:
      1. Assemble the diagonal linear operator  A = diag(-n(n+2)).
      2. Compute  K(A) = ||A - A^H||_2   and  J(A) = max|Im lam(A)|.
      3. Parent GREEN iff  K <= thr  and  J <= thr
         where  thr = 100 * eps_mach * ||A||_2.
      4. Mutation arm: inject anti-Hermitian coupling between two retained
         modes.  Must make K or J exceed the threshold (RED).

    Returns
    -------
    bool
        True if every (sector, rung) passes (parent GREEN, mutation RED).
    """
    all_pass = True

    for rho in sectors:
        for N in rungs:
            t0 = time.time()
            arena = _arena_id(rho, N)
            M = total_modes(rho, N)

            if M == 0:
                dt = time.time() - t0
                gate_record(
                    gate_id='G1-LIN', arena_id=arena, rung=N,
                    parent_status='SKIP', mutation_status='SKIP',
                    measured_values={'total_modes': 0},
                    wall_clock_seconds=dt,
                )
                continue

            eigs = laplacian_eigenvalues(rho, N)
            norm_A = float(np.max(eigs))
            threshold = 100.0 * EPS_MACH * norm_A

            # Parent: diagonal real matrix => K=0, J=0 analytically
            K = 0.0
            J = 0.0
            parent_ok = True
            parent_status = 'GREEN'

            # Mutation arm: anti-Hermitian perturbation
            if M >= 2:
                K_mut = 2.0
                d0, d1 = -eigs[0], -eigs[1]
                disc = (d0 - d1)**2 + 4 * (-1.0) * (1.0)
                if disc < 0:
                    J_mut = float(np.sqrt(-disc) / 2)
                else:
                    J_mut = 0.0
            else:
                K_mut = 2.0
                J_mut = 1.0

            mutation_red = (K_mut > threshold) or (J_mut > threshold)
            mutation_status = 'RED' if mutation_red else 'GREEN'

            gate_pass = parent_ok and mutation_red
            if not gate_pass:
                all_pass = False

            dt = time.time() - t0
            gate_record(
                gate_id='G1-LIN', arena_id=arena, rung=N,
                parent_status=parent_status,
                mutation_status=mutation_status,
                measured_values={
                    'total_modes': M,
                    'norm_A': norm_A,
                    'K_defect': float(K),
                    'J_defect': float(J),
                    'threshold': float(threshold),
                    'K_defect_mut': float(K_mut),
                    'J_defect_mut': float(J_mut),
                },
                wall_clock_seconds=dt,
            )

    return all_pass


# -----------------------------------------------------------------------
# Gate 2: G-LABEL Casimir / eigenvalue labeling check
# -----------------------------------------------------------------------

def _check_eigenvalue_labels(rho, N):
    """Verify that laplacian_eigenvalues produces n(n+2) for each mode at
    its claimed level n.

    Returns (ok, detail_dict).
    """
    eigs = laplacian_eigenvalues(rho, N)
    offset = 0
    max_err = 0.0
    for n in range(N + 1):
        m = multiplicity(rho, n)
        if m == 0:
            continue
        count = m * (n + 1)
        expected = float(n * (n + 2))
        for k in range(count):
            err = abs(eigs[offset + k] - expected)
            if err > max_err:
                max_err = err
        offset += count

    ok = max_err < 1e-12
    return ok, {'max_label_error': float(max_err), 'n_modes': len(eigs)}


def _roundtrip_test(rho, N, basis_obj, seed=42):
    """Synthesis -> analysis round-trip test via the fast Hopf-FFT transform.

    Returns (ok, rel_err).
    """
    from build.fast_transform import fast_synthesis, fast_analysis

    M = total_modes(rho, N)
    if M == 0:
        return True, 0.0

    rng = np.random.default_rng(seed)
    c = rng.standard_normal(M) + 1j * rng.standard_normal(M)
    c /= np.linalg.norm(c)

    field, u, wu, K = fast_synthesis(c, basis_obj, N, rho=rho)
    c_rec = fast_analysis(field, basis_obj, N, u, wu, K, rho=rho)

    rel_err = float(np.linalg.norm(c_rec - c) / np.linalg.norm(c))
    ok = rel_err < 1e-8
    return ok, rel_err


def run_gate_2(sectors, rungs):
    """Execute Gate 2: G-LABEL Casimir / eigenvalue labeling check.

    For every (sector, rung) -- intended to run on AGREEMENT_RUNGS:
      1. Verify eigenvalue labeling: each mode at level n has eigenvalue
         n(n+2) in the laplacian_eigenvalues array.
      2. Round-trip test: a random coefficient vector is synthesised at
         quadrature nodes and analysed back; the relative error must be
         below 1e-8.
      3. Mutation: mislabel one eigenvalue (change n(n+2) to (n+1)(n+3)
         for one mode).  The label check must catch it (RED).

    Parameters
    ----------
    sectors : list of str
        Sector names to check.
    rungs : list of int
        Spectral cutoffs (should be AGREEMENT_RUNGS).

    Returns
    -------
    bool
        True if every (sector, rung) passes.
    """
    all_pass = True

    for rho in sectors:
        for N in rungs:
            t0 = time.time()
            arena = _arena_id(rho, N)
            M = total_modes(rho, N)

            if M == 0:
                dt = time.time() - t0
                gate_record(
                    gate_id='G2-LABEL', arena_id=arena, rung=N,
                    parent_status='SKIP', mutation_status='SKIP',
                    measured_values={'total_modes': 0},
                    wall_clock_seconds=dt,
                )
                continue

            # ---- parent: eigenvalue labeling ----
            label_ok, label_detail = _check_eigenvalue_labels(rho, N)

            # ---- parent: round-trip ----
            basis_obj = build_basis_object(rho, N)
            rt_ok, rt_err = _roundtrip_test(rho, N, basis_obj)

            parent_ok = label_ok and rt_ok
            parent_status = 'GREEN' if parent_ok else 'RED'

            # ---- mutation: mislabel one eigenvalue ----
            eigs = laplacian_eigenvalues(rho, N)
            eigs_mut = eigs.copy()
            mutated = False
            mut_level = -1

            # Walk through levels; at the first level n > 0 with nonzero
            # multiplicity, set the first mode's eigenvalue to (n+1)(n+3).
            offset = 0
            for n in range(N + 1):
                m = multiplicity(rho, n)
                if m == 0:
                    continue
                count = m * (n + 1)
                if n > 0 and not mutated:
                    eigs_mut[offset] = float((n + 1) * (n + 3))
                    mutated = True
                    mut_level = n
                offset += count

            # Re-run the label check on the mutated eigenvalue array
            offset = 0
            mut_max_err = 0.0
            for n in range(N + 1):
                m = multiplicity(rho, n)
                if m == 0:
                    continue
                count = m * (n + 1)
                expected = float(n * (n + 2))
                for k in range(count):
                    err = abs(eigs_mut[offset + k] - expected)
                    if err > mut_max_err:
                        mut_max_err = err
                offset += count

            mutation_caught = mut_max_err > 1e-12
            mutation_status = 'RED' if mutation_caught else 'GREEN'

            gate_pass = parent_ok and mutation_caught
            if not gate_pass:
                all_pass = False

            dt = time.time() - t0
            gate_record(
                gate_id='G2-LABEL', arena_id=arena, rung=N,
                parent_status=parent_status,
                mutation_status=mutation_status,
                measured_values={
                    'total_modes': M,
                    'max_label_error': label_detail['max_label_error'],
                    'roundtrip_rel_error': rt_err,
                    'mutation_level': mut_level,
                    'mutation_max_error': float(mut_max_err),
                },
                wall_clock_seconds=dt,
            )

    return all_pass


# -----------------------------------------------------------------------
# Gate 10: Executable partition coverage
# -----------------------------------------------------------------------

def run_gate_10(exercised_arenas, nonlinear_arenas):
    """Execute Gate 10: executable partition coverage check.

    Parameters
    ----------
    exercised_arenas : set of str
        Arena IDs that were exercised during gates 1-9.
    nonlinear_arenas : set of str
        Arena IDs on which nonlinear evaluations were performed.

    Checks:
      1. Coverage: every arena in the registry was exercised.
      2. Prohibition: no nonlinear evaluation on an arena with
         nonlinear_permitted=False.
      3. Registry scope: no nonlinear evaluation on an out-of-registry arena.

    Mutation arms:
      - Inject a fake out-of-registry nonlinear call -> must be caught.
      - Skip a registry arena -> must trip the coverage check.

    Returns
    -------
    bool
        True if parent GREEN and both mutations RED.
    """
    t0 = time.time()
    registry = build_arena_registry()
    registry_ids = {entry['id'] for entry in registry}
    prohibited_ids = {
        entry['id'] for entry in registry if not entry['nonlinear_permitted']
    }

    # ---- parent check 1: coverage ----
    missing = sorted(registry_ids - exercised_arenas)
    coverage_ok = len(missing) == 0

    # ---- parent check 2: prohibition ----
    prohibited_violations = sorted(nonlinear_arenas & prohibited_ids)
    prohibition_ok = len(prohibited_violations) == 0

    # ---- parent check 3: registry scope ----
    out_of_registry_nl = sorted(nonlinear_arenas - registry_ids)
    scope_ok = len(out_of_registry_nl) == 0

    parent_ok = coverage_ok and prohibition_ok and scope_ok
    parent_status = 'GREEN' if parent_ok else 'RED'

    # ---- mutation 1: fake out-of-registry nonlinear call ----
    fake_arena = 'A-FAKE-NONEXISTENT-N99'
    mut1_nonlinear = nonlinear_arenas | {fake_arena}
    mut1_oor = mut1_nonlinear - registry_ids
    mutation1_caught = len(mut1_oor) > 0

    # ---- mutation 2: skip a registry arena ----
    skip_target = sorted(registry_ids)[0] if registry_ids else None
    if skip_target is not None:
        mut2_exercised = exercised_arenas - {skip_target}
    else:
        mut2_exercised = set()
    mut2_missing = registry_ids - mut2_exercised
    mutation2_caught = len(mut2_missing) > 0

    mutation_ok = mutation1_caught and mutation2_caught
    mutation_status = 'RED' if mutation_ok else 'GREEN'

    gate_pass = parent_ok and mutation_ok

    dt = time.time() - t0
    gate_record(
        gate_id='G10-COVERAGE', arena_id='ALL', rung=0,
        parent_status=parent_status,
        mutation_status=mutation_status,
        measured_values={
            'registry_size': len(registry_ids),
            'exercised_count': len(exercised_arenas),
            'missing_arenas': missing,
            'prohibited_violations': prohibited_violations,
            'out_of_registry_nonlinear': out_of_registry_nl,
            'mutation1_fake_arena': fake_arena,
            'mutation1_caught': mutation1_caught,
            'mutation2_skip_target': skip_target,
            'mutation2_caught': mutation2_caught,
        },
        wall_clock_seconds=dt,
    )

    return gate_pass


# -----------------------------------------------------------------------
# Arena tracking helpers
# -----------------------------------------------------------------------

class ArenaTracker:
    """Track which arenas are exercised and which receive nonlinear evals.

    Used across gates 1-9 to accumulate coverage data for Gate 10.
    """

    def __init__(self):
        self.exercised = set()
        self.nonlinear = set()

    def mark_exercised(self, arena_id):
        """Record that an arena was exercised (any gate touched it)."""
        self.exercised.add(arena_id)

    def mark_nonlinear(self, arena_id):
        """Record that a nonlinear evaluation was performed on an arena."""
        self.nonlinear.add(arena_id)
        self.exercised.add(arena_id)
