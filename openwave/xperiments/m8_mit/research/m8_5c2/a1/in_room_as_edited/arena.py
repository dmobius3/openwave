"""§ 10 arena-constructor registry.

Every arena identifier, its constructor, and the nonlinear-permitted flag.
The registry keys on ARENA alone; the law is global c1=+1 with ONE named
exception: § 7's control-(i) full-width wiring run under c1=0 on A-R0-N36.
"""

AGREEMENT_RUNGS = [24, 32, 40, 48]
CONTROL_B_RUNGS = [36, 44, 52, 60]
ALL_RUNGS = sorted(set(AGREEMENT_RUNGS + CONTROL_B_RUNGS))
NONTRIVIAL_SECTORS = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8']


def build_arena_registry():
    """Build the § 10 arena-constructor registry as a list of dicts."""
    registry = []

    for N in ALL_RUNGS:
        registry.append({
            'id': f'A-R0-N{N}',
            'arena': f'E_R0 retained space at cutoff N={N}',
            'nonlinear_permitted': True,
            'constructor': f'R0 section basis levels 0..{N}, Galerkin projector P_N^R0',
        })

    for N in ALL_RUNGS:
        registry.append({
            'id': f'A-R0C2-N{N}',
            'arena': f'E_R0 ⊗ C² manufactured extension at cutoff N={N}',
            'nonlinear_permitted': True,
            'constructor': f'R0 section basis levels 0..{N} tensored with C², projector P_N^R0 ⊗ I_2',
        })

    registry.append({
        'id': 'A-CTRLA',
        'arena': 'Control A two-level full-S³ arena (levels 2 and 6), § 6',
        'nonlinear_permitted': True,
        'constructor': 'full-S³ scalar harmonics at levels {2, 6}, manufactured projector P_{2,6}',
    })

    registry.append({
        'id': 'A-CTRLI',
        'arena': 'control (i) c1=0 linear operator on H_{R0,12}, § 7',
        'nonlinear_permitted': True,
        'constructor': 'H_{R0,12} ≅ V_12^* (R0 level 12, dim 13), linear law c1=0',
    })

    for rho in NONTRIVIAL_SECTORS:
        for N in ALL_RUNGS:
            registry.append({
                'id': f'A-SECTOR-{rho}-N{N}',
                'arena': f'nontrivial W_{rho}-valued bases and projectors at N={N}',
                'nonlinear_permitted': False,
                'constructor': f'{rho} section basis levels 0..{N}, linear operations only',
            })

    return registry
