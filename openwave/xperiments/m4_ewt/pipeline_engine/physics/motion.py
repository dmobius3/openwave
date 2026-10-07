"""
Wave-centre motion (plan Section 6, item 1.6).

WCMotionProcessor integrates each active WC's position forward under a
WCDriftRule. Two rules ship here: NoOpRule (F = 0, control) and
GradRhoRule (F = -grad rho, plan 2.7 B7a).

Integration is symplectic Euler:
    v += a * dt
    x += v * dt
It is not the field's leapfrog: a WC position is a coarse-grained
quantity, not a field value, and the plan does not tie WC motion to
the field stepper.

Rule contract:
    rule(ctx, wc_index) -> (fx, fy, fz)
The return value is interpreted as acceleration, grid units per
sim-time-squared. It is applied uniformly to every active WC by the
processor; the rule decides what acceleration each index gets.

Rules must be pure: no ctx.data mutation, no instance state written
inside __call__. A rule that writes to ctx.data or to itself breaks
the "every mutation goes through a processor" rule from the engine's
Z2. The tests in _test_motion.py pin the two shipped rules as pure.

Efficiency note: GradRhoRule reads the rho feature via to_numpy()
once per WC per step. That is fine for a handful of WCs and is the
kind of thing that gets optimised once a real run hits it. It is not
optimised here; there is no real run yet.
"""

from ..pipeline import BaseProcessor, Stage
from .features import EMCDensityField, WaveGrid
from .wc_types import WCState


class WCDriftRule:
    """
    Callable contract for a per-WC acceleration rule.

    Concrete rules override __call__. The processor invokes
    rule(ctx, wc_index) once per active WC per step.
    """

    def __call__(self, ctx, wc_index):
        raise NotImplementedError(
            f"{type(self).__name__} must implement __call__(ctx, wc_index)"
        )


class NoOpRule(WCDriftRule):
    """
    F = 0. Control rule.

    Every WC keeps its starting velocity. If the starting velocity is
    zero, the position does not move. This is the default rule on
    WCMotionProcessor.
    """

    def __call__(self, ctx, wc_index):
        return (0.0, 0.0, 0.0)


class GradRhoRule(WCDriftRule):
    """
    F = -grad rho, evaluated at the WC's nearest voxel.

    rho is the EMC packing density. Inside a soliton it is lower than
    the statutory background, so -grad rho points toward the deficit
    centre. A WC sitting outside the deficit drifts in; a WC sitting
    at the centre has zero gradient and stops.

    Reads EMCDensityField.rho. Reads nothing else. The nearest-voxel
    lookup uses round() and clamps the index to the interior so the
    central difference stays in range; the clamp is a test-only
    convenience, a production rule would want a proper boundary
    treatment.
    """

    def __call__(self, ctx, wc_index):
        grid = ctx.data.require(WaveGrid)
        wc_state = ctx.data.require(WCState)
        rho_field = ctx.data.require(EMCDensityField).rho
        wc = wc_state.centers[wc_index]

        i = int(round(wc.x))
        j = int(round(wc.y))
        k = int(round(wc.z))
        i = max(1, min(grid.nx - 2, i))
        j = max(1, min(grid.ny - 2, j))
        k = max(1, min(grid.nz - 2, k))

        rho = rho_field.to_numpy()
        gx = 0.5 * (rho[i + 1, j, k] - rho[i - 1, j, k])
        gy = 0.5 * (rho[i, j + 1, k] - rho[i, j - 1, k])
        gz = 0.5 * (rho[i, j, k + 1] - rho[i, j, k - 1])
        return (-float(gx), -float(gy), -float(gz))


class WCMotionProcessor(BaseProcessor):
    """
    Integrate WC motion under a drift rule. Stage.POST_UPDATE, order 5
    (before BoundaryProcessor at 10, so a boundary can clamp a WC if
    one is ever written; none is written here).

    Mutates WCState in place: per active WC,
        velocity += rule * dt
        position += velocity * dt
    Inactive WCs are skipped entirely; their velocity is not updated.

    The processor does not enforce a domain boundary; a WC under a
    persistent rule can leave the grid.

    Requires: WaveGrid, WCState, EMCDensityField. The last is required
    even for NoOpRule, so that rule swapping does not require rebuilding
    the pipeline. It is a union of what the shipped rules need; a rule
    that reads no field would want a lighter declaration, and if that
    becomes common the declaration can move to per-rule checks at
    setup() time. Not done yet; the pipeline's build-time validation is
    simpler with a fixed tuple.
    """

    name = "WCMotionProcessor"
    stage = Stage.POST_UPDATE
    order = 5
    requires = (WaveGrid, WCState, EMCDensityField)
    provides = ()

    def __init__(self, rule=None):
        self.rule = rule if rule is not None else NoOpRule()

    def process(self, ctx):
        wc_state = ctx.data.require(WCState)
        dt = ctx.sim.dt
        for idx, wc in enumerate(wc_state.centers):
            if not wc.active:
                continue
            ax, ay, az = self.rule(ctx, idx)
            vx, vy, vz = wc.velocity
            vx += ax * dt
            vy += ay * dt
            vz += az * dt
            wc.velocity = (vx, vy, vz)
            wc.x += vx * dt
            wc.y += vy * dt
            wc.z += vz * dt
            wc.force = (ax, ay, az)