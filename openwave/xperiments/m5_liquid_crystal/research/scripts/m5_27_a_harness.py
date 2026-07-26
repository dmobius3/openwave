"""M5.27 arm A: the driven-V4 harness (phase A, prescribed uniform background).

THE DRIVE (host-side, no engine surgery). The verified-L spectral potential is

    V4(M) = w Sum_{p=1..4} ( tr((M eta)^p) - C_p )^2 ,   C_p = sg^p + 1 + delta^p

and ALL of its g-dependence sits in the trace targets C_p through `sg`, which
`v4_of` / `dv4_of` / `evolve_M_eta_finish` take as a per-call ARGUMENT
(engine2_pde.py). The prescribed background scalar is therefore imposed by
passing a time-dependent sample every step:

    chi(t) = A cos(om_bar t),  A = 1        (only the product kappa*A enters)
    sg(t)  = g + kappa chi(t) = g (1 + eps cos(om_bar t)),   eps = kappa/g

Production engine files are consumed READ-ONLY; everything new lives here.

Contents:
  - `Harness`: field + state loading + the driven leapfrog step
  - the drive-power ledger  P(t) = sum_cells dV4/dsg * dsg/dt   (analytic kernel)
  - the clock-phase instrument (unwrapped, apolar mod-pi handled)  [blindspot B9]
  - boundary handling for the driven box (adiabatic vacuum tracking) [B7]
  - the staged-4x4 mixed-block projection with leak booking

Imported by the gate / scan / re-read scripts; not run directly.
"""
# NOTE: no `from __future__ import annotations` here — PEP 563 turns annotations
# into strings, which breaks taichi's ti.template() kernel-argument typing.
import math
import os
import sys

import numpy as np
import taichi as ti

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


_TI_ARCH = None


def init_taichi(prefer_gpu=True):
    """Init taichi once, GPU (Metal on the M4) when available: the scan needs
    ~1e5 steps per run point, so the arch choice sets the whole task budget."""
    global _TI_ARCH
    if _TI_ARCH is not None:
        return _TI_ARCH
    arch = ti.gpu if prefer_gpu else ti.cpu
    try:
        ti.init(arch=arch, log_level=ti.WARN, random_seed=0, default_fp=ti.f32)
        _TI_ARCH = str(ti.cfg.arch)
    except Exception:
        ti.init(arch=ti.cpu, log_level=ti.WARN, random_seed=0, default_fp=ti.f32)
        _TI_ARCH = "cpu(fallback)"
    return _TI_ARCH


# ---- registered constants (frozen at go; see task_details TASK PLANNING) ----
H_RES = 1.5          # research grid unit h (the m5_21_2b/3/9 arena)
DT = 0.005           # dt_eff, the certified canonical timestep
RENV = 10.0          # clock-envelope radius for set_fixed_j / read_carried_j
G0 = 8.0             # medium.LC_G, the undriven time-axis eigenvalue
DELTA = 0.5          # medium.LC_DELTA
N_GRID = 31          # 32^3 source npz -> 31^3 crop (load_npz_M crop path)

FIXJ_OMSTAR = {
    "conj_om0.2": 0.1992308983677509,
    "conj_om0.5": 0.49655289701674865,
    "conj_om1": 0.9820390008334264,
}


def fixj_path(tag: str) -> str:
    return os.path.join(DATA, f"m5_21_9_fixedj_{tag}_end.npz")


# ================================================================
# the drive (host-side)
# ================================================================
RAMP_T = 60.0   # drive-amplitude ramp DURATION in time units (~2 clock periods
#                 at om* = 0.199, whose period is 31.6). Registered in TIME, not
#                 in drive cycles: a cycle-based ramp scales with 1/om_bar, so at
#                 the low-om_bar end it consumed most of the run (deviation
#                 logged 2026-07-24 after the P1 probe).


def sg_of(t: float, eps: float, om_bar: float, ramp_t: float = RAMP_T) -> float:
    """sg(t) = g (1 + eps_eff(t) cos(om_bar t)), with an ADIABATIC RAMP of the
    amplitude over `ramp_t` TIME UNITS (blindspot B8: a step switch-on shocks
    C_p and launches a transient that can fake or mask a capture)."""
    if om_bar <= 0.0 or eps == 0.0:
        return G0
    env = 1.0 if (t >= ramp_t or ramp_t <= 0.0) else \
        0.5 * (1.0 - math.cos(math.pi * t / ramp_t))
    return G0 * (1.0 + eps * env * math.cos(om_bar * t))


def dsg_dt(t: float, eps: float, om_bar: float, ramp_t: float = RAMP_T,
           hstep: float = 1e-4) -> float:
    """d sg/dt by central difference on the SAME sg_of the stepper samples
    (keeps the ledger consistent with the ramp; no separate analytic form)."""
    return (sg_of(t + hstep, eps, om_bar, ramp_t)
            - sg_of(t - hstep, eps, om_bar, ramp_t)) / (2.0 * hstep)


# ================================================================
# engine imports (read-only consumption)
# ================================================================
from openwave.xperiments.m5_liquid_crystal import engine1_seeds as seeds  # noqa: E402
from openwave.xperiments.m5_liquid_crystal import engine2_pde as pde  # noqa: E402
from openwave.xperiments.m5_liquid_crystal.medium import TensorField  # noqa: E402

W1 = pde.W1_SPECTRAL
eta_right = pde.eta_right


# ================================================================
# new taichi kernels (accumulate atomically into 0-D fields)
# ================================================================
@ti.kernel
def dv4_dsg_sum_k(tf: ti.template(), sg: ti.f32, delta: ti.f32,  # type: ignore
                  w1: ti.f32, out: ti.template()):  # type: ignore
    """Sum over cells of dV4/dsg, exact:
        dV4/dsg = -2 w Sum_p p sg^(p-1) (t_p - C_p)
    (the traces t_p = tr((M eta)^p) carry no sg dependence)."""
    c1 = sg + 1.0 + delta
    c2 = sg * sg + 1.0 + delta * delta
    c3 = sg * sg * sg + 1.0 + delta * delta * delta
    c4 = sg * sg * sg * sg + 1.0 + delta * delta * delta * delta
    for i, j, k in tf.M_am:
        me = eta_right(tf.M_am[i, j, k])
        p2 = me @ me
        p3 = p2 @ me
        t1 = me.trace()
        t2 = p2.trace()
        t3 = p3.trace()
        t4 = (p3 @ me).trace()
        out[None] += -2.0 * w1 * (
            (t1 - c1)
            + 2.0 * sg * (t2 - c2)
            + 3.0 * sg * sg * (t3 - c3)
            + 4.0 * sg * sg * sg * (t4 - c4)
        )


@ti.kernel
def v4_sum_k(tf: ti.template(), sg: ti.f32, delta: ti.f32,  # type: ignore
             w1: ti.f32, out: ti.template()):  # type: ignore
    """Sum over cells of V4 (the potential-energy ledger line)."""
    c1 = sg + 1.0 + delta
    c2 = sg * sg + 1.0 + delta * delta
    c3 = sg * sg * sg + 1.0 + delta * delta * delta
    c4 = sg * sg * sg * sg + 1.0 + delta * delta * delta * delta
    for i, j, k in tf.M_am:
        me = eta_right(tf.M_am[i, j, k])
        p2 = me @ me
        p3 = p2 @ me
        out[None] += w1 * (
            (me.trace() - c1) ** 2
            + (p2.trace() - c2) ** 2
            + (p3.trace() - c3) ** 2
            + ((p3 @ me).trace() - c4) ** 2
        )


@ti.kernel
def kinetic_sum_k(tf: ti.template(), dt_eff: ti.f32, out: ti.template()):  # type: ignore
    """Sum of 1/2 ||Mdot||_F^2, Mdot = (M - M_prev)/dt."""
    inv = 1.0 / dt_eff
    for i, j, k in tf.M_am:
        md = (tf.M_am[i, j, k] - tf.M_prev_am[i, j, k]) * inv
        out[None] += 0.5 * (md * md).sum()


@ti.kernel
def project_mixed_k(tf: ti.template(), leak: ti.template()):  # type: ignore
    """The STAGED 4x4: zero the (0,i) mixed block in M and M_prev, recording the
    largest removed magnitude (atomic max) so the deferral is BOOKED, never
    silently dropped. The time-time entry (0,0) and the spatial block stay live."""
    for i, j, k in tf.M_am:
        m = tf.M_am[i, j, k]
        mp = tf.M_prev_am[i, j, k]
        for a in ti.static(range(1, 4)):
            ti.atomic_max(leak[None], ti.abs(m[0, a]))
            m[0, a] = 0.0
            m[a, 0] = 0.0
            mp[0, a] = 0.0
            mp[a, 0] = 0.0
        tf.M_am[i, j, k] = m
        tf.M_prev_am[i, j, k] = mp


@ti.kernel
def snap_time_entry_k(tf: ti.template(), store: ti.template()):  # type: ignore
    """Store the t = 0 time-time entry (the boundary tracker's reference)."""
    for i, j, k in tf.M_am:
        store[i, j, k] = tf.M_am[i, j, k][0, 0]


@ti.kernel
def track_boundary_k(tf: ti.template(), store: ti.template(), ratio: ti.f32):  # type: ignore
    """B7 guard: the interior stepper never updates boundary cells, so under a
    global drive a pinned boundary stays at the UNDRIVEN vacuum while the
    interior tracks -sg(t) — a mismatch layer that radiates. Track the drive on
    the boundary shell analytically by scaling ONLY its time-time entry by
    sg(t)/sg(0), preserving whatever structure the loaded state has there."""
    nx, ny, nz = tf.nx, tf.ny, tf.nz
    for i, j, k in ti.ndrange(nx, ny, nz):
        if i == 0 or j == 0 or k == 0 or i == nx - 1 or j == ny - 1 or k == nz - 1:
            v = store[i, j, k] * ratio
            tf.M_am[i, j, k][0, 0] = v
            tf.M_prev_am[i, j, k][0, 0] = v
            tf.M_new_am[i, j, k][0, 0] = v


@ti.kernel
def probe_block_k(tf: ti.template(), iv: ti.i32, jv: ti.i32, kv: ti.i32,  # type: ignore
                  out: ti.template()):  # type: ignore
    """Copy one voxel's spatial 3x3 block to a small field (avoids a full
    to_numpy every probe: the phase instrument samples every few steps)."""
    for a, b in ti.ndrange(3, 3):
        out[a, b] = tf.M_am[iv, jv, kv][a + 1, b + 1]


@ti.kernel
def max_abs_m_k(tf: ti.template(), out: ti.template()):  # type: ignore
    """max |M| over the grid (the bounded-energy guard, launcher pattern).

    NOTE: this CANNOT detect NaN. atomic_max compares, and every comparison
    against NaN is false, so a NaN field still reports a finite max. Use
    `sum_abs_m_k` for the finiteness test (addition propagates NaN)."""
    for i, j, k in tf.M_am:
        m = tf.M_am[i, j, k]
        for a, b in ti.static(ti.ndrange(4, 4)):
            ti.atomic_max(out[None], ti.abs(m[a, b]))


@ti.kernel
def sum_abs_m_k(tf: ti.template(), out: ti.template()):  # type: ignore
    """sum |M| over the grid: the NaN-propagating finiteness test (see above)."""
    for i, j, k in tf.M_am:
        m = tf.M_am[i, j, k]
        out[None] += ti.abs(m).sum()


# ================================================================
# the harness
# ================================================================
class Harness:
    """A driven 4x4 eta-stack run on the certified canonical integrator."""

    def __init__(self, n=N_GRID, h=H_RES, dt=DT, delta=DELTA,
                 boundary="track", project_mixed=True):
        self.n, self.h, self.dt, self.delta = n, h, dt, delta
        self.boundary = boundary            # "track" (B7 guard) | "pin"
        self.project_mixed = project_mixed  # staged 4x4 (task doc SS 2)
        self.tf = TensorField([n * h * 1e-18] * 3, n**3)
        self.t0_00 = ti.field(ti.f32, shape=(n, n, n))
        self._acc = ti.field(ti.f32, shape=())
        self._leak = ti.field(ti.f32, shape=())
        self._probe = ti.field(ti.f32, shape=(3, 3))
        self.t = 0.0
        self.step_i = 0

    # ---- state ----
    def load_fixedj(self, tag):
        info = seeds.load_npz_M(self.tf, fixj_path(tag))
        if not info["ok"]:
            raise SystemExit(f"load failed: {info['warn']}")
        snap_time_entry_k(self.tf, self.t0_00)
        self.t = 0.0
        self.step_i = 0
        self._leak[None] = 0.0
        return info

    def load_vacuum(self):
        """A defect-free covariant vacuum box: diag(-g, 1, delta, 0) everywhere
        (the G-vac gate arena)."""
        m = np.zeros((self.n, self.n, self.n, 4, 4), np.float32)
        m[..., 0, 0] = -G0
        m[..., 1, 1] = 1.0
        m[..., 2, 2] = self.delta
        m[..., 3, 3] = 0.0
        self.tf.M_am.from_numpy(m)
        self.tf.M_prev_am.from_numpy(m)
        self.tf.M_new_am.from_numpy(m)
        snap_time_entry_k(self.tf, self.t0_00)
        self.t = 0.0
        self.step_i = 0
        self._leak[None] = 0.0

    def set_fixed_j(self, om_target):
        """SET-J kick (engine2_pde), then RELEASE: the constraint is NOT
        re-imposed during the driven run (that is the whole point)."""
        return pde.set_fixed_j(self.tf, self.h, RENV, om_target, self.dt, shell=1)

    # ---- ledger reads ----
    def _sum(self, kern, *args):
        self._acc[None] = 0.0
        kern(self.tf, *args, self._acc)
        return float(self._acc[None])

    def v4_total(self, sg):
        return self._sum(v4_sum_k, sg, self.delta, W1)

    def dv4_dsg(self, sg):
        return self._sum(dv4_dsg_sum_k, sg, self.delta, W1)

    def kinetic(self):
        return self._sum(kinetic_sum_k, self.dt)

    def max_abs_m(self):
        self._acc[None] = 0.0
        max_abs_m_k(self.tf, self._acc)
        return float(self._acc[None])

    def carried_j(self):
        return pde.read_carried_j(self.tf, self.h, RENV, self.dt)

    def is_finite(self):
        """NaN-safe health check (max_abs_m cannot see NaN: see sum_abs_m_k)."""
        self._acc[None] = 0.0
        sum_abs_m_k(self.tf, self._acc)
        v = float(self._acc[None])
        return bool(np.isfinite(v)) and v < 1e12

    @property
    def mixed_leak(self):
        return float(self._leak[None])

    # ---- the driven step ----
    def step(self, eps, om_bar, ramp_t=RAMP_T):
        """One canonical leapfrog step at the CURRENT drive sample.

        Sampling convention (blindspot B2): sg is sampled at the MIDPOINT
        t + dt/2 of the step (the leapfrog-consistent choice); the ledger uses
        the same sample, so drive work and drive force always agree."""
        t_mid = self.t + 0.5 * self.dt
        sg = sg_of(t_mid, eps, om_bar, ramp_t)
        tf, h, dt = self.tf, self.h, self.dt

        if self.boundary == "track":
            track_boundary_k(tf, self.t0_00, sg / G0)

        pde.compute_eta_flux(tf, 0, h)
        pde.evolve_M_eta_start(tf, dt, h)
        pde.compute_eta_flux(tf, 1, h)
        pde.evolve_M_eta_finish(tf, dt, h, sg, self.delta, W1)
        tf.swap_matrix_buffers()

        if self.project_mixed:
            project_mixed_k(tf, self._leak)

        self.t += dt
        self.step_i += 1
        return sg

    def drive_power(self, eps, om_bar, ramp_t=RAMP_T):
        """Instantaneous drive power P(t) = dV4/dsg * dsg/dt: the rate at which
        the prescribed background does work on the field. At injection lock its
        TIME AVERAGE is zero (task doc SS 3 / SS 4, the Adler row)."""
        t_mid = self.t - 0.5 * self.dt
        sg = sg_of(t_mid, eps, om_bar, ramp_t)
        return self.dv4_dsg(sg) * dsg_dt(t_mid, eps, om_bar, ramp_t)

    # ---- the phase probe ----
    def probe_axis(self, probe, which=1):
        """An eigen-axis at a probe voxel from the spatial 3x3 block.

        which = 1 -> the MIDDLE eigenvector = the delta clock-hand (the
        M5.23.2 D1-proven probe); which = 2 -> the LEADING eigenvector, which is
        the local clock ROTATION AXIS (the a0_conj construction rotates about
        it, m5_21_9_d_fixedj.a0_conj). Signs are arbitrary (apolar n = -n, the
        spin-1/2 result) — treat mod pi downstream."""
        probe_block_k(self.tf, probe[0], probe[1], probe[2], self._probe)
        m = self._probe.to_numpy().astype(np.float64)
        if not np.isfinite(m).all():
            return None            # caller decides; eigh would raise instead
        m = 0.5 * (m + m.T)
        _, V = np.linalg.eigh(m)
        return V[:, which]


# ================================================================
# the clock-phase instrument (blindspot B9: apolar mod-pi)
# ================================================================
class PhaseTracker:
    """Unwrapped clock phase from an APOLAR axis.

    The axis rotates in a plane; project it on a fixed orthonormal basis
    (e1, e2) of that plane and take phi = atan2(<a,e2>, <a,e1>). Because the
    axis is apolar the raw phi is defined only mod pi, so each increment is
    folded into (-pi/2, pi/2] BEFORE accumulating: a naive 2-pi unwrap would
    double or halve the measured rate (blindspot B9). Validated on the fixed-J
    live hold at known omega* before any verdict consumes it.
    """

    def __init__(self, a0, axis_hint=None):
        a0 = np.asarray(a0, float)
        a0 = a0 / np.linalg.norm(a0)
        if axis_hint is None:
            tmp = np.array([0.0, 0.0, 1.0])
            if abs(float(np.dot(tmp, a0))) > 0.9:
                tmp = np.array([1.0, 0.0, 0.0])
            axis_hint = np.cross(a0, tmp)
        n = np.asarray(axis_hint, float)
        n = n / np.linalg.norm(n)
        e1 = a0 - float(np.dot(a0, n)) * n
        if np.linalg.norm(e1) < 1e-9:
            raise ValueError("probe axis is parallel to the rotation axis")
        e1 = e1 / np.linalg.norm(e1)
        e2 = np.cross(n, e1)
        self.e1, self.e2, self.n = e1, e2, n
        self.phi_raw_prev = self._raw(a0) or 0.0
        self.phi_acc = 0.0
        self.hist = []

    def _raw(self, a):
        a = np.asarray(a, float)
        a = a - float(np.dot(a, self.n)) * self.n
        nn = float(np.linalg.norm(a))
        if nn < 1e-12:
            return None
        a = a / nn
        return math.atan2(float(np.dot(a, self.e2)), float(np.dot(a, self.e1)))

    def update(self, a, t=None):
        if a is None:
            return self.phi_acc
        r = self._raw(a)
        if r is None:
            return self.phi_acc
        d = r - self.phi_raw_prev
        while d > 0.5 * math.pi:       # APOLAR fold
            d -= math.pi
        while d <= -0.5 * math.pi:
            d += math.pi
        self.phi_acc += d
        self.phi_raw_prev = r
        if t is not None:
            self.hist.append((float(t), float(self.phi_acc)))
        return self.phi_acc

    def rate(self, t_from=None):
        """Least-squares dphi/dt over the recorded history (optionally from a
        start time: use it to skip the drive ramp)."""
        if len(self.hist) < 3:
            return float("nan")
        h = np.array(self.hist, float)
        if t_from is not None:
            h = h[h[:, 0] >= t_from]
            if len(h) < 3:
                return float("nan")
        A = np.vstack([h[:, 0], np.ones(len(h))]).T
        return float(np.linalg.lstsq(A, h[:, 1], rcond=None)[0][0])


def adiabatic_vacuum_m00(sg):
    """The uniform-field analytic reference: the exact V4 zero is
    M_vac = diag(-sg, 1, delta, 0), so the time-time entry tracks -sg(t)
    quasi-statically (om_bar << the stiff M00 mode ~ 78 at g = 8)."""
    return -float(sg)
