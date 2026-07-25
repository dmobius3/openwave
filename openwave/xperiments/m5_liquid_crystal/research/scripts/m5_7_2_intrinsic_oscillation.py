"""
M5.7.2 — is the defect ITSELF the metastable particle? (intrinsic-oscillation measurement)

M5.7.1 found that a *seeded* l=1 perturbation disperses (Close's "mostly dispersion"). But
M5.6.2b found the biaxial hedgehog **sources its own twist** (`C_μν≠0` drives a ψ-independent
force — released from rest, the defect oscillates on its own). So the synthesis (the maintainer's
2026-05-28 call): stop seeding from outside — **measure the defect's OWN oscillation** and ask
whether it is the long-lived localized object we're hunting.

Two independent questions, two metrics:
  1. PARTICLE?  Does the **dynamical** energy (kinetic + curvature, excluding the static V well)
     stay gathered at the core, or radiate away?  → `E_core_frac(t) = Σ_{r<R_core}H_dyn / Σ H_dyn`.
  2. CLOCK?     Does the core director oscillate at a **coherent frequency** (a sharp FFT peak)
     or broadband-decay (incoherent dissolve)?  → FFT a near-core director probe; peak prominence
     `P = peak / median(spectrum)`: P≫1 ⇒ a clock, P~1 ⇒ broadband.

Caveat (honest): the *clean* Zitterbewegung clock (`ω=2mc²/ℏ`) is Duda's 4D-Lorentz negative-energy
mechanism (M5.8). In pure 3D the defect's intrinsic oscillation is the `C_μν`-sourced twist — it may
or may not be a single clean frequency. M5.7.2 measures what the 3D dynamics actually do; a coherent
peak is a proto-clock (great), broadband means the clean clock needs 4D (M5.8) and the 3D "particle"
claim rests on localization persistence alone.

Substrate + engine: identical to M5.7.1 (the numpy mirror of production `evolve_M` + V-on `b=0`
amplitude well). Released from rest. V-on vs V-off contrast (does the amplitude well help the
oscillation persist / cohere? — recall V confines amplitude but is rotation-flat, M5.6.5c).

USAGE:
    python -m openwave.xperiments.m5_liquid_crystal.research.sandbox_v7.m5_7_2_intrinsic_oscillation
    M57_N=32 M57_STEPS=3000 python -m ...m5_7_2_intrinsic_oscillation   # fast low-res smoke
"""
import os

import numpy as np

# ----- substrate (matches _topo_biaxial1_von.py + M5.7.1) ---------------------------
DTYPE = np.float32
DELTA = 0.30
S2_STAR = 1.0 + DELTA**2
N = int(os.environ.get("M57_N", 44))
L = 5.0
RC = 0.9                                 # radial core melt = r₀
RHOC = 0.9                               # z-disclination melt
C_WAVE = 0.3
C_COEF = 2.0                             # V-on well stiffness (K=1 analog)
A_COEF = -2.0 * C_COEF * S2_STAR
DT = float(os.environ.get("M57_DT", 0.004))
N_STEPS = int(os.environ.get("M57_STEPS", 4000))   # long: need many periods for a clean FFT
SAMPLE_EVERY = 5                                    # fine time sampling for the spectrum
R_CORE = 3.0 * RC                        # localization aperture
PROBE_R = 2.5 * RC                       # director-probe radius (active shell, off frozen core)


# ----- matrix helpers (numpy mirror of the @ti.func ops; identical to M5.7.1) -------
def matmul(A, B):
    return np.matmul(A, B)


def commf(A, B):
    return np.matmul(A, B) - np.matmul(B, A)


def frob2(A):
    return np.einsum("...ab,...ab->...", A, A)


def tr(A):
    return np.einsum("...aa->...", A)


def central(f, axis, h):
    out = np.zeros_like(f)
    sp, sm, so = [slice(None)] * f.ndim, [slice(None)] * f.ndim, [slice(None)] * f.ndim
    sp[axis], sm[axis], so[axis] = slice(2, None), slice(0, -2), slice(1, -1)
    out[tuple(so)] = (f[tuple(sp)] - f[tuple(sm)]) / (2 * h)
    return out


def V_M(M, a, b, c):
    M2 = matmul(M, M)
    tr2 = tr(M2)
    return a * tr2 - b * tr(matmul(M2, M)) + c * tr2 * tr2


def dV_M(M, a, b, c):
    M2 = matmul(M, M)
    tr2 = tr(M2)[..., None, None]
    return 2.0 * a * M - 3.0 * b * M2 + 4.0 * c * tr2 * M


def force(M, h, a, b, c):
    Mx, My, Mz = central(M, 0, h), central(M, 1, h), central(M, 2, h)
    cc = lambda A, B: commf(commf(A, B), B)
    Gx = 8.0 * (cc(Mx, My) + cc(Mx, Mz))
    Gy = 8.0 * (cc(My, Mx) + cc(My, Mz))
    Gz = 8.0 * (cc(Mz, Mx) + cc(Mz, My))
    divG = central(Gx, 0, h) + central(Gy, 1, h) + central(Gz, 2, h)
    return C_WAVE**2 * divG - dV_M(M, a, b, c)


def energy_parts(M, M_prev, h):
    """Returns (H_dyn, H_total). H_dyn = kinetic + curvature (the MOTION energy, V excluded);
    H_total adds the static V well (for the conservation check)."""
    Mdot = (M - M_prev) / DT
    kinetic = 0.5 * frob2(Mdot)
    Mx, My, Mz = central(M, 0, h), central(M, 1, h), central(M, 2, h)
    curv = 4.0 * (frob2(commf(Mx, My)) + frob2(commf(Mx, Mz)) + frob2(commf(My, Mz)))
    H_dyn = kinetic + C_WAVE**2 * curv
    return H_dyn, H_dyn + V_M(M, A_COEF, 0.0, C_COEF)


# ----- biaxial hedgehog background (mirror of seed_biaxial_hedgehog_M) ----------------
def build_biaxial_M():
    xs = np.linspace(-L, L, N)
    h = xs[1] - xs[0]
    X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
    r = np.sqrt(X**2 + Y**2 + Z**2) + 1e-9
    rho = np.sqrt(X**2 + Y**2) + 1e-9
    rhat = np.stack([X / r, Y / r, Z / r], -1)
    azim = np.stack([-Y / rho, X / rho, np.zeros_like(r)], -1)
    sdisc = np.clip(rho / RHOC, 0, 1)
    shrink = (sdisc**2 * (3 - 2 * sdisc))[..., None]
    ephi = azim * shrink
    ephi = ephi - np.einsum("...a,...a->...", ephi, rhat)[..., None] * rhat
    etheta = np.cross(ephi, rhat)
    srad = (r / np.sqrt(r**2 + RC**2))[..., None, None]
    diso = (1.0 + DELTA) / 3.0
    d0 = diso + srad[..., 0, 0] * (1.0 - diso)
    d1 = diso + srad[..., 0, 0] * (DELTA - diso)
    d2 = diso + srad[..., 0, 0] * (0.0 - diso)
    oprod = lambda v: v[..., :, None] * v[..., None, :]
    M = (d0[..., None, None] * oprod(rhat) + d1[..., None, None] * oprod(etheta)
         + d2[..., None, None] * oprod(ephi))
    interior = np.zeros(r.shape, bool)
    interior[2:-2, 2:-2, 2:-2] = True
    active = (r > 2 * RC) & (rho > RHOC) & interior
    in_core = active & (r < R_CORE)
    floor = in_core.sum() / active.sum()                 # uniform-dispersion floor for E_core_frac
    # probe voxel: nearest active voxel to (PROBE_R, 0, 0) on the x-axis (off the z-disclination)
    ix = int(round((PROBE_R + L) / (2 * L) * (N - 1)))
    ic = (N - 1) // 2
    probe = (ix, ic, ic)
    return dict(M=M.astype(DTYPE), h=DTYPE(h), r=r.astype(DTYPE), active=active, probe=probe, floor=floor)


def principal_director(m3):
    """Director = eigenvector of the largest eigenvalue of the symmetric 3×3 m3 (numpy eigh)."""
    w, v = np.linalg.eigh(m3.astype(np.float64))
    return v[:, np.argmax(w)]


# ----- evolve + record intrinsic-oscillation diagnostics ----------------------------
def evolve_intrinsic(bg, v_on, n_steps):
    """Evolve the unperturbed biaxial hedgehog, released from rest. Record per sample:
    E_core_frac (dynamical-energy localization), H_total (conservation), and the core
    director probe n̂(t) (3 comps) for the spectral/coherence analysis."""
    h, active, r, probe = bg["h"], bg["active"], bg["r"], bg["probe"]
    mask = active[..., None, None]
    in_core = active & (r < R_CORE)
    a, b, c = (A_COEF, 0.0, C_COEF) if v_on else (0.0, 0.0, 0.0)
    M, M_prev = bg["M"].copy(), bg["M"].copy()

    Ecore, Htot, ndir, tsteps = [], [], [], []

    def sample():
        H_dyn, H_total = energy_parts(M, M_prev, h)
        H_dyn = (H_dyn.astype(np.float64)) * active
        tot = H_dyn.sum() + 1e-30
        Ecore.append(H_dyn[in_core].sum() / tot)
        Htot.append(float((H_total.astype(np.float64) * active).sum()))
        ndir.append(principal_director(M[probe]))
        tsteps.append(len(tsteps) * SAMPLE_EVERY)

    sample()
    for step in range(1, n_steps + 1):
        f = force(M, h, a, b, c)
        M_new = np.where(mask, 2 * M - M_prev + DT**2 * f, bg["M"])
        M_prev, M = M, M_new
        if step % SAMPLE_EVERY == 0:
            sample()
            if step % 800 == 0:
                print(f"      step {step:>5}/{n_steps}  E_core={Ecore[-1]:.3f}", flush=True)
    return dict(Ecore=np.array(Ecore), Htot=np.array(Htot),
                ndir=np.array(ndir), finite=np.isfinite(M).all())


# ----- spectral coherence: FFT the probe director, robust concentration metrics --------
def spectral_coherence(ndir):
    """Hann-windowed FFT of each director component (DC removed). Returns the component whose
    spectrum is most CONCENTRATED, by two robust measures (the old peak/median blows up because
    most bins are ~0 float-noise — median≈0 — so it is replaced):
      • band_frac = fraction of AC power within ±2 bins of the peak (a clean tone → near 1)
      • prom = peak / MEAN(spectrum) (bounded ~n_bins for a pure tone, ~few for broadband)."""
    dt_s = DT * SAMPLE_EVERY
    n = ndir.shape[0]
    win = np.hanning(n)
    best = dict(band_frac=0.0, prom=0.0, f=0.0, peak_frac=0.0, comp=-1)
    freqs = np.fft.rfftfreq(n, d=dt_s)
    for comp in range(3):
        sig = ndir[:, comp] - ndir[:, comp].mean()
        if np.std(sig) < 1e-9:
            continue
        spec = np.abs(np.fft.rfft(sig * win)) ** 2
        spec[0] = 0.0                                    # drop DC
        total = spec.sum()
        if total < 1e-30:
            continue
        k = int(np.argmax(spec))
        lo, hi = max(0, k - 2), min(len(spec), k + 3)
        band_frac = spec[lo:hi].sum() / total            # power within ±2 bins of the peak
        if band_frac > best["band_frac"]:
            best = dict(band_frac=float(band_frac), prom=float(spec[k] / (spec.mean() + 1e-30)),
                        f=float(freqs[k]), peak_frac=float(spec[k] / total), comp=comp)
    return best


# ====================================================================================
def run(bg, v_on, label):
    res = evolve_intrinsic(bg, v_on, N_STEPS)
    E = res["Ecore"]; H = res["Htot"]; floor = bg["floor"]
    Hdrift = (H.max() - H.min()) / (abs(H[0]) + 1e-30)
    sp = spectral_coherence(res["ndir"])
    # localization vs the uniform-dispersion floor: (E_end − floor)/(E₀ − floor) — 1=fully localized, 0=dispersed
    E_loc = (E[-1] - floor) / (E[0] - floor + 1e-30)
    print(f"  [{label}]  finite={res['finite']}  H-drift={100*Hdrift:.2f}%")
    print(f"      E_core(dyn): {E[0]:.3f} → min {E.min():.3f} → end {E[-1]:.3f}   "
          f"(floor {floor:.3f}; localization-excess {100*E_loc:.0f}%)")
    print(f"      core osc: f={sp['f']:.3f}/t  band±2-frac={100*sp['band_frac']:.0f}%  "
          f"peak-bin={100*sp['peak_frac']:.0f}%  prom(peak/mean)={sp['prom']:.0f}×  (comp {sp['comp']})")
    return dict(label=label, E=E, E_loc=E_loc, Hdrift=Hdrift, sp=sp, finite=res["finite"])


def main():
    print("=" * 80)
    print("M5.7.2 — intrinsic oscillation: is the defect ITSELF the metastable particle?")
    print(f"  grid {N}³  L={L}  c={C_WAVE}  dt={DT}  steps={N_STEPS}  δ={DELTA}  R_core={R_CORE:.2f}  f32")
    print("=" * 80)
    bg = build_biaxial_M()
    print(f"  active {100*bg['active'].mean():.0f}%   E_core floor={bg['floor']:.3f}   "
          f"probe {bg['probe']} (r≈{PROBE_R:.2f}, x-axis)", flush=True)

    print("\n  V-ON (b=0 amplitude well) — does confinement help the defect persist/cohere?")
    von = run(bg, True, "V-on")
    print("\n  V-OFF (control — M5.5.4 'dilutes') —")
    voff = run(bg, False, "V-off")

    print("\n" + "=" * 80)
    pipeline_ok = von["finite"] and voff["finite"] and von["Hdrift"] < 0.05
    # PARTICLE: dynamical energy stays above the uniform floor — 3-way (YES ≥0.6, PARTIAL 0.3–0.6, NO <0.3)
    el = von["E_loc"]
    particle = el >= 0.6
    p_label = "YES" if el >= 0.6 else ("PARTIAL" if el >= 0.3 else "NO")
    # CLOCK: a genuinely concentrated spectrum — most AC power in a narrow band around one frequency
    clock = von["sp"]["band_frac"] >= 0.5 and von["sp"]["peak_frac"] >= 0.3
    print(f"M5.7.2  (N={N}, {N_STEPS} steps, dt={DT}):  pipeline={pipeline_ok}  "
          f"(V-on H-drift {100*von['Hdrift']:.2f}%)")
    print(f"  PARTICLE? (dyn energy stays localized): {p_label} "
          f"— V-on localization-excess {100*von['E_loc']:.0f}% (vs V-off {100*voff['E_loc']:.0f}%, "
          f"H-drift {100*voff['Hdrift']:.1f}% {'untrustworthy' if voff['Hdrift'] >= 0.05 else 'ok'})")
    print(f"  CLOCK? (coherent core oscillation): {'YES' if clock else 'NOT CRISP'} "
          f"— band±2 {100*von['sp']['band_frac']:.0f}%, peak-bin {100*von['sp']['peak_frac']:.0f}% "
          f"at f={von['sp']['f']:.3f}/t")
    if particle and clock:
        print("  → THE DEFECT IS A LONG-LIVED COHERENT OSCILLATOR — the metastable 'particle' is the")
        print("    defect's own intrinsic clock (proto-Zitterbewegung; the clean ω=2mc²/ℏ is M5.8/4D).")
    elif particle:
        print("  → LOCALIZED but not a single-frequency clock in 3D — energy stays gathered (a particle);")
        print("    the coherent clock needs the 4D negative-energy mechanism (M5.8). Still the metastable seed.")
    elif el >= 0.3 and clock:
        print(f"  → PARTIALLY localized ({100*el:.0f}% excess, {von['E'][-1]:.2f} vs floor {bg['floor']:.2f} = "
              f"{von['E'][-1]/bg['floor']:.1f}× uniform) + a coherent core oscillation (f={von['sp']['f']:.3f}/t).")
        print("    A metastable oscillator that sheds some energy then settles — promising; needs the biaxiality-")
        print("    stabilizing V term (Q7) / stronger confinement to tighten, and 4D (M5.8) for the clean clock.")
    else:
        print("  → The defect DELOCALIZES under Evolve-PDE — V confines amplitude but not orientation")
        print("    (M5.6.5c); a biaxiality-stabilizing V term (Q7) or stronger confinement is needed.")
    print("PASS" if pipeline_ok else "FAIL — pipeline/conservation issue")
    print("=" * 80)
    return 0 if pipeline_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
