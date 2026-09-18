# METHOD (written before any computation; not edited afterwards — later changes are appended as dated notes)

## Conventions I fix
- Quaternion -> SU(2): i -> -i sigma_x, j -> -i sigma_y, k -> -i sigma_z, i.e.
  q = w+xi+yj+zk -> [[w - i z, -y - i x], [y - i x, w + i z]] = [[a, b], [-conj b, conj a]].
  Checked to be a homomorphism (i j = k, etc.) numerically and exactly.
- D^j(g): the representation Sym^{2j}(C^2) on binary forms, basis v_m = x^{j+m} y^{j-m} / sqrt((j+m)!(j-m)!/(2j)!),
  with g.x = g11 x + g21 y, g.y = g12 x + g22 y. Then J+ v_m = sqrt((j-m)(j+m+1)) v_{m+1} (Condon-Shortley phases),
  D^{1/2}(g) = g. I verify: homomorphism, unitarity, conj(D_{mk}) = (-1)^{m-k} D_{-m,-k} (integer j).
- Clebsch-Gordan: my own implementation of the Racah closed formula, exact (each coefficient is +-sqrt(rational)).
  Cross-check by a second construction: highest-weight vector of V_J inside V_{j1} (x) V_{j2} as the kernel of J+
  in the weight-J subspace (exact rational linear algebra), phase fixed by <j1 j1; j2 (J-j1) | J J> > 0, then repeated
  application of J- and normalization. Compare all coefficients for all spins used (j1,j2 <= 9 as needed).
- Exact arithmetic: a small multi-quadratic number class MQ = sum_s c_s sqrt(s), s squarefree positive integer,
  c_s in Q(i) (pairs of Fractions). Products sqrt(s)sqrt(t) = gcd(s,t) sqrt(st/gcd^2). Zero test is exact because
  square roots of distinct squarefree integers are linearly independent over Q(i). Inverse by successive
  Galois conjugation over each prime radical, then over i. All values reported "exactly" come from this class
  (or from rational/Q(sqrt5) arithmetic); nothing exact is identified from floats unless stated.

## Group (item 0, 1)
- Enumerate Gamma by closure of {q1, q2} under multiplication with exact Q(sqrt5) quaternions. Report |Gamma|.
- Derived subgroup: closure of all commutators; compare orders.
- Conjugacy classes by exact conjugation.
- Isotypic projectors of V_3 without tables: class sum Z of the class of q1 in D^3; verify Z^2 = c Z exactly, so
  P = Z/c and I - P are the two isotypic projectors; ranks (traces) identify the 3- and 4-dim sectors; verify
  P^2 = P, P commutes with D^3(q1), D^3(q2); irreducibility via character norm (1/|G|) sum |tr(P D^3(h))|^2 = 1.
- Characters: chi_sigma(h) = tr(P_sigma D^3(h)) exactly; chi_j(h) = U_{2j}(Re h) (Chebyshev 2nd kind) exactly.
  dim Hom(sigma, V_j) = (1/|G|) sum conj(chi_sigma) chi_j, for j = 0..9 (n = 2j <= 18), plus half-integer j
  (should vanish since -1 acts by (-1)^{2j} on V_j and trivially on sigma; verified).
  Declaration planned: DERIVED (from generators, by the above computation). I recognize the group by its
  order and generators, but no table is used.
- Explicit intertwiners (numerical, mpmath 50 digits, and float64 for quadrature): eta = orthonormal basis of
  range P_sigma (so sigma(h) := eta^dag D^3(h) eta is the one fixed matrix representation); at level j, a basis of
  Hom(sigma,V_j) by group averaging E = (1/|G|) sum D^j(h) X D^3(h)^dag applied to eta, Gram-Schmidt with
  <eta,eta'> = tr(eta^dag eta')/d; check D^j(h)eta = eta sigma(h) at h = q1, q2, eta^dag eta = I, eta^dag eta' = 0;
  report max residuals (required < 1e-40 at 50 digits).

## The "P-trick" (exact route)
A section Phi(g) = u^T D^3(g) eta is replaced by Phi~ = Phi eta^dag = u^T D^3(g) P (values in C^7). Right
multiplication by the isometry eta^dag preserves pointwise norms and commutes with N, DN, Delta, level projections
and inner products, so every quantity asked is computed from P alone (no square roots from normalizing eta).
Level-J components are kept as sums of rank-one terms x^T D^J(g) Y (x in C^{2J+1}, Y (2J+1)x7). Rules:
- conj(x^T D y) = (Theta x)^T D (Theta y), (Theta x)_m = (-1)^m conj(x_{-m});
- (x1^T D y1)(x2^T D y2) = sum_L [x1 (x) x2]_L^T D^L [y1 (x) y2]_L;
- trilinear T(f1,f2,f3)_a = sum_b conj(f1_b) f2_b f3_a: level-J part = sum_L A_L^T D^J Z_L with
  A_L = [[Theta x1 (x) x2]_L (x) x3]_J, W_L = sum_b [Theta Y1_b (x) Y2_b]_L, Z_L[:,a] = [W_L (x) Y3_a]_J;
- <x^T D^J Y, x'^T D^J' Y'> = delta_JJ' (x^dag x') tr(Y^dag Y')/(2J+1).
- A level-3 (block) function has Y proportional to P; its fibre vector is sum_r x_r tr(P Y_r)/d.
N(Phi) = T(Phi,Phi,Phi); DN_Phi[h] = T(h,Phi,Phi)+T(Phi,h,Phi)+T(Phi,Phi,h).
Normalized Phi: fibre u with ||u||^2 = 7/d.

## Items
- 2: F(u) := fibre of Pi_6 N(Phi_u); test F(u) - (u^dag F / u^dag u) u == 0 exactly; Q = <Phi,N(Phi)> = int |Phi|^4.
  r6hat(u) = ||[u (x) Theta u]_6||^2/||u||^4 exactly; stationarity on the unit sphere: exact real gradient
  (first derivative along all 14 real directions v_m, i v_m after projecting out u), all zero or not.
- 3: stabilizers of the line [u] under the action r.Phi(g) = Phi(r^-1 g), i.e. u -> conj(D(r)) u
  (= D(r)u-line for these real u up to conjugation): continuous part from the Lie algebra (exact kernel of
  X -> X.u mod i u), finite part from the Majorana roots of the binary sextic of u (rotations permuting roots,
  enumerated by matching root pairs, then verified exactly/numerically). Character: eigenvalue on Phi. Dimension of
  the block subspace with the same character: character projection by averaging over the (finite part of the)
  stabilizer, plus weight analysis for the torus parts.
- 4: xi = -g sum_{n != 6} Pi_n N(Phi)/(n(n+2)-48); ||Pi_n xi||^2/g^2 exact per level from the rank-one pairs.
  Character check: apply the stabilizer generators to xi's fibre data (numerically and via the equivariance argument).
- 5: Pi_6 DN_Phi[xi]/g as a fibre vector (g = 1); component along Phi = <Phi,.>, perp norm; U5: <e_t,.>_R,
  <i e_t,.>_R; U6: <tau_x,.>_R, <tau_y,.>_R. Also compare <Phi,DN xi> with -3 sum ||Pi_n N||^2/mu_n (identity).
- 6: L := real-linear operator E -> Pi_perp Pi_6 DN_Phi[E] - Q E on the tangent space (fibre coordinates, exact
  14x14 real matrix); quadratic form values on e_t, i e_t, tau_x, tau_y, cross term. r6hat second derivative by
  the exact expansion of ||rho_6(cos s u + sin s w)||^2.
- 7: solve L kappa = -[Pi_6 DN xi]_perp exactly on the tangent space; ker L computed exactly and identified
  (compare with the span of infinitesimal rotations X_a.Phi projected perp to Phi); kappa taken orthogonal (real
  inner product) to ker L; report residual component (projection of the RHS onto ker L) = the equation's perp
  component with kappa; without kappa = item 5 perp norm.
- 8: lambda_4 = g<Phi, DN_Phi[xi+kappa]> - gQ<Phi,kappa>, both with and without kappa; exact reason for equality:
  Re/Im parts of <Phi, DN_Phi kappa> reduce to 3Re<N(Phi),kappa> and Im<N(Phi),kappa>, which vanish because
  Pi_6 N(Phi) = QPhi and <Phi,kappa> = 0. Sign: lambda_4 = -3 g^2 sum_{n != 6} ||Pi_n N(Phi)||^2/(n(n+2)-48),
  so the sign is fixed by the weighted level sum (levels below 6 contribute with + sign, above with -); give the
  argument with the exact values.
- 9: tabulate every zero; arguments: (a) left rotations (commute with Delta, with N, preserve the Gamma-sector
  because left and right translations commute) and the stabilizer character: when the character subspace of the
  block is C.Phi, all block quantities covariant with that character are multiples of Phi, so perp parts vanish;
  (b) an antiunitary symmetry: conj composed with an intertwiner W in Hom(sigma, conj sigma) (exists since
  chi_sigma is real; verified) and a pi-rotation about the y-axis, acting on fibres as u -> +-conj(u); preserves
  Delta, N, lambda real, g real, and the sector; forces the i e_t / tau_y components of real-covariant quantities
  to vanish; (c) algebraic identities (U(1) phase, self-adjointness of Delta). Any zero without such argument is
  reported unresolved.
- 11: the same pipeline at sin^2 t = 1/4.
- 12 (numerical route, independent of the exact route): float64 Haar quadrature on SU(2): Euler angles,
  40 x 40 uniform in alpha, gamma and 20-node Gauss-Legendre in cos beta, exact for integrands with spin <= 18
  (after alpha,gamma integration only D^J_00 = P_J(cos beta) survive). Evaluate Phi and N(Phi) pointwise with
  explicit eta; project onto every level n <= 18 (and n = 20, 22 as a check) keeping every intertwiner copy; build
  -Delta xi by 4 * Casimir (Jx^2+Jy^2+Jz^2 from explicit spin matrices, acting on the left index as C^T) applied to
  item 4's xi (converted from the exact route); residual = max |(-Delta-48)xi + g(N - Pi_6 N)| over all
  coefficients, required < 1e-10. Why no component above 18: N(Phi) is a sum of products of three spin-3 matrix
  coefficients, whose CG series stops at spin 9; confirmed by Parseval (quadrature int|N|^2 = sum of level norms).

## Item 10 argument (outline, committed now)
Setting: real Hilbert space of L^2 sections in the sector; the equation is the Euler-Lagrange equation of
E(psi) = 1/2<psi,-Delta psi> + g/4 int|psi|^4 with multiplier lambda on ||psi||^2; symmetry group
G = SU(2)_left x U(1)_phase (plus the antiunitary symmetry), all verified to preserve equation and sector.
Step 1 (Lyapunov-Schmidt with blow-up): write psi = a(Phi + a^2 w), Phi in the unit sphere of the block, w
perp block; divide by a^3; the range equation (-Delta-48)w = -g(1-Pi_6)N(Phi + a^2 w) + (lambda-48)a^{-2}... is
solved by the IFT near a = 0 since -Delta-48 is invertible off the block (spectral gap: nearest levels 4 and 8),
N is smooth (cubic, H^s algebra for s > 3/2). This handles the degeneracy at a = 0 (the trivial branch psi = 0
and the 14-real-dimensional kernel): after blow-up the problem at a = 0 is regular off the block.
Step 2 (symmetry reduction): let H be the stabilizer of [Phi] with character chi; H-hat = {(r, chi(r)^-1)} in G.
Restrict to Fix(H-hat), a closed subspace preserved by the equation (principle of symmetric criticality / equivariance).
If Fix(H-hat) meets the block in C.Phi only (item 3 dimension 1), the reduced bifurcation equation is a single real
equation for lambda, solved by lambda = 48 + gQ a^2 + ... with Q = int|Phi|^4 > 0; IFT gives existence, real-
analyticity in a, psi odd / lambda even in a, and uniqueness within Fix(H-hat) near a Phi (for a > 0, modulo
the residual phase, fixed by a = <Phi,psi> real positive).
If the block part of Fix(H-hat) is 2-dimensional (U5, U6): the residual symmetry reduces the block sphere to a
1-real-parameter curve (t for U5; for U6 the curve u(z) with the complex parameter z and the residual symmetry
to be determined); the reduced equation is the derivative of the reduced functional
f_a = a^4 W4 + a^6 W6 + O(a^8) along the curve. Existence near the given point requires: dW4 = 0 there
(item 2), and either d^2W4 != 0 along the free direction(s) (item 6: nondegenerate -> IFT, unique branch) or,
if d^2W4 = 0 identically along the curve, dW6 = 0 at the point (item 5 component along e_t / tau) and
d^2W6 != 0 (needs a further computation) -- otherwise no solution with that leading term.
Global uniqueness modulo G (outside Fix) needs ker L = tangent of the G-orbit (item 7 kernel identification); I will
state which variant is proved for each U. Hypotheses to be listed with item numbers where verified; for each step,
what fails if omitted (e.g. without the spectral gap the range equation is not solvable; without Q != 0 the
lambda equation degenerates; without nondegeneracy mod symmetry, persistence of the particular Phi fails).

---
## Appended note, 2026-09-18 (after computing; the text above is unchanged)
1. Spectral gap: the "nearest levels 4 and 8" guess above was wrong. Item 1 shows that neither sector has a level
   below 6. The nearest levels are n = 10 (sector 3) and n = 8 (sector 4). So every n(n+2)-48 on the xi levels is
   positive. This is what fixes the sign of lambda_4.
2. Item 10 route: the exact kernel computation (items 6/7) shows ker L = tangent space of the SU(2) orbit at all
   six points, in both sectors. So U5 and U6 are nondegenerate modulo symmetry in the full block, not only along
   their curves. The degenerate branch of the outline (d^2 W4 = 0) never occurs, and W6 is not needed. The final
   argument uses Lyapunov-Schmidt plus an equivariant slice (Palais symmetric criticality) for all six points. For
   U1-U4 the Fix(H-hat) route of the outline gives the same result; it is kept as an alternative proof.
3. New route for the item-2 zeros at U5/U6: I proved the exact identity int|Phi|^4 = 1 + beta_sigma r6hat on the
   unit block sphere, with beta_3 = 28/39 and beta_4 = 21/52. The proof is a rank-4 evaluation argument on the
   4-dimensional space of SU(2)-invariant Hermitian quartics; the script is s_verify.py. With this identity, the
   stationarity of int|Phi|^4 reduces to the explicit one-variable r6hat curves.
4. Level-16 zeros at U4 (and U2, U3): these are not explained by stabilizer characters alone. The explanation is
   that Theta u = +-u there and V_8 does not occur in Sym^3 V_3.
5. Stabilizers: the finite parts were enumerated in float64 from the Majorana roots, and the characters were checked
   numerically, with residuals around 1e-15. They are not certified exactly. The exact route does not depend on
   them; they enter only the arguments.
6. The sign of an exact MQ value is decided by 60-digit evaluation, after the exact zero test has shown the value
   is nonzero.
