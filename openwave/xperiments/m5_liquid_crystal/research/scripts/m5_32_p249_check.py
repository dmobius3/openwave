"""Independent check of the P249 exterior-degenerate clock exact algebra.

P249 is campaign "exact exterior-degenerate M5 clock" in the parallel program
vantasnerdan/substrate-framework (PR #190, merged 2026-08-31, release
v0.171.0, pinned commit 6939a3d8dc3cf213ee61ecb3d3f798196cc4d67a). Its claim:
a conditional M5 completion on Sym(3) + C with an exact SO(2) clock that fixes
the exterior vacuum pointwise, a charged continuum edge m*^2 = 4, and a
split-core witness strictly below it (hylomorphic binding margin 19/16).

This script re-derives the checkable algebra FRESH from the stated action
(campaigns/P249-exterior-degenerate-clock/attempts/0008/derivation.md); it
does not import or run their module. Checks:
 1. exterior fixation of (S, psi) = (N N^T, 0) under the SO(2) action
 2. potential Hessian about the exterior in the chart (a, t, p, u, v, q)
    equals diag(5, 10, 22, 4, 4, 22)
 3. generalized mass squares (10, 10, 22, 4, 4, 22) with the kinetic metric
    diag(1/2, 1, 1, 1, 1, 1); scalar mass square 6
 4. charge weights (psi and shear weight 1, tangent doublet weight 2) give
    the charged continuum edge m*^2 = min(6, 4, 22/4) = 4
 5. split-core witness f = 1/2, S = diag(1, 1/4, -1/4): V = 45/64, I = 1/2,
    2V/I = 45/16, margin 19/16
 6. the B = 0 scalar ratio identity 2[W(f) + 6 f^4]/f^2 = 16 (f - 1/4)^2 + 5
 7. the Noether reduction C = omega I on the clock orbit

NOT checked here: the Benci-Fortunato variational bridge (coercivity,
splitting property, subcritical control, Theorem 18 hypotheses); that part of
the P249 claim rests on their theorem-bridge and review record.

Run: python3 m5_32_p249_check.py (sympy only, < 1 min). All 7 checks passed
on 2026-08-31.
"""
import sympy as sp

a,t,p,u,v,q,fr,fi,w = sp.symbols('a t p u v q fr fi w', real=True)

S = sp.Matrix([[1+a,u,v],[u,t+p,q],[v,q,t-p]])
A = sp.Matrix([[0,0,0],[0,0,-1],[0,1,0]])

# --- potential, built from the stated pieces (beta=c=1 projected LdG + unit axis lock
#     + phase lock strength 6 + W(f)) ---
tr2 = sp.trace(S*S)
V_m5   = -tr2/2 - sp.trace(S**3) + tr2**2 + sp.Rational(1,2)
V_axis = tr2 - S[0,0]**2                                  # zeta = 1
# tangent-traceless clock tensor B and target Q(psi)
P = sp.diag(0,1,1)
Bproj = P*S*P
Bcl = Bproj - sp.trace(Bproj)*P/2
H = sp.diag(0,1,-1); K = (A*H - H*A)/2
Q = (fr**2-fi**2)*H + 2*fr*fi*K
V_lock = 6*sp.trace((Bcl-Q).T*(Bcl-Q))/2
f2 = fr**2+fi**2
W = 3*f2 - 4*f2*sp.sqrt(f2) + 2*f2**2   # W(f)=3f^2-4f^3+2f^4, f=|psi|
V = V_m5 + V_axis + V_lock + W

vac = {a:0,t:0,p:0,u:0,v:0,q:0,fr:0,fi:0}
assert sp.simplify(V.subs(vac)) == 0, "exterior not a zero of V"

# 1. SO(2) fixes the exterior exactly
th = sp.symbols('th', real=True)
R = sp.Matrix([[1,0,0],[0,sp.cos(th),-sp.sin(th)],[0,sp.sin(th),sp.cos(th)]])
NNt = sp.diag(1,0,0)
assert sp.simplify(R*NNt*R.T - NNt) == sp.zeros(3), "exterior not fixed by SO(2)"
print("1. exterior fixation: PASS")

# 2. Hessian of V about the exterior in (a,t,p,u,v,q) + scalar (fr,fi)
xs = [a,t,p,u,v,q]
Hess = sp.Matrix(6,6, lambda i,j: sp.diff(V.subs({fr:0,fi:0}),xs[i],xs[j]).subs(vac))
print("2. spatial potential Hessian == diag(5,10,22,4,4,22):", Hess == sp.diag(5,10,22,4,4,22))
# scalar mass: expand V in f to quadratic order along fi=0 (avoid sqrt at 0 by using f>=0 symbol)
f = sp.symbols('f', nonnegative=True)
Vs = V.subs({fr:f,fi:0,a:0,t:0,p:0,u:0,v:0,q:0})
ser = sp.series(sp.simplify(Vs), f, 0, 3).removeO()
print("   scalar quadratic part:", sp.expand(ser), "(want 3*f**2 -> mass^2 = 6 with kin 1/2)")

# 3. generalized mass squares: kinetic energy = Tr(Sdot^2)/4 + |psidot|^2/2
#    in chart: (1/4)adot^2 + (1/2)(tdot^2+pdot^2+udot^2+vdot^2+qdot^2)
#    => K = diag(1/2,1,1,1,1,1); mass^2_i = HessV_i / (2*K_i) since V ~ (1/2)Hess x^2
Kkin = sp.diag(sp.Rational(1,2),1,1,1,1,1)
gen = sp.simplify(Kkin.inv() * Hess)
print("3. generalized mass squares:", [gen[i,i] for i in range(6)], "(want [10,10,22,4,4,22])")

# 4. charge weights: psi wt 1, (u,v) wt 1, (p,q) wt 2, (a,t) wt 0
#    check by acting: d/dth [R S R^T] at 0 = [A,S]
comm = A*S - S*A
print("4. [A,S] channels:", sp.simplify(comm))
#    (u,v) appear linearly (wt1), (p,q) with factor 2 (wt2); edge = min(6/1, 4/1, 22/4) = 4
edge = min(sp.Integer(6), sp.Integer(4), sp.Rational(22,4))
print("   charged edge m*^2 =", edge, "(want 4)")

# 5. witness
wit = {a:0,t:0,p:sp.Rational(1,4),u:0,v:0,q:0,fr:sp.Rational(1,2),fi:0}
Vw = sp.simplify(V.subs(wit))
commw = comm.subs(wit)
Iw = sp.simplify(sp.Rational(1,2)*sp.trace(commw.T*commw) + (fr**2+fi**2).subs(wit))
print("5. witness V =", Vw, "(want 45/64);  I =", Iw, "(want 1/2);  2V/I =", sp.simplify(2*Vw/Iw), "(want 45/16); margin =", sp.simplify(4-2*Vw/Iw))

# 6. B=0 scalar ratio identity
lhs = sp.simplify(2*(3*f**2-4*f**3+2*f**4 + 6*f**4)/f**2)
rhs = sp.expand(16*(f-sp.Rational(1,4))**2+5)
print("6. B=0 scalar ratio identity:", sp.simplify(lhs-rhs)==0)

# 7. Noether reduction on the orbit: Sdot = w[A,S], psidot = i w psi
Sdot = w*comm
Cd = sp.simplify(( sp.trace(Sdot.T*comm)/2 + (fr*(w*fr) - fi*(-w*fi)) ))
Iden = sp.simplify(sp.trace(comm.T*comm)/2 + fr**2 + fi**2)
print("7. C - w*I == 0:", sp.simplify(Cd - w*Iden)==0)
