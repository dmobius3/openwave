"""Resolve the two simplify() failures in s8: they are simplifier failures or
real defects, and the difference matters."""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import sympy as sp
import core as C

HERE = pathlib.Path(__file__).parent
out = []
P = out.append

# --- (A) the span{v3,v-3} formula, with alpha,beta split into real parts ----
p, q, r, s = sp.symbols('p q r s', real=True)
al = p + sp.I * q
be = r + sp.I * s
u = [al, 0, 0, 0, 0, 0, be]
val = sp.simplify(C.rhat6_exact(u))
A2 = p ** 2 + q ** 2
B2 = r ** 2 + s ** 2
tgt = (2 * A2 * B2 + (A2 + B2) ** 2 / 924) / (A2 + B2) ** 2
d = sp.simplify(sp.together(sp.expand(sp.simplify(val - tgt))))
P("span{v3,v-3} family, alpha = p+iq, beta = r+is:")
P(f"  rhat6 - [2|a|^2|b|^2 + (|a|^2+|b|^2)^2/924]/(|a|^2+|b|^2)^2 = {d}")
P(f"  identically zero ? {sp.simplify(d) == 0}")
P(f"  numeric spot: p,q,r,s = 1,2,-3,0.5 -> diff = "
  f"{sp.N(d.subs({p:1,q:2,r:-3,s:sp.Rational(1,2)}),25)}")
P("  => for a UNIT vector in that plane, rhat6 = 2|alpha|^2|beta|^2 + 1/924,")
P("     sweeping [1/924, 463/924] as |alpha||beta| goes from 0 to 1/2.")

# --- (B) the vector-operator relation for f_y under a z-rotation -----------
P("")
P("D(z,t)^dag f_i D(z,t) = R_ij f_j  for a z-rotation, residuals printed in full:")
t = sp.symbols('t', real=True)
Dz = sp.diag(*[sp.exp(-sp.I * m * t) for m in C.MS])
R = sp.Matrix([[sp.cos(t), -sp.sin(t), 0], [sp.sin(t), sp.cos(t), 0], [0, 0, 1]])
for i in range(3):
    lhs = Dz.conjugate().T * C.JOPS_EXACT[i] * Dz
    rhs = R[i, 0] * C.JOPS_EXACT[0] + R[i, 1] * C.JOPS_EXACT[1] + R[i, 2] * C.JOPS_EXACT[2]
    D = sp.Matrix(7, 7, lambda a, b: sp.simplify(sp.expand(
        sp.rewrite(lhs[a, b] - rhs[a, b], sp.exp).rewrite(sp.cos))) if False else
        sp.simplify(sp.expand((lhs[a, b] - rhs[a, b]).rewrite(sp.cos))))
    P(f"  i={i}: max |residual entry| symbolic = {set(map(str, set(D)))}")
    P(f"        zero matrix ? {D == sp.zeros(7,7)}")

txt = "\n".join(map(str, out))
print(txt)
(HERE / "out_s9.txt").write_text(txt + "\n")
