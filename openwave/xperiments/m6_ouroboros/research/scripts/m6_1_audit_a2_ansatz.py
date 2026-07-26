"""AUDIT A2: is A = rhat x grad(phi(r)) cos(wt) identically zero for radial phi?"""
import sympy as sp

x, y, z = sp.symbols("x y z", real=True)
r = sp.sqrt(x ** 2 + y ** 2 + z ** 2)
phi = sp.Function("phi")

grad = [sp.diff(phi(r), v) for v in (x, y, z)]
rhat = [x / r, y / r, z / r]
cross = [rhat[1] * grad[2] - rhat[2] * grad[1],
         rhat[2] * grad[0] - rhat[0] * grad[2],
         rhat[0] * grad[1] - rhat[1] * grad[0]]
print("rhat x grad(phi(r)) components, simplified:", [sp.simplify(cmp) for cmp in cross])

# show grad(phi(r)) = phi'(r) * rhat (parallel), the reason the cross product vanishes
s = sp.Symbol("s", positive=True)
gradcheck = [sp.simplify(grad[i] - phi(s).diff(s).subs(s, r) * [x, y, z][i] / r) for i in range(3)]
print("grad(phi(r)) - phi'(r)*rhat components:", gradcheck)
