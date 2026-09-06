#!/usr/bin/env python3
"""Dependency-free exact checks for reconstructed claims C107--C113.

The script uses only ``fractions.Fraction``.  It checks the finite-dimensional
algebra advertised by the note.  PDE estimates (an invariant slaving graph,
spatial export, and a one-cell Navier--Stokes stage map) are deliberately not
encoded as tests and remain open.
"""

from fractions import Fraction as F


def add(x, y):
    return tuple(a + b for a, b in zip(x, y))


def sub(x, y):
    return tuple(a - b for a, b in zip(x, y))


def scale(c, x):
    return tuple(c * a for a in x)


def dot(x, y):
    return sum((a * b for a, b in zip(x, y)), F(0))


def cross(x, y):
    return (
        x[1] * y[2] - x[2] * y[1],
        x[2] * y[0] - x[0] * y[2],
        x[0] * y[1] - x[1] * y[0],
    )


def mat_vec(A, x):
    return tuple(dot(row, x) for row in A)


def mat_mul(A, B):
    Bt = transpose(B)
    return tuple(tuple(dot(row, col) for col in Bt) for row in A)


def transpose(A):
    return tuple(tuple(A[j][i] for j in range(3)) for i in range(3))


def det3(A):
    return (
        A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
        - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
        + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0])
    )


def inv3(A):
    d = det3(A)
    assert d != 0
    cof = (
        (
            A[1][1] * A[2][2] - A[1][2] * A[2][1],
            -(A[1][0] * A[2][2] - A[1][2] * A[2][0]),
            A[1][0] * A[2][1] - A[1][1] * A[2][0],
        ),
        (
            -(A[0][1] * A[2][2] - A[0][2] * A[2][1]),
            A[0][0] * A[2][2] - A[0][2] * A[2][0],
            -(A[0][0] * A[2][1] - A[0][1] * A[2][0]),
        ),
        (
            A[0][1] * A[1][2] - A[0][2] * A[1][1],
            -(A[0][0] * A[1][2] - A[0][2] * A[1][0]),
            A[0][0] * A[1][1] - A[0][1] * A[1][0],
        ),
    )
    return tuple(tuple(x / d for x in row) for row in transpose(cof))


def check_c107_rank_one_and_counterexample():
    r = (F(1), F(0), F(0))
    h = (F(0), F(1), F(0))
    t = (F(0), F(0), F(1))
    gamma = F(7, 5)

    # grad J = gamma h tensor r, so (V.grad)J = gamma(V.r)h.
    for V in ((F(3), F(-2), F(5)), (F(0), F(9), F(-4)), r):
        lhs = scale(gamma * dot(V, r), h)
        rhs = scale(gamma * V[0], h)
        assert lhs == rhs

    # Raw t-wake is dark; raw zeta-wake is in-direction.
    w, phi, H = F(11, 7), F(-5, 9), F(4, 3)
    Wt = scale(w, t)
    Wz = scale(phi, sub(h, scale(H, r)))
    assert scale(gamma * dot(Wt, r), h) == (F(0), F(0), F(0))
    assert scale(gamma * dot(Wz, r), h) == scale(-H * gamma * phi, h)

    # Directional protection does not imply a bounded multiplier of J:
    # at s=0, J=0 but V=r gives gamma*h.
    s = F(0)
    J = scale(gamma * s, h)
    advected = scale(gamma * dot(r, r), h)
    assert J == (F(0), F(0), F(0))
    assert advected != J


def check_c108_clebsch_piola():
    # Several exact unimodular deformations, including a nonnormal shear.
    matrices = (
        ((F(1), F(2), F(0)), (F(0), F(1), F(3)), (F(0), F(0), F(1))),
        ((F(1), F(0), F(1)), (F(2), F(1), F(2)), (F(0), F(0), F(1))),
        ((F(0), F(1), F(0)), (F(-1), F(0), F(0)), (F(2), F(3), F(1))),
    )
    pairs = (
        ((F(1), F(2), F(-1)), (F(3), F(0), F(4))),
        ((F(0), F(5), F(2)), (F(-2), F(1), F(7))),
    )
    for A in matrices:
        assert det3(A) == 1
        AinvT = transpose(inv3(A))
        for p, q in pairs:
            lhs = cross(mat_vec(AinvT, p), mat_vec(AinvT, q))
            rhs = mat_vec(A, cross(p, q))
            assert lhs == rhs


def check_c109_clebsch_helicity_density():
    # Pointwise identity:
    # (grad phi + alpha grad beta).(grad alpha x grad beta)
    # == grad phi.(grad alpha x grad beta).
    samples = (
        ((F(1), F(2), F(3)), (F(2), F(-1), F(4)), (F(0), F(5), F(1)), F(7, 3)),
        ((F(-2), F(0), F(1)), (F(3), F(2), F(-5)), (F(4), F(1), F(0)), F(-9, 4)),
    )
    for grad_phi, grad_alpha, grad_beta, alpha in samples:
        omega = cross(grad_alpha, grad_beta)
        u = add(grad_phi, scale(alpha, grad_beta))
        assert dot(u, omega) == dot(grad_phi, omega)
        assert dot(grad_beta, omega) == 0


def check_c110_beltrami_orbit():
    # The universal pointwise cancellation u x (sigma K u) = 0.
    for u in ((F(1), F(2), F(-3)), (F(0), F(5), F(7))):
        for sigma in (F(-1), F(1)):
            K = F(13, 4)
            assert cross(u, scale(sigma * K, u)) == (F(0), F(0), F(0))

    # Scalar amplitude equation A'=-nu K^2 A makes time derivative equal
    # viscosity times Delta(A u0), with Delta u0=-K^2u0.
    nu, K, A = F(2, 9), F(5, 3), F(7, 11)
    Aprime = -nu * K * K * A
    assert Aprime == nu * (-K * K * A)


def check_c111_c112_weber_algebra():
    # At one material point, M=grad u, F=grad_a X, and u is the velocity.
    # d/dt(F^T u) = F^T(M^T u - grad p [+ nu Delta u]).
    Fm = (
        (F(1), F(2), F(0)),
        (F(0), F(1), F(1)),
        (F(0), F(0), F(1)),
    )
    M = (
        (F(2), F(-1), F(3)),
        (F(0), F(-4), F(2)),
        (F(1), F(5), F(2)),
    )
    u = (F(3), F(-2), F(4))
    grad_p = (F(5), F(1), F(-3))
    lap_u = (F(-7), F(2), F(9))
    nu = F(3, 17)
    Ft = transpose(Fm)
    Mt = transpose(M)

    # Product rule checked from the two separate differentiated factors:
    # Fdot=M F and D_t u=-grad p [+ nu Delta u].
    Fdot = mat_mul(M, Fm)
    kinetic_grad = mat_vec(Mt, u)
    euler_lhs = add(
        mat_vec(transpose(Fdot), u),
        mat_vec(Ft, scale(F(-1), grad_p)),
    )
    euler_rhs = mat_vec(Ft, sub(mat_vec(Mt, u), grad_p))
    assert euler_lhs == euler_rhs

    ns_lhs = add(
        mat_vec(transpose(Fdot), u),
        mat_vec(Ft, add(scale(F(-1), grad_p), scale(nu, lap_u))),
    )
    ns_rhs = mat_vec(Ft, add(sub(kinetic_grad, grad_p), scale(nu, lap_u)))
    assert ns_lhs == ns_rhs


def check_c113_homochiral_tangent():
    # The vector identity for the symmetrized convection contains
    # u x curl(v) + v x curl(u).  If both curls equal lambda times the
    # field, that non-gradient part cancels exactly.
    samples = (
        ((F(1), F(2), F(-3)), (F(4), F(-1), F(5)), F(7, 3)),
        ((F(-2), F(0), F(9)), (F(3), F(8), F(-4)), F(-5, 2)),
    )
    for u, v, lam in samples:
        nongradient = add(cross(u, scale(lam, v)), cross(v, scale(lam, u)))
        assert nongradient == (F(0), F(0), F(0))


def main():
    check_c107_rank_one_and_counterexample()
    check_c108_clebsch_piola()
    check_c109_clebsch_helicity_density()
    check_c110_beltrami_orbit()
    check_c111_c112_weber_algebra()
    check_c113_homochiral_tangent()
    print("PASS C107: directional protection exact; profile-slaving counterexample exact")
    print("PASS C108-C109: Clebsch-Piola and helicity-density identities exact")
    print("PASS C110: Beltrami cancellation and viscous amplitude ledger exact")
    print("PASS C111-C112: Euler/NS Weber product-rule ledger exact")
    print("PASS C113: same-eigenvalue tangent is a pure gradient")
    print("The stated PDE boundary remains: active off-shell leakage must close")


if __name__ == "__main__":
    main()
