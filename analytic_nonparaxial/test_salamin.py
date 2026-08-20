"""Reference test for a radially polarized non-paraxial Gaussian beam.

Reference:
Y. I. Salamin, Opt. Lett. 31, 2619-2621 (2006).

The same field through O(eps**5) is used in
M. Jirka and S. V. Bulanov, Phys. Rev. Lett. 133, 125001 (2024),
with the propagation axis relabelled.

Salamin starts from an axially polarized vector potential

    A = z_hat A0 Psi exp[i(omega t - k z)]

with

    Psi = Psi0 + eps**2 Psi2 + eps**4 Psi4 + ...

Here the common carrier and time factor are omitted. They do not affect the
benchmark because the potential has only a longitudinal component and
B_theta = -dA_z/dr uses a transverse derivative.
"""

from functools import lru_cache

import numpy as np
import sympy as sp

from lax_beams import compile_vector_field, curl, divergence, truncate, truncate_vector

I = sp.I

X, Y, Z = sp.symbols("X Y Z", real=True)
R = sp.symbols("R", nonnegative=True, real=True)
eps, k, A0 = sp.symbols("eps k A0", positive=True, real=True)

# For X=x/w0, Y=y/w0, Z=z/zR with eps=2/(k*w0):
SCALES = (k * eps / 2, k * eps / 2, k * eps**2 / 2)


def salamin_vector_potential():
    """Longitudinal vector-potential envelope through Psi4."""
    rho2 = X**2 + Y**2
    f = I / (I + Z)
    gaussian = sp.exp(-f * rho2)

    psi0 = f * gaussian

    psi2 = (
        sp.Rational(1, 2) * f**2
        - sp.Rational(1, 4) * rho2**2 * f**4
    ) * gaussian

    psi4 = (
        sp.Rational(3, 8) * f**3
        - sp.Rational(3, 16) * rho2**2 * f**5
        - sp.Rational(1, 8) * rho2**3 * f**6
        + sp.Rational(1, 32) * rho2**4 * f**7
    ) * gaussian

    psi = psi0 + eps**2 * psi2 + eps**4 * psi4
    return (sp.Integer(0), sp.Integer(0), A0 * psi)


def salamin_reference_btheta():
    """Published complex B_theta envelope normalized by k*A0 = E0/c."""
    f = I / (I + Z)
    gaussian = sp.exp(-f * R**2)

    b1 = R * f**2

    b3 = (
        sp.Rational(1, 2) * R * f**3
        + sp.Rational(1, 2) * R**3 * f**4
        - sp.Rational(1, 4) * R**5 * f**5
    )

    b5 = (
        sp.Rational(3, 8) * R * f**4
        + sp.Rational(3, 8) * R**3 * f**5
        + sp.Rational(3, 16) * R**5 * f**6
        - sp.Rational(1, 4) * R**7 * f**7
        + sp.Rational(1, 32) * R**9 * f**8
    )

    return gaussian * (eps * b1 + eps**3 * b3 + eps**5 * b5)


@lru_cache(maxsize=1)
def full_constructed_B():
    """Compute the curl once and reuse it in all tests."""
    return curl(salamin_vector_potential(), (X, Y, Z), SCALES)


def constructed_B(order):
    """Truncate the already-constructed magnetic field at a chosen Lax order."""
    return truncate_vector(full_constructed_B(), eps, order)


@lru_cache(maxsize=1)
def full_constructed_btheta_on_x_axis():
    """Return B_theta/(k*A0) along Y=0, X=R.

    On the positive X axis, e_theta=e_y, hence B_theta=B_y. Rotational
    symmetry then identifies the radial coefficient without introducing a
    symbolic sqrt(X**2 + Y**2).
    """
    B = full_constructed_B()
    return sp.simplify((B[1] / (k * A0)).subs({X: R, Y: 0}))


def test_reference_orders():
    """Reproduce the published B_theta at orders 1, 3, and 5."""
    ours_full = full_constructed_btheta_on_x_axis()
    reference = salamin_reference_btheta()

    for order in (1, 3, 5):
        ours = truncate(ours_full, eps, order)
        expected = truncate(reference, eps, order)
        assert sp.simplify(ours - expected) == 0


def test_divergence_through_order_5():
    """Verify div(B)=0 through the retained O(eps**5) field."""
    B5 = constructed_B(5)
    residual = truncate(divergence(B5, (X, Y, Z), SCALES), eps, 5)
    assert sp.simplify(residual) == 0


def test_lambdified_field():
    """Smoke-test NumPy evaluation of the symbolic result."""
    B5 = constructed_B(5)
    B_numeric = tuple(component.subs({k: 1, A0: 1}) for component in B5)
    B_fn = compile_vector_field(B_numeric, (X, Y, Z, eps))

    x = np.array([0.0, 0.25, 0.5])
    y = np.zeros_like(x)
    z = np.zeros_like(x)
    bx, by, bz = B_fn(x, y, z, 0.3)

    for component in (bx, by, bz):
        assert np.all(np.isfinite(np.asarray(component, dtype=complex)))


if __name__ == "__main__":
    ours_full = full_constructed_btheta_on_x_axis()
    reference = salamin_reference_btheta()

    print("Symbolic comparison with published B_theta:")
    for order in (1, 3, 5):
        difference = sp.simplify(
            truncate(ours_full, eps, order)
            - truncate(reference, eps, order)
        )
        print(f"  order {order}: difference = {difference}")

    B5 = constructed_B(5)
    residual = sp.simplify(
        truncate(divergence(B5, (X, Y, Z), SCALES), eps, 5)
    )
    print(f"  div(B) through order 5: {residual}")

    test_lambdified_field()
    print("  lambdify smoke test: passed")
