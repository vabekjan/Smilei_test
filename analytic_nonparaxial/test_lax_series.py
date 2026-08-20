
import sympy as sp

from lax_series import (
    lax_expand_scalar,
    lax_expand_vector,
    lax_coefficients,
    paraxial_residual,
    wave_residual,
    truncate,
    magnetic_field_from_envelope,
)

I = sp.I
X, Y, Z = sp.symbols("X Y Z", real=True)
eps, k, Aamp = sp.symbols("eps k Aamp", nonzero=True)
px, py = sp.symbols("px py", real=True)

p2 = px**2 + py**2

# A simple exact paraxial seed: one transverse Fourier mode.
a0 = sp.exp(I*px*X + I*py*Y - I*p2*Z/4)


def test_seed_is_paraxial():
    assert sp.simplify(paraxial_residual(a0, (X, Y, Z))) == 0


def test_generated_series_satisfies_wave_equation_through_order_6():
    a = lax_expand_scalar(a0, (X, Y, Z), eps, order=6)
    residual = truncate(
        wave_residual(a, (X, Y, Z), eps),
        eps,
        6,
    )
    assert sp.simplify(residual) == 0


def test_generated_corrections_vanish_at_focus():
    a = lax_expand_scalar(a0, (X, Y, Z), eps, order=6)
    coeff = lax_coefficients(a, eps, 6)

    assert sp.simplify(coeff[0] - a0) == 0

    for n in (2, 4, 6):
        assert sp.simplify(coeff[n].subs(Z, 0)) == 0


def test_second_order_operator_coefficient():
    a = lax_expand_scalar(a0, (X, Y, Z), eps, order=2)
    coeff = lax_coefficients(a, eps, 2)

    lap2_a0 = (
        sp.diff(a0, X, 4)
        + 2*sp.diff(a0, X, 2, Y, 2)
        + sp.diff(a0, Y, 4)
    )

    expected = -I * Z * lap2_a0 / 64
    assert sp.simplify(coeff[2] - expected) == 0


def test_carrier_term_for_uniform_transverse_potential():
    A = lax_expand_vector(
        (Aamp, 0, 0),
        (X, Y, Z),
        eps,
        order=4,
    )

    B = magnetic_field_from_envelope(
        A,
        (X, Y, Z),
        eps,
        order=4,
        k=k,
    )

    assert sp.simplify(B[0]) == 0
    assert sp.simplify(B[1] - I*k*Aamp) == 0
    assert sp.simplify(B[2]) == 0
