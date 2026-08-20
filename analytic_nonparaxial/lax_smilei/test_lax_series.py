"""Unit tests for the standalone Lax/Smilei backend.

Run with:
    python -m unittest -v test_lax_series.py
"""

import unittest

import numpy as np
import sympy as sp

from lax_series import (
    build_lax_field,
    evaluate_compiled_vector,
    lax_expand_scalar,
    lax_expand_vector,
    lax_hierarchy_residuals,
    vector_lax_coefficients,
    wave_residual,
    truncate,
)


class LaxBackendTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.I = sp.I
        cls.X, cls.Y, cls.Z, cls.tau = sp.symbols(
            "X Y Z tau", real=True
        )
        cls.eps, cls.k, cls.Bamp, cls.sigma = sp.symbols(
            "eps k Bamp sigma", positive=True, real=True
        )
        cls.coords = (cls.X, cls.Y, cls.Z)

    def test_focus_preserving_wave_equation(self):
        X, Y, Z = self.coords
        p2 = sp.Integer(5)
        seed = sp.exp(self.I * X + 2 * self.I * Y - self.I * p2 * Z / 4)
        series = lax_expand_scalar(seed, self.coords, self.eps, order=4)
        residual = truncate(
            wave_residual(series, self.coords, self.eps), self.eps, 4
        )
        self.assertEqual(sp.simplify(residual), 0)

    def test_specified_orders_and_hierarchy(self):
        X, Y, Z = self.coords
        rho2 = X**2 + Y**2
        # Complex-conjugated Salamin convention, matching this backend's
        # carrier exp(i*k*z-i*omega*t) and +4*i*d_Z paraxial operator.
        f = 1 / (1 + self.I * Z)
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

        A = lax_expand_vector(
            (0, 0, psi0),
            self.coords,
            self.eps,
            order=4,
            construction="specified_orders",
            order_coefficients={
                2: (0, 0, psi2),
                4: (0, 0, psi4),
            },
        )
        coefficients = vector_lax_coefficients(A, self.eps, 4)
        residuals = lax_hierarchy_residuals(coefficients, self.coords)
        for residual in residuals.values():
            for component in residual:
                self.assertEqual(sp.simplify(component), 0)

    def _linear_gaussian_build(self, pulse_application):
        X, Y, Z = self.coords
        rho2 = X**2 + Y**2
        f = 1 / (1 + self.I * Z)
        psi0 = f * sp.exp(-f * rho2)
        A0 = (0, self.I * self.Bamp * psi0 / self.k, 0)
        scales = (
            self.k * self.eps / 2,
            self.k * self.eps / 2,
            self.k * self.eps**2 / 2,
        )
        return build_lax_field(
            A0=A0,
            coords=self.coords,
            eps=self.eps,
            k=self.k,
            order=4,
            derivative_scales=scales,
            pulse=sp.exp(-self.tau**2 / self.sigma),
            pulse_application=pulse_application,
            retarded_time=self.tau,
            parameter_subs={
                self.eps: 0.3,
                self.k: 1.0,
                self.Bamp: 1.0,
                self.sigma: 20.0,
            },
        )

    def test_single_A_order_sets_complete_B_order(self):
        result = self._linear_gaussian_build("vector_potential")
        self.assertEqual(result.order, 4)
        self.assertEqual(result.magnetic_order, 5)

    def test_pulse_options_coincide_at_pulse_centre(self):
        from_A = self._linear_gaussian_build("vector_potential")
        on_B = self._linear_gaussian_build("magnetic_field")
        args = (0.2, 0.3, -0.1, 0.0)
        A_value = evaluate_compiled_vector(from_A.compiled_field, *args)
        B_value = evaluate_compiled_vector(on_B.compiled_field, *args)
        np.testing.assert_allclose(A_value, B_value, rtol=1e-12, atol=1e-12)

    def test_vectorised_evaluation(self):
        result = self._linear_gaussian_build("magnetic_field")
        values = np.linspace(-0.2, 0.2, 7)
        field = evaluate_compiled_vector(
            result.compiled_field, values, 0.5 * values, 0.0, values
        )
        self.assertEqual(field.shape, (3, 7))
        self.assertTrue(np.all(np.isfinite(field)))


if __name__ == "__main__":
    unittest.main()
