"""Minimal symbolic utilities for Lax-expanded vector potentials.

The core deliberately knows nothing about Gaussian waists, Bessel cones,
gauges, normalization, or pulse envelopes. It only

1. takes a symbolic vector potential A,
2. computes B = curl(A),
3. truncates the final B field in a formal parameter eps,
4. optionally lambdifies the result for NumPy evaluation.

Coordinates may be physical Cartesian coordinates, in which case
``derivative_scales=(1, 1, 1)``, or normalized coordinates q_i with

    d/dx_i = derivative_scales[i] * d/dq_i.

The second form is useful for Lax expansions because derivative scales can
carry powers of eps.
"""

from __future__ import annotations

from typing import Iterable, Sequence

import sympy as sp

Vector3 = tuple[sp.Expr, sp.Expr, sp.Expr]


def _as_vector3(v: Sequence[sp.Expr]) -> Vector3:
    """Convert a length-three sequence to a SymPy vector tuple."""
    if len(v) != 3:
        raise ValueError("Expected exactly three vector components.")
    return tuple(sp.sympify(component) for component in v)  # type: ignore[return-value]


def curl(
    A: Sequence[sp.Expr],
    coords: Sequence[sp.Symbol],
    derivative_scales: Sequence[sp.Expr] = (1, 1, 1),
) -> Vector3:
    """Return the Cartesian curl with optionally scaled derivatives.

    If ``coords = (q1, q2, q3)``, physical derivatives are interpreted as

        D_i = derivative_scales[i] * d/dq_i.

    Physical Cartesian coordinates use the default scales.
    """
    Ax, Ay, Az = _as_vector3(A)
    if len(coords) != 3 or len(derivative_scales) != 3:
        raise ValueError("coords and derivative_scales must each have length 3.")

    x, y, z = coords
    sx, sy, sz = map(sp.sympify, derivative_scales)

    Dx = lambda expr: sx * sp.diff(expr, x)
    Dy = lambda expr: sy * sp.diff(expr, y)
    Dz = lambda expr: sz * sp.diff(expr, z)

    return (
        Dy(Az) - Dz(Ay),
        Dz(Ax) - Dx(Az),
        Dx(Ay) - Dy(Ax),
    )


def divergence(
    F: Sequence[sp.Expr],
    coords: Sequence[sp.Symbol],
    derivative_scales: Sequence[sp.Expr] = (1, 1, 1),
) -> sp.Expr:
    """Return the Cartesian divergence with the same derivative convention."""
    Fx, Fy, Fz = _as_vector3(F)
    if len(coords) != 3 or len(derivative_scales) != 3:
        raise ValueError("coords and derivative_scales must each have length 3.")

    x, y, z = coords
    sx, sy, sz = map(sp.sympify, derivative_scales)

    return (
        sx * sp.diff(Fx, x)
        + sy * sp.diff(Fy, y)
        + sz * sp.diff(Fz, z)
    )


def truncate(expr: sp.Expr, eps: sp.Symbol, order: int) -> sp.Expr:
    """Keep terms through ``eps**order`` in an expression."""
    if order < 0:
        raise ValueError("order must be non-negative.")
    return sp.series(sp.sympify(expr), eps, 0, order + 1).removeO()


def truncate_vector(
    F: Sequence[sp.Expr],
    eps: sp.Symbol,
    order: int,
) -> Vector3:
    """Apply ``truncate`` component by component."""
    Fx, Fy, Fz = _as_vector3(F)
    return (
        truncate(Fx, eps, order),
        truncate(Fy, eps, order),
        truncate(Fz, eps, order),
    )


def magnetic_field(
    A: Sequence[sp.Expr],
    coords: Sequence[sp.Symbol],
    eps: sp.Symbol,
    order: int,
    derivative_scales: Sequence[sp.Expr] = (1, 1, 1),
) -> Vector3:
    """Construct ``B = curl(A)`` and truncate the final B at the requested order."""
    return truncate_vector(curl(A, coords, derivative_scales), eps, order)


def compile_vector_field(
    F: Sequence[sp.Expr],
    args: Iterable[sp.Symbol],
):
    """Return a NumPy-compatible callable for a symbolic vector field.

    Example
    -------
    ``B_fn = compile_vector_field(B, (x, y, z, k, eps))``
    ``Bx, By, Bz = B_fn(x_values, y_values, z_values, k_value, eps_value)``
    """
    field = _as_vector3(F)
    return sp.lambdify(tuple(args), field, modules="numpy", cse=True)
