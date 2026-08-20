
"""
Generic monochromatic Lax-series construction from a paraxial vector-potential seed.

Propagation convention
----------------------
The full vector potential is

    A_full(r,t) = a(X,Y,Z) * exp(i*k*z - i*omega*t)

and the normalized envelope equation is

    (Delta_perp + 4 i d_Z + eps^2 d_Z^2) a = 0.

The user supplies only the zeroth-order paraxial envelope a0, which must satisfy

    (Delta_perp + 4 i d_Z) a0 = 0.

Boundary convention
-------------------
The higher-order terms are chosen as the Taylor expansion of the forward
Helmholtz propagator that preserves the supplied transverse profile at Z=0:

    a^(2j)(X,Y,0) = 0,  j > 0.

This fixes the otherwise arbitrary homogeneous paraxial additions.
"""

from __future__ import annotations

from typing import Iterable, Sequence
import sympy as sp


Vector3 = tuple[sp.Expr, sp.Expr, sp.Expr]


def _as_vector3(v: Sequence[sp.Expr]) -> Vector3:
    if len(v) != 3:
        raise ValueError("Expected exactly three vector components.")
    return tuple(sp.sympify(c) for c in v)


def truncate(expr: sp.Expr, eps: sp.Symbol, order: int) -> sp.Expr:
    """Keep terms through eps**order."""
    if order < 0:
        raise ValueError("order must be non-negative.")
    return sp.expand(sp.series(sp.sympify(expr), eps, 0, order + 1).removeO())


def truncate_vector(
    field: Sequence[sp.Expr],
    eps: sp.Symbol,
    order: int,
) -> Vector3:
    return tuple(truncate(c, eps, order) for c in _as_vector3(field))


def transverse_laplacian(
    expr: sp.Expr,
    X: sp.Symbol,
    Y: sp.Symbol,
) -> sp.Expr:
    """Dimensionless transverse Laplacian Delta_perp."""
    return sp.diff(expr, X, 2) + sp.diff(expr, Y, 2)


def paraxial_residual(
    expr: sp.Expr,
    coords: Sequence[sp.Symbol],
) -> sp.Expr:
    """Return (Delta_perp + 4 i d_Z) expr."""
    X, Y, Z = coords
    return transverse_laplacian(expr, X, Y) + 4 * sp.I * sp.diff(expr, Z)


def wave_residual(
    expr: sp.Expr,
    coords: Sequence[sp.Symbol],
    eps: sp.Symbol,
) -> sp.Expr:
    """Return the full normalized envelope-wave-equation residual."""
    X, Y, Z = coords
    return (
        transverse_laplacian(expr, X, Y)
        + 4 * sp.I * sp.diff(expr, Z)
        + eps**2 * sp.diff(expr, Z, 2)
    )


def _relative_propagator_series(
    eps: sp.Symbol,
    Z: sp.Symbol,
    order: int,
) -> tuple[sp.Expr, sp.Symbol]:
    r"""
    Return the forward Helmholtz correction operator as a polynomial in q.

    q is a placeholder for Delta_perp.

    Exact forward envelope propagator:
        exp[ 2 i Z/eps^2 * (sqrt(1 + eps^2 q/4) - 1) ]

    Paraxial propagator:
        exp[ i Z q/4 ]

    Their ratio is expanded in eps.
    """
    if order < 0:
        raise ValueError("order must be non-negative.")

    even_order = order if order % 2 == 0 else order - 1
    q = sp.Dummy("q")

    sqrt_series = sp.series(
        sp.sqrt(1 + eps**2 * q / 4),
        eps,
        0,
        even_order + 4,
    ).removeO()

    exact_exponent = 2 * sp.I * Z / eps**2 * (sqrt_series - 1)
    paraxial_exponent = sp.I * Z * q / 4

    relative_exponent = truncate(
        sp.expand(exact_exponent - paraxial_exponent),
        eps,
        even_order,
    )

    correction = truncate(
        sp.exp(relative_exponent),
        eps,
        even_order,
    )

    return sp.expand(correction), q


def lax_operator_series(
    eps: sp.Symbol,
    Z: sp.Symbol,
    order: int,
) -> tuple[sp.Expr, sp.Symbol]:
    """Public view of the forward-Helmholtz correction operator."""
    return _relative_propagator_series(eps, Z, order)


def lax_expand_scalar(
    a0: sp.Expr,
    coords: Sequence[sp.Symbol],
    eps: sp.Symbol,
    order: int,
    *,
    simplify: bool = False,
) -> sp.Expr:
    r"""
    Generate a = a^(0) + eps^2 a^(2) + ... from one paraxial seed a0.

    The seed must be independent of eps and should satisfy the paraxial
    equation. The generated higher-order corrections vanish at Z=0.
    """
    X, Y, Z = coords
    a0 = sp.sympify(a0)

    if a0.has(eps):
        raise ValueError("The zeroth-order seed a0 must not contain eps.")

    correction, q = _relative_propagator_series(eps, Z, order)
    polynomial = sp.Poly(correction, eps, q)

    lap_cache = {0: a0}

    def lap_power(n: int) -> sp.Expr:
        if n not in lap_cache:
            highest = max(lap_cache)
            value = lap_cache[highest]
            for p in range(highest + 1, n + 1):
                value = transverse_laplacian(value, X, Y)
                lap_cache[p] = value
        return lap_cache[n]

    result = sp.Integer(0)

    for (eps_power, q_power), coefficient in polynomial.terms():
        if eps_power <= order:
            result += (
                coefficient
                * eps**eps_power
                * lap_power(q_power)
            )

    result = sp.expand(result)
    return sp.simplify(result) if simplify else result


def lax_expand_vector(
    A0: Sequence[sp.Expr],
    coords: Sequence[sp.Symbol],
    eps: sp.Symbol,
    order: int,
    *,
    simplify: bool = False,
) -> Vector3:
    """Apply the scalar Lax generator independently to Ax0, Ay0, Az0."""
    return tuple(
        lax_expand_scalar(
            component,
            coords,
            eps,
            order,
            simplify=simplify,
        )
        for component in _as_vector3(A0)
    )


def lax_coefficients(
    expr: sp.Expr,
    eps: sp.Symbol,
    order: int,
) -> dict[int, sp.Expr]:
    """Return the nonzero coefficients of eps through the requested order."""
    expanded = sp.expand(expr)
    out: dict[int, sp.Expr] = {}

    for n in range(order + 1):
        coefficient = sp.simplify(expanded.coeff(eps, n))
        if coefficient != 0:
            out[n] = coefficient

    return out


def curl_envelope(
    A: Sequence[sp.Expr],
    coords: Sequence[sp.Symbol],
    derivative_scales: Sequence[sp.Expr] = (1, 1, 1),
) -> Vector3:
    r"""
    Curl of the spatial envelope.

    derivative_scales means

        d/dx = sx d/dX,
        d/dy = sy d/dY,
        d/dz = sz d/dZ.
    """
    Ax, Ay, Az = _as_vector3(A)
    X, Y, Z = coords
    sx, sy, sz = map(sp.sympify, derivative_scales)

    Dx = lambda expr: sx * sp.diff(expr, X)
    Dy = lambda expr: sy * sp.diff(expr, Y)
    Dz = lambda expr: sz * sp.diff(expr, Z)

    return (
        Dy(Az) - Dz(Ay),
        Dz(Ax) - Dx(Az),
        Dx(Ay) - Dy(Ax),
    )


def magnetic_field_from_envelope(
    A: Sequence[sp.Expr],
    coords: Sequence[sp.Symbol],
    eps: sp.Symbol,
    order: int,
    k: sp.Expr,
    derivative_scales: Sequence[sp.Expr] = (1, 1, 1),
) -> Vector3:
    r"""
    Magnetic-field envelope for carrier exp(i*k*z - i*omega*t).

    Since

        A_full = A * exp(i*k*z),

    the field envelope is

        B = curl(A) + i*k*z_hat x A.

    The second term is essential for generic transverse vector potentials.
    """
    Ax, Ay, Az = _as_vector3(A)
    Bcurl = curl_envelope(A, coords, derivative_scales)

    carrier_term = (
        -sp.I * k * Ay,
        sp.I * k * Ax,
        sp.Integer(0),
    )

    B = tuple(
        Bcurl[j] + carrier_term[j]
        for j in range(3)
    )

    return truncate_vector(B, eps, order)


def full_divergence_residual(
    F: Sequence[sp.Expr],
    coords: Sequence[sp.Symbol],
    k: sp.Expr,
    derivative_scales: Sequence[sp.Expr] = (1, 1, 1),
) -> sp.Expr:
    r"""
    Divergence residual for F_full = F * exp(i*k*z).

    Returns

        div_envelope(F) + i*k*Fz.
    """
    Fx, Fy, Fz = _as_vector3(F)
    X, Y, Z = coords
    sx, sy, sz = map(sp.sympify, derivative_scales)

    return (
        sx * sp.diff(Fx, X)
        + sy * sp.diff(Fy, Y)
        + sz * sp.diff(Fz, Z)
        + sp.I * k * Fz
    )


def compile_vector_field(
    field: Sequence[sp.Expr],
    args: Iterable[sp.Symbol],
):
    """Compile a symbolic vector field to a NumPy-compatible callable."""
    return sp.lambdify(
        tuple(args),
        _as_vector3(field),
        modules="numpy",
        cse=True,
    )
