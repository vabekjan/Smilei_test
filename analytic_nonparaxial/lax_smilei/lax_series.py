"""Symbolic Lax-series backend for custom Smilei laser profiles.

The user supplies the symbolic vector-potential seed (and optionally explicit
higher-order coefficients) in the Smilei namelist.  This module constructs and
checks the Lax series, builds the magnetic field, applies one of two pulse
prescriptions, lambdifies the result, and writes an optional readable report.

Carrier convention
------------------
    A_full = A(X, Y, Z, tau) * exp(i*k*z - i*omega*t)

Normalized Lax equation
-----------------------
    (Delta_perp + 4*i*d_Z + eps**2*d_Z**2) A = 0

For an A series complete through even order N, the magnetic series is complete
through order N+1.  Formal order-N+2 pieces produced by longitudinal
derivatives are deliberately discarded because A^(N+2) is unavailable.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable, Mapping, Sequence

import numpy as np
import sympy as sp


Vector3 = tuple[sp.Expr, sp.Expr, sp.Expr]
COMPONENT_NAMES = ("x", "y", "z")


class LaxValidationError(ValueError):
    """Raised after a symbolic report has been assembled for a failed build."""

    def __init__(self, message: str, report: str = ""):
        super().__init__(message)
        self.report = report


@dataclass(frozen=True)
class LaxFieldBuild:
    """Result of :func:`build_lax_field`."""

    vector_potential: Vector3
    magnetic_envelope: Vector3
    compiled_field: Callable
    report: str
    order: int
    magnetic_order: int
    nonzero_magnetic_orders: tuple[int, ...]


def _as_vector3(v: Sequence[sp.Expr]) -> Vector3:
    if len(v) != 3:
        raise ValueError("Expected exactly three vector components.")
    return tuple(sp.sympify(component) for component in v)  # type: ignore[return-value]


def _require_even_order(order: int) -> None:
    if order < 0 or order % 2:
        raise ValueError("The vector-potential Lax order must be non-negative and even.")


def truncate(expr: sp.Expr, eps: sp.Symbol, order: int) -> sp.Expr:
    """Keep formal powers through ``eps**order``."""
    if order < 0:
        raise ValueError("order must be non-negative.")
    return sp.expand(sp.series(sp.sympify(expr), eps, 0, order + 1).removeO())


def truncate_vector(field: Sequence[sp.Expr], eps: sp.Symbol, order: int) -> Vector3:
    return tuple(truncate(component, eps, order) for component in _as_vector3(field))


def transverse_laplacian(expr: sp.Expr, X: sp.Symbol, Y: sp.Symbol) -> sp.Expr:
    return sp.diff(expr, X, 2) + sp.diff(expr, Y, 2)


def paraxial_residual(expr: sp.Expr, coords: Sequence[sp.Symbol]) -> sp.Expr:
    X, Y, Z = coords
    return transverse_laplacian(expr, X, Y) + 4 * sp.I * sp.diff(expr, Z)


def wave_residual(
    expr: sp.Expr,
    coords: Sequence[sp.Symbol],
    eps: sp.Symbol,
) -> sp.Expr:
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
    """Forward-Helmholtz correction relative to paraxial propagation."""
    _require_even_order(order)
    q = sp.Dummy("q")

    sqrt_series = sp.series(
        sp.sqrt(1 + eps**2 * q / 4), eps, 0, order + 4
    ).removeO()
    exact_exponent = 2 * sp.I * Z / eps**2 * (sqrt_series - 1)
    paraxial_exponent = sp.I * Z * q / 4
    relative_exponent = truncate(
        sp.expand(exact_exponent - paraxial_exponent), eps, order
    )
    correction = truncate(sp.exp(relative_exponent), eps, order)
    return sp.expand(correction), q


def lax_operator_series(
    eps: sp.Symbol,
    Z: sp.Symbol,
    order: int,
) -> tuple[sp.Expr, sp.Symbol]:
    return _relative_propagator_series(eps, Z, order)


def _focus_preserving_scalar(
    a0: sp.Expr,
    coords: Sequence[sp.Symbol],
    eps: sp.Symbol,
    order: int,
    *,
    simplify: bool,
) -> sp.Expr:
    X, Y, Z = coords
    correction, q = _relative_propagator_series(eps, Z, order)
    polynomial = sp.Poly(correction, eps, q)
    lap_cache: dict[int, sp.Expr] = {0: a0}

    def lap_power(n: int) -> sp.Expr:
        if n not in lap_cache:
            highest = max(lap_cache)
            value = lap_cache[highest]
            for power in range(highest + 1, n + 1):
                value = transverse_laplacian(value, X, Y)
                lap_cache[power] = value
        return lap_cache[n]

    result = sp.Integer(0)
    for (eps_power, q_power), coefficient in polynomial.terms():
        result += coefficient * eps**eps_power * lap_power(q_power)

    result = sp.expand(result)
    return sp.simplify(result) if simplify else result


def lax_expand_scalar(
    a0: sp.Expr,
    coords: Sequence[sp.Symbol],
    eps: sp.Symbol,
    order: int,
    *,
    construction: str = "focus_preserving",
    order_coefficients: Mapping[int, sp.Expr] | None = None,
    simplify: bool = False,
) -> sp.Expr:
    """Construct one scalar Lax series using either supported prescription."""
    _require_even_order(order)
    a0 = sp.sympify(a0)
    if a0.has(eps):
        raise ValueError("A^(0) must not contain eps.")

    if construction == "focus_preserving":
        if order_coefficients:
            raise ValueError(
                "order_coefficients are not used by focus_preserving construction."
            )
        return _focus_preserving_scalar(
            a0, coords, eps, order, simplify=simplify
        )

    if construction != "specified_orders":
        raise ValueError(
            "construction must be 'focus_preserving' or 'specified_orders'."
        )

    coefficients = dict(order_coefficients or {})
    expected = set(range(2, order + 1, 2))
    actual = set(coefficients)
    if actual != expected:
        raise ValueError(
            "specified_orders requires exactly the coefficients "
            f"{sorted(expected)}; received {sorted(actual)}."
        )

    result = a0
    for power in sorted(coefficients):
        coefficient = sp.sympify(coefficients[power])
        if coefficient.has(eps):
            raise ValueError(f"A^({power}) must not contain eps.")
        result += eps**power * coefficient

    result = sp.expand(result)
    return sp.simplify(result) if simplify else result


def lax_expand_vector(
    A0: Sequence[sp.Expr],
    coords: Sequence[sp.Symbol],
    eps: sp.Symbol,
    order: int,
    *,
    construction: str = "focus_preserving",
    order_coefficients: Mapping[int, Sequence[sp.Expr]] | None = None,
    simplify: bool = False,
) -> Vector3:
    """Construct a three-component Lax vector-potential series."""
    A0 = _as_vector3(A0)

    if construction == "focus_preserving":
        if order_coefficients:
            raise ValueError(
                "order_coefficients are not used by focus_preserving construction."
            )
        return tuple(
            lax_expand_scalar(
                component,
                coords,
                eps,
                order,
                construction=construction,
                simplify=simplify,
            )
            for component in A0
        )

    vector_coefficients = {
        power: _as_vector3(vector)
        for power, vector in dict(order_coefficients or {}).items()
    }
    return tuple(
        lax_expand_scalar(
            A0[index],
            coords,
            eps,
            order,
            construction=construction,
            order_coefficients={
                power: vector[index]
                for power, vector in vector_coefficients.items()
            },
            simplify=simplify,
        )
        for index in range(3)
    )


def lax_coefficients(
    expr: sp.Expr,
    eps: sp.Symbol,
    order: int,
) -> dict[int, sp.Expr]:
    expanded = sp.expand(expr)
    return {
        power: coefficient
        for power in range(order + 1)
        if (coefficient := sp.simplify(expanded.coeff(eps, power))) != 0
    }


def vector_lax_coefficients(
    field: Sequence[sp.Expr],
    eps: sp.Symbol,
    order: int,
) -> dict[int, Vector3]:
    vector = _as_vector3(field)
    return {
        power: tuple(
            sp.simplify(sp.expand(component).coeff(eps, power))
            for component in vector
        )
        for power in range(0, order + 1, 2)
    }


def lax_hierarchy_residuals(
    coefficients: Mapping[int, Sequence[sp.Expr]],
    coords: Sequence[sp.Symbol],
) -> dict[int, Vector3]:
    """Return coefficient-wise residuals of the normalized Lax hierarchy."""
    powers = sorted(coefficients)
    if not powers or powers[0] != 0 or any(power < 0 or power % 2 for power in powers):
        raise ValueError("coefficients must contain order 0 and only even orders.")
    expected = list(range(0, powers[-1] + 1, 2))
    if powers != expected:
        raise ValueError(f"coefficients must contain every order in {expected}.")

    vectors = {power: _as_vector3(coefficients[power]) for power in powers}
    _, _, Z = coords
    residuals: dict[int, Vector3] = {}
    for power in powers:
        residuals[power] = tuple(
            sp.simplify(
                paraxial_residual(vectors[power][index], coords)
                + (
                    sp.diff(vectors[power - 2][index], Z, 2)
                    if power >= 2
                    else 0
                )
            )
            for index in range(3)
        )
    return residuals


def curl_envelope(
    A: Sequence[sp.Expr],
    coords: Sequence[sp.Symbol],
    derivative_scales: Sequence[sp.Expr] = (1, 1, 1),
) -> Vector3:
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
    """Complex magnetic envelope for carrier ``exp(i*k*z-i*omega*t)``."""
    Ax, Ay, _ = _as_vector3(A)
    Bcurl = curl_envelope(A, coords, derivative_scales)
    carrier_term = (-sp.I * k * Ay, sp.I * k * Ax, sp.Integer(0))
    return truncate_vector(
        tuple(Bcurl[index] + carrier_term[index] for index in range(3)),
        eps,
        order,
    )


def magnetic_field_from_pulsed_potential(
    A: Sequence[sp.Expr],
    pulse: sp.Expr,
    retarded_time: sp.Symbol,
    group_velocity: sp.Expr,
    coords: Sequence[sp.Symbol],
    eps: sp.Symbol,
    order: int,
    k: sp.Expr,
    derivative_scales: Sequence[sp.Expr] = (1, 1, 1),
) -> Vector3:
    """Magnetic envelope when ``pulse(tau)`` multiplies A before curl.

    With tau=t-z/v_g, the product rule gives

        B_pulse = pulse*B_mono - pulse'/v_g * z_hat cross A.
    """
    Ax, Ay, _ = _as_vector3(A)
    pulse = sp.sympify(pulse)
    velocity = sp.sympify(group_velocity)
    if velocity == 0:
        raise ValueError("group_velocity must be nonzero.")

    Bmono = magnetic_field_from_envelope(
        A, coords, eps, order, k, derivative_scales
    )
    pulse_derivative = sp.diff(pulse, retarded_time)
    z_cross_A = (-Ay, Ax, sp.Integer(0))
    result = tuple(
        pulse * Bmono[index]
        - pulse_derivative * z_cross_A[index] / velocity
        for index in range(3)
    )
    return truncate_vector(result, eps, order)


def apply_pulse_to_magnetic_field(
    B: Sequence[sp.Expr],
    pulse: sp.Expr,
) -> Vector3:
    pulse = sp.sympify(pulse)
    return tuple(pulse * component for component in _as_vector3(B))


def full_divergence_residual_retarded(
    F: Sequence[sp.Expr],
    coords: Sequence[sp.Symbol],
    k: sp.Expr,
    derivative_scales: Sequence[sp.Expr],
    retarded_time: sp.Symbol,
    group_velocity: sp.Expr,
) -> sp.Expr:
    """Divergence for an envelope depending on tau=t-z/v_g and carrier z."""
    Fx, Fy, Fz = _as_vector3(F)
    X, Y, Z = coords
    sx, sy, sz = map(sp.sympify, derivative_scales)
    velocity = sp.sympify(group_velocity)
    return (
        sx * sp.diff(Fx, X)
        + sy * sp.diff(Fy, Y)
        + sz * sp.diff(Fz, Z)
        - sp.diff(Fz, retarded_time) / velocity
        + sp.I * k * Fz
    )


def compile_vector_field(field: Sequence[sp.Expr], args: Iterable[sp.Symbol]):
    return sp.lambdify(
        tuple(args), _as_vector3(field), modules="numpy", cse=True
    )


def evaluate_compiled_vector(compiled_field: Callable, *args) -> np.ndarray:
    """Evaluate a lambdified vector and broadcast symbolic scalar zeros."""
    arrays = tuple(np.asarray(argument) for argument in args)
    shape = np.broadcast(*arrays).shape
    raw = compiled_field(*args)
    return np.stack(
        [
            np.broadcast_to(np.asarray(component, dtype=complex), shape)
            for component in raw
        ],
        axis=0,
    )


def _format_vector(vector: Sequence[sp.Expr], indent: str = "    ") -> list[str]:
    lines: list[str] = []
    for name, component in zip(COMPONENT_NAMES, _as_vector3(vector)):
        pretty = sp.pretty(component, use_unicode=False, wrap_line=False)
        lines.append(f"{indent}{name}:")
        lines.extend(f"{indent}    {line}" for line in pretty.splitlines())
    return lines


def _write_report(path: str | Path | None, report: str) -> None:
    if path is not None:
        Path(path).write_text(report, encoding="utf-8")


def build_lax_field(
    *,
    A0: Sequence[sp.Expr],
    coords: Sequence[sp.Symbol],
    eps: sp.Symbol,
    k: sp.Expr,
    order: int,
    derivative_scales: Sequence[sp.Expr],
    construction: str = "focus_preserving",
    order_coefficients: Mapping[int, Sequence[sp.Expr]] | None = None,
    pulse: sp.Expr = sp.Integer(1),
    pulse_application: str = "magnetic_field",
    retarded_time: sp.Symbol | None = None,
    group_velocity: sp.Expr = sp.Integer(1),
    parameter_subs: Mapping[sp.Symbol, sp.Expr] | None = None,
    validate: bool = True,
    report_path: str | Path | None = None,
) -> LaxFieldBuild:
    """Build, validate, report, and compile one pulsed Lax magnetic field."""
    _require_even_order(order)
    if len(coords) != 3 or len(derivative_scales) != 3:
        raise ValueError("coords and derivative_scales must each have length 3.")
    if pulse_application not in {"vector_potential", "magnetic_field"}:
        raise ValueError(
            "pulse_application must be 'vector_potential' or 'magnetic_field'."
        )
    if retarded_time is None:
        raise ValueError("retarded_time must be supplied explicitly.")

    A0 = _as_vector3(A0)
    A = lax_expand_vector(
        A0,
        coords,
        eps,
        order,
        construction=construction,
        order_coefficients=order_coefficients,
    )
    A_coefficients = vector_lax_coefficients(A, eps, order)
    magnetic_order = order + 1
    Bmono = magnetic_field_from_envelope(
        A,
        coords,
        eps,
        magnetic_order,
        k,
        derivative_scales,
    )

    if pulse_application == "vector_potential":
        Bpulse = magnetic_field_from_pulsed_potential(
            A,
            pulse,
            retarded_time,
            group_velocity,
            coords,
            eps,
            magnetic_order,
            k,
            derivative_scales,
        )
    else:
        Bpulse = truncate_vector(
            apply_pulse_to_magnetic_field(Bmono, pulse),
            eps,
            magnetic_order,
        )

    nonzero_orders = tuple(
        power
        for power in range(magnetic_order + 1)
        if any(
            sp.simplify(sp.expand(component).coeff(eps, power)) != 0
            for component in Bpulse
        )
    )

    errors: list[str] = []
    warnings: list[str] = []
    hierarchy: dict[int, Vector3] | None = None
    focus_checks: dict[int, Vector3] = {}
    divergence = sp.Symbol("validation_skipped")

    if validate:
        hierarchy = lax_hierarchy_residuals(A_coefficients, coords)
        for power, residual in hierarchy.items():
            for index, component in enumerate(residual):
                if sp.simplify(component) != 0:
                    errors.append(
                        f"Lax hierarchy failed at A^({power}) component "
                        f"{COMPONENT_NAMES[index]}."
                    )

        if construction == "focus_preserving":
            Z = coords[2]
            for power in range(2, order + 1, 2):
                focus_checks[power] = tuple(
                    sp.simplify(component.subs(Z, 0))
                    for component in A_coefficients[power]
                )
                if any(component != 0 for component in focus_checks[power]):
                    errors.append(f"Focus condition failed at A^({power}).")

        divergence = sp.simplify(
            truncate(
                full_divergence_residual_retarded(
                    Bpulse,
                    coords,
                    k,
                    derivative_scales,
                    retarded_time,
                    group_velocity,
                ),
                eps,
                magnetic_order,
            )
        )
        if divergence != 0:
            message = (
                "The pulsed magnetic field has a nonzero divergence residual "
                f"through order {magnetic_order}."
            )
            if pulse_application == "vector_potential":
                errors.append(message)
            else:
                warnings.append(message + " This is expected for post-B pulsing.")

    substitutions = dict(parameter_subs or {})
    # Substitution is intentionally not followed by global simplify: the
    # expressions are already expanded/truncated, and simplify can dominate
    # namelist start-up time on every MPI process.
    Bnumeric = tuple(component.subs(substitutions) for component in Bpulse)
    compile_args = tuple(coords) + (retarded_time,)
    unresolved = set().union(*(component.free_symbols for component in Bnumeric))
    unresolved -= set(compile_args)
    if unresolved:
        errors.append(
            "Unresolved symbols before compilation: "
            + ", ".join(sorted(str(symbol) for symbol in unresolved))
        )

    compiled = None
    smoke_status = "SKIPPED"
    if not unresolved:
        try:
            compiled = compile_vector_field(Bnumeric, compile_args)
            if validate:
                smoke_status = "PASS"
                scalar_raw = compiled(0.0, 0.0, 0.0, 0.0)
                vector_raw = compiled(
                    np.asarray([0.0, 0.1]),
                    np.asarray([0.0, 0.1]),
                    np.asarray([0.0, 0.1]),
                    np.asarray([0.0, 0.1]),
                )
                for raw in (scalar_raw, vector_raw):
                    for component in raw:
                        if not np.all(np.isfinite(np.asarray(component, dtype=complex))):
                            errors.append("Compiled-field smoke test returned non-finite values.")
                            smoke_status = "FAIL"
                            break
        except Exception as exc:  # report a friendly error before stopping Smilei
            errors.append(f"Compilation or smoke test failed: {type(exc).__name__}: {exc}")
            smoke_status = "FAIL"

    # Only rank zero requests validation/reporting in the Smilei namelist.
    # Avoid formatting large symbolic expressions on all other MPI processes.
    make_full_report = validate or report_path is not None
    status = "FAIL" if errors else ("WARNING" if warnings else "PASS")
    lines = [
        "LAX SYMBOLIC INPUT REPORT",
        "=" * 72,
        f"Generated (UTC): {datetime.now(timezone.utc).isoformat()}",
        f"Status: {status}",
        "",
        "Configuration",
        "-" * 72,
        f"construction:       {construction}",
        f"pulse application:  {pulse_application}",
        f"A expansion order:  {order}",
        f"B complete through: {magnetic_order}",
        f"nonzero B orders:   {list(nonzero_orders)}",
        f"incomplete order {order + 2} terms: discarded",
        "",
        "Pulse",
        "-" * 72,
        f"g(tau)  = {sp.sstr(pulse)}",
        f"g'(tau) = {sp.sstr(sp.diff(pulse, retarded_time))}",
        f"group velocity = {sp.sstr(group_velocity)}",
        "",
        "Vector-potential coefficients",
        "-" * 72,
    ]
    if make_full_report:
        for power in range(0, order + 1, 2):
            lines.append(f"A^({power})")
            lines.extend(_format_vector(A_coefficients[power]))
    else:
        lines.append("Full symbolic formatting skipped on this process.")

    lines.extend(["", "Validation", "-" * 72])
    if not validate:
        lines.append("Symbolic validation skipped on this process.")
    else:
        assert hierarchy is not None
        for power, residual in hierarchy.items():
            passed = all(sp.simplify(component) == 0 for component in residual)
            lines.append(f"Lax hierarchy A^({power}): {'PASS' if passed else 'FAIL'}")
            if not passed:
                lines.extend(_format_vector(residual))
        for power, values in focus_checks.items():
            passed = all(component == 0 for component in values)
            lines.append(f"Focus condition A^({power}): {'PASS' if passed else 'FAIL'}")
        lines.append(f"div(B) through order {magnetic_order}: {sp.sstr(divergence)}")
        lines.append(f"Compiled scalar/vector smoke tests: {smoke_status}")

    if warnings:
        lines.extend(["", "Warnings", "-" * 72])
        lines.extend(f"- {message}" for message in warnings)
    if errors:
        lines.extend(["", "Errors", "-" * 72])
        lines.extend(f"- {message}" for message in errors)

    report = "\n".join(lines) + "\n"
    _write_report(report_path, report)

    if errors or compiled is None:
        raise LaxValidationError(
            "Lax field construction failed; see the symbolic report.", report
        )

    return LaxFieldBuild(
        vector_potential=A,
        magnetic_envelope=Bpulse,
        compiled_field=compiled,
        report=report,
        order=order,
        magnetic_order=magnetic_order,
        nonzero_magnetic_orders=nonzero_orders,
    )
