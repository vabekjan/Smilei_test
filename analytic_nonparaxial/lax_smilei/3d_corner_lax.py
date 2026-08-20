"""Smilei input: vectorised oblique Lax-series Gaussian laser injection.

The symbolic vector potential and pulse are intentionally defined in this
namelist.  ``lax_series.py`` supplies construction, checks, compilation, and a
rank-zero report named ``lax_symbolic_report.log`` in the run directory.
"""

from math import pi
from pathlib import Path
import traceback

import numpy as np
import sympy as sp

from lax_series import (
    LaxValidationError,
    build_lax_field,
    evaluate_compiled_vector,
)


# -----------------------------------------------------------------------------
# Smilei and laser parameters
# -----------------------------------------------------------------------------
sim_scale = np.asarray([10.0, 10.0, 10.0])

l0 = 2.0 * np.pi
t0 = l0
Lsim = list(l0 * sim_scale)
Tsim = 2.0 * 25.0 * t0
resx = 12.0
rest = 22.0

a0 = 1.0
omega = 1.0
k0 = omega  # vacuum, c=1
waist = 2.5 * l0
fwhm = 6.0 * t0
group_velocity = 1.0

# Time at which the pulse centre reaches the focus.
t_offset = 11.0 + 5.0 * fwhm

ang = [pi / 4.0, -pi / 4.0]
focus = [0.5 * Lsim[0], 0.5 * Lsim[1], 0.5 * Lsim[2]]


# -----------------------------------------------------------------------------
# Readable Lax configuration
# -----------------------------------------------------------------------------
lax_order = 4
lax_construction = "focus_preserving"       # or "specified_orders"
pulse_application = "magnetic_field"        # or "vector_potential"


# -----------------------------------------------------------------------------
# Symbolic input supplied by the user
# Backend axes (X,Y,Z): transverse, transverse, propagation.
# -----------------------------------------------------------------------------
I = sp.I
X, Y, Z, tau = sp.symbols("X Y Z tau", real=True)
eps, k, Bamp, sigma = sp.symbols(
    "eps k Bamp sigma", positive=True, real=True
)

rho2 = X**2 + Y**2
f = 1 / (1 + I * Z)
psi0 = f * sp.exp(-f * rho2)

# Linearly polarized Gaussian seed.  At order zero Bamp is the peak magnetic
# amplitude and reproduces the original B_Gauss_lin field.
A0 = (
    sp.Integer(0),
    I * Bamp * psi0 / k,
    sp.Integer(0),
)

# In specified_orders mode, replace None by complete vector coefficients:
# A2 = (A2x, A2y, A2z)
# A4 = (A4x, A4y, A4z)
# A_orders = {2: A2, 4: A4}
A_orders = None

# This is a field-amplitude FWHM, matching the original input convention.
pulse = sp.exp(-(tau**2) / sigma)
pulse_sigma = (0.5 * fwhm) ** 2 / np.log(2.0)

zR = k0 * waist**2 / 2.0
eps0 = 2.0 / (k0 * waist)
derivative_scales = (
    k * eps / 2,
    k * eps / 2,
    k * eps**2 / 2,
)

# Every MPI process needs the compiled callable.  Expensive symbolic checks and
# file output are performed only by rank zero.
smilei_rank = globals().get("smilei_mpi_rank", 0)
report_path = "lax_symbolic_report.log" if smilei_rank == 0 else None

try:
    lax_field = build_lax_field(
        A0=A0,
        order_coefficients=A_orders,
        construction=lax_construction,
        pulse=pulse,
        pulse_application=pulse_application,
        coords=(X, Y, Z),
        retarded_time=tau,
        eps=eps,
        k=k,
        order=lax_order,
        derivative_scales=derivative_scales,
        group_velocity=group_velocity,
        parameter_subs={
            eps: eps0,
            k: k0,
            Bamp: a0,
            sigma: pulse_sigma,
        },
        validate=(smilei_rank == 0),
        report_path=report_path,
    )
except Exception as error:
    # Physics-validation failures already contain the full backend report.
    # This fallback ensures that structural or import-time failures are also
    # recorded beside the Smilei output before the namelist stops.
    if smilei_rank == 0 and not isinstance(error, LaxValidationError):
        Path("lax_symbolic_report.log").write_text(
            "LAX SYMBOLIC INPUT REPORT\n"
            + "=" * 72
            + "\nStatus: FAIL\n\n"
            + traceback.format_exc(),
            encoding="utf-8",
        )
    raise


# -----------------------------------------------------------------------------
# Vectorised local field evaluator
# Local x[0] is the propagation coordinate.
# -----------------------------------------------------------------------------
def Bfield(x, t):
    x = np.asarray(x)
    local_time = t - t_offset

    Xn = x[1] / waist
    Yn = x[2] / waist
    Zn = x[0] / zR
    taun = local_time - x[0] / group_velocity

    B_backend = evaluate_compiled_vector(
        lax_field.compiled_field, Xn, Yn, Zn, taun
    )

    # Backend (X,Y,Z) -> local (x1,x2,x0).
    B_local_envelope = np.stack(
        [
            B_backend[2],
            B_backend[0],
            B_backend[1],
        ],
        axis=0,
    )

    carrier = np.exp(1j * (k0 * x[0] - omega * local_time))
    return np.real(B_local_envelope * carrier)


# -----------------------------------------------------------------------------
# Tested coordinate and vector transformation from 3d_corner.py
# -----------------------------------------------------------------------------
def RotM(axis, angle):
    if axis == "x":
        return np.asarray(
            [
                [1, 0, 0],
                [0, np.cos(angle), -np.sin(angle)],
                [0, np.sin(angle), np.cos(angle)],
            ]
        )
    if axis == "y":
        return np.asarray(
            [
                [np.cos(angle), 0, np.sin(angle)],
                [0, 1, 0],
                [-np.sin(angle), 0, np.cos(angle)],
            ]
        )
    if axis == "z":
        return np.asarray(
            [
                [np.cos(angle), -np.sin(angle), 0],
                [np.sin(angle), np.cos(angle), 0],
                [0, 0, 1],
            ]
        )
    raise ValueError("axis must be 'x', 'y', or 'z'.")


def transform_vector_field(
    R, offset, vector, selection=(0, 1, 2), vectorised=True
):
    if not vectorised:
        def vector_field_transformed(x_, t_):
            return np.asarray(
                [
                    (R.T @ vector(R @ (x_ - offset), t_))[index]
                    for index in selection
                ]
            )
    else:
        def vector_field_transformed(x_, t_):
            off = np.asarray(offset).reshape((3,) + (1,) * (x_.ndim - 1))
            x_rot = np.einsum("ij,j...->i...", R, x_ - off)
            vec = vector(x_rot, t_)
            vec_rot = np.einsum("ij,j...->i...", R.T, vec)
            return np.asarray([vec_rot[index] for index in selection])

    return vector_field_transformed


Main(
    geometry="3Dcartesian",
    interpolation_order=2,
    cell_length=[l0 / resx, l0 / resx, l0 / resx],
    grid_length=Lsim,
    number_of_patches=[4, 4, 4],
    timestep=t0 / rest,
    simulation_time=Tsim,
    EM_boundary_conditions=[["PML"]],
)


# -----------------------------------------------------------------------------
# Boundary profiles.  The selected components follow Smilei's Laser syntax.
# -----------------------------------------------------------------------------
Rbeam = RotM("z", ang[1]) @ RotM("y", ang[0])


def Bfield_xplane(xfix_):
    def B1(y_, z_, t_):
        return transform_vector_field(
            Rbeam, focus, Bfield, selection=(1,)
        )(np.asarray([xfix_, y_, z_]), t_)

    def B2(y_, z_, t_):
        return transform_vector_field(
            Rbeam, focus, Bfield, selection=(2,)
        )(np.asarray([xfix_, y_, z_]), t_)

    return [B1, B2]


Laser(box_side="xmin", space_time_profile=Bfield_xplane(0.0))


def Bfield_yplane(yfix_):
    def B1(x_, z_, t_):
        return transform_vector_field(
            Rbeam, focus, Bfield, selection=(0,)
        )(np.asarray([x_, yfix_, z_]), t_)

    def B2(x_, z_, t_):
        return transform_vector_field(
            Rbeam, focus, Bfield, selection=(2,)
        )(np.asarray([x_, yfix_, z_]), t_)

    return [B1, B2]


Laser(box_side="ymin", space_time_profile=Bfield_yplane(0.0))


def Bfield_zplane(zfix_):
    def B1(x_, y_, t_):
        return transform_vector_field(
            Rbeam, focus, Bfield, selection=(0,)
        )(np.asarray([x_, y_, zfix_]), t_)

    def B2(x_, y_, t_):
        return transform_vector_field(
            Rbeam, focus, Bfield, selection=(1,)
        )(np.asarray([x_, y_, zfix_]), t_)

    return [B1, B2]


Laser(box_side="zmin", space_time_profile=Bfield_zplane(0.0))


globalEvery = int(rest)

DiagScalar(every=globalEvery)

DiagFields(
    every=globalEvery,
    fields=["Ex", "Ey", "Ez", "Bx", "By", "Bz"],
)
