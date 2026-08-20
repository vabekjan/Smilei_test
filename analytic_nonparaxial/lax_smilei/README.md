# Lax-series Smilei prototype

This folder is independent of the earlier notebooks. It contains:

- `lax_series.py` — symbolic Lax construction, validation, pulse handling,
  NumPy compilation and reporting;
- `3d_corner_lax.py` — a Smilei namelist with a readable symbolic linearly
  polarised Gaussian input;
- `test_lax_series.py` — backend unit tests.
- `requirements.txt` — Python dependencies used by the symbolic backend.

## User choices

The namelist exposes one expansion order, referring only to the vector
potential:

```python
lax_order = 4
```

The backend then retains the complete magnetic series through order 5 and
records the consequence in `lax_symbolic_report.log`.

Two independent choices are available:

```python
lax_construction = "focus_preserving"   # or "specified_orders"
pulse_application = "magnetic_field"    # or "vector_potential"
```

`focus_preserving` generates the higher Lax coefficients with
`A^(2j)(Z=0)=0`. `specified_orders` uses complete user expressions supplied in
the namelist:

```python
A_orders = {
    2: A2,
    4: A4,
}
```

All supplied coefficients must use the backend convention
`exp(i*k*z-i*omega*t)`. Expressions written with the conjugate carrier must be
complex-conjugated before use.

With `magnetic_field`, the monochromatic magnetic field is multiplied by the
pulse envelope. With `vector_potential`, the envelope multiplies A before the
curl, so the derivative of the envelope is included.

## Report

Rank zero writes `lax_symbolic_report.log` in the Smilei run directory. The
report contains the configuration, symbolic coefficients, retained magnetic
orders, Lax-hierarchy checks, focus checks, divergence result and compilation
smoke tests. Other ranks still compile the evaluator but skip the expensive
checks and file output.

## Test

The Python environment used by Smilei must provide NumPy and SymPy. From this
folder, run:

```bash
python -m unittest -v test_lax_series.py
```

For the first Smilei comparison, keep
`pulse_application="magnetic_field"`. At order zero this follows the legacy
Gaussian convention. Then compare with `vector_potential`, first at the pulse
centre and subsequently on the rising and falling edges.
