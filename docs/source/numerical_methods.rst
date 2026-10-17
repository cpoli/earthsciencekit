Numerical methods and sampled results
=====================================

Static APIs apply NumPy broadcasting. Invalid real quantities and singular or
overflowing arithmetic raise ValueError; very small exponentials may underflow
to zero. A 100% statement coverage score does not bound floating-point error.

Initial-value integration
-------------------------

:py:func:`earthsciencekit.integrators.integrate_ode` uses SciPy solve_ivp with
adaptive internal timesteps. The supplied time array chooses output samples,
not fixed solver steps. A nonzero first time is valid: the initial state belongs
to that first sample. Domain simulation wrappers require nonnegative times.
Defaults are rtol=1e-8, atol=1e-10 and RK45. State variables with very different
scales should be nondimensionalized or checked against independent solutions.
The generic wrapper currently accepts scalar tolerances only.

Solver failure and invalid derivative shapes raise ValueError. The wrapper
checks labels before calling user callbacks. Callbacks may still have their own
side effects and should be pure if repeatability is required.

Periodic thermal diffusion
--------------------------

:py:func:`earthsciencekit.plate_tectonics.diffuse_temperature` diagonalizes the
centered periodic discrete Laplacian with an FFT. A grid Fourier mode decays
with eigenvalue -4 kappa sin(pi j/n)^2 / dx^2. Evolution of the discretized
system is exact in time to roundoff; it still approximates a continuous PDE in
space. Periodic boundaries conserve the mean. Compare with the fixed-surface
half-space analytical model only after matching boundary conditions.

Result objects
--------------

:py:class:`earthsciencekit.results.TimeSeries` stores time [s], an (n, m) state
matrix and explicit state names/units. Arrays are copied and marked read-only.
``summary()`` reports extrema and final values. ``plot()`` returns a figure and
one axis per state. ``animate()`` returns a Matplotlib FuncAnimation; keep it
alive for display or save it with a writer such as Pillow. Display/export timing
does not change simulated physical time.
