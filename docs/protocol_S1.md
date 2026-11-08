# Protocol S1: energy stability of learned closures in a Galerkin truncation (side track, frozen 2026-11-08)

Question. For a 2D Kolmogorov flow truncated to a Fourier-Galerkin system, does a learned closure keep the energy from growing? This is a physics check of one law. It makes no claim about blow-up or about the Navier-Stokes Millennium problem.

## Physics fixed in advance
- Domain [0, 2pi]^2, incompressible 2D Navier-Stokes in vorticity form, viscosity nu = 0.02, forcing F = sin(4y).
- Pseudo-spectral DNS, 64 x 64 grid, 2/3 dealiasing, RK4, dt = 1e-3.
- Galerkin truncation keeps all real Fourier modes with |k_x|, |k_y| <= K, with K in {1, 2, 3}. The nonlinear advection term is projected, so it does no net work on the truncated energy (checked numerically in step S1.2).

## Steps and gates
- **S1.1 DNS validation.** Energy budget dE/dt = P - D must match, with relative residual below 1e-3 over the run. Gate S-G1: pass or stop.
- **S1.2 Advection energy neutrality.** The projected advection term must have zero net energy rate (relative error below 1e-10). Gate S-G2: pass or stop.
- **S1.3 Closure.** Each truncation level gets a learned closure: linear in the modes (eddy-viscosity type), a(t) -> closure = C a, with C fitted by least squares to the discarded-mode energy flux from DNS.
- **S1.4 Exact energy-stability test.** The energy rate from the closure is a quadratic form a^T Q a with Q = (C + C^T)/2. Test Q negative semidefinite with exact rational LDL in `fractions.Fraction`, and cross-check against the eigenvalues. Gate S-G3: the two tests agree on every level.

## Outcome to report, whatever it is
For each K, whether the fitted closure is energy-stable, and the smallest eigenvalue of Q. If every closure is stable, report that. The side track is a check, not a claim of novelty.

## Limits written down before running
- Two-dimensional flow. No vortex stretching, so this does not stand in for 3D.
- Linear closures only in this protocol. Nonlinear closures need cubic-form tests and are a later protocol.
- No statement about finite-time singularities.
