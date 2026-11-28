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

## Addendum 1 (2026-11-28, before the first run)

The 2026-11-08 text left some run details open. They are fixed here before anything runs. Nothing above is changed.

- Forcing. "F = sin(4y)" is read literally as the forcing of the vorticity equation, w_t + u . grad w = nu lap w + F. In velocity form that is f = (cos(4y)/4, 0). Energy input P = <u_x cos(4y)> / 4, dissipation D = nu <w^2>.
- Run. Seed 0. The initial vorticity is random in the modes with |k| <= 8, rescaled to sd 0.1. Spin-up 50 time units, then 100 time units of sampling, with a snapshot every 0.1 (1000 snapshots). RK4 at dt = 1e-3 on 64^2 with 2/3 dealiasing, as frozen. One run.
- S-G1 is computed as |E(end) - E(start) - integral of (P - D)| / integral of (|P| + |D|), with a trapezoid over every step, for the whole run and for every window of one time unit. Both must be below 1e-3.
- S-G2 uses 1000 random resolved states per K. Pass if every |energy rate| / (|a| |rate vector|) is below 1e-10.
- Energy coordinates. Half-plane Fourier modes, a = sqrt2 (Re, Im) c_k / |k|, so the resolved energy is |a|^2 / 2. Closure target per snapshot: tau = -P_K N(w) + P_K N(w_K) with N = u . grad w (the subgrid flux). Linear closure b ~ C a by least squares over the 1000 snapshots; energy rate a' C a, Q = (C + C') / 2.
- Code: `src/ns_dns.py`, `src/ns_galerkin.py`, `src/run_S1.py`, `src/analyze_S1.py`; run on Kaggle (mode S1).

Prediction written before the run. The forcing sits at |k| = 4, outside every truncation K <= 3. The resolved modes therefore get no direct forcing, and in a statistically steady state their only energy source is the closure term. Its true mean energy rate must then be positive, balancing resolved viscosity. A closure that fits well should inherit that, so I expect Q to have positive eigenvalues at every K, i.e. no fitted closure is energy-stable. If so, the energy-stability law is the wrong law for this truncation, because here the subgrid scales are the source. That is the result to report, not a failure of the closure.

Seen before this addendum: a pipeline test with 1 time unit of spin-up and 10 snapshots (outputs deleted). Budget residual 1e-7, advection neutrality 5e-16, all three Q indefinite, mean closure rate positive at K = 2 and 3 and slightly negative at K = 1. With 10 snapshots against up to 48 modes that fit is meaningless (R^2 = 1), so I take nothing from it beyond the pipeline working.
