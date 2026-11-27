# A sharp condition for recovering a law sign with an integrating meter

## Observation model

Let μ be the time-normalized occupation probability measure of a fixed measurement window on a metric state space, with finite first moment. Let R be an isometric involution. The meter reports a scalar Gaussian reading with mean sA + ∫q dμ and known standard deviation σ>0. Here s is ±1, A>0 is calibrated, and q is any real L-Lipschitz nuisance satisfying q(Rz)=-q(z). The law sign is the sign of the nominal component; at fixed points of R the nuisance vanishes. The observation path and μ are known and shared by both hypotheses. Full state trajectories observed as additional evidence, sign-dependent dynamics, uncalibrated even terms, and unknown paths require separate analysis.

## Theorem

Write ν=R#μ for the reflected measure. On Euclidean state space the attainable additive bias is precisely [-B,B], with

    B = (L/2) W₁(μ,ν).

Consequently the optimal worst-case probability of identifying the sign from one reading, allowing randomized decisions, is Φ(max(A-B,0)/σ). In particular a positive uniform advantage over guessing exists exactly when A>B. A required probability p>1/2 is attainable exactly when A≥B+σ Φ⁻¹(p).

## Proof

For odd q, ∫q dμ = (1/2)∫q d(μ−ν). Kantorovich–Rubinstein duality bounds its absolute value by B. Conversely, for any 1-Lipschitz potential f, the antisymmetrization (f−f∘R)/2 is odd and 1-Lipschitz, and its integral under μ equals half the integral of f against μ−ν. Taking an optimal dual potential gives equality. Existence follows also by anchoring a maximizing sequence at the origin, extracting a locally uniformly convergent subsequence, and applying dominated convergence via the finite first moments. Negating and scaling the extremizer fills the interval.

The possible Gaussian means are intervals centered at ±A with radius B. If disjoint, the sign-of-reading test achieves the displayed probability. The closest two means ±(A−B) give a simple symmetric Gaussian testing problem whose equal-prior optimum is that same number, bounding minimax accuracy from above. If the intervals touch or overlap, both signs can generate the identical distribution, and a fair independent guess attains the upper bound 1/2.

## A one-window construction and a measurable robustness bound

For x''=−V'(x) with differentiable even V and locally Lipschitz force, put t=0 at an x=0 crossing. By ODE uniqueness x(t)=−x(−t) and v(t)=v(−t). With R(x,v)=(−x,v), a uniform window centered at zero has μ=R#μ. Thus any integrable position-odd contribution cancels exactly, including nonpolynomial terms. A single forward trajectory logged through the crossing suffices. Even terms survive.

For a possibly asymmetric path z on [-h,h], pairing time t with −t yields an explicit transport plan and the computable upper bound
B ≤ L/(4h) ∫ ||z(t)−Rz(−t)|| dt over [-h,h].
State units and the metric must be fixed before specifying L. The bound can be conservative. In particular, accurate evaluation of discretized paths needs sampling and measurement-error bounds. A physical apparatus must calibrate A, σ, and a defensible L; this theorem does not infer them from the test window.

## Prospective numerical check

Run python exploratory/reflection.py. Each odd power term is scaled to have unit integrated shift on its own conservative forward window for that oscillator. The centered window yielded zero at printed precision for all six odd term/oscillator combinations (harmonic and Duffing). The even control retained 0.291 and 0.303 of its respective forward signal. A timing offset of 0.025 on a window of duration 0.5 gave absolute leakage 0.059–0.123 across the odd terms. Centered windows at damping 0.2 gave 0.021–0.071, relative to the conservative forward reference. These simulations test cancellation and sensitivity, not statistical calibration or attainment of the transport extremizer.

The theorem combines classical optimal-transport duality, symmetry and Gaussian testing. No priority claim is made for those ingredients or their combination. The proposed scientific test is prospective prediction of sign recoverability from the measured trajectory imbalance and instrument calibration.
