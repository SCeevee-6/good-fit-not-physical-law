# When one run cannot resolve a physical sign

## Meter

Let T≥(1+1/√2)m/L, h=m/L. The underlying accumulated reading on [-T,0] is Y(t)=∫p+ηW(t+T) with standard Brownian W and integral starting at -T. Unknown deterministic power is L-Lipschitz with p(0)≥m or p(0)≤−m. Querying Y adds independent N(0,σ²) noise per query. Repeated timestamps are independent readouts of the same analog value, sharing the Brownian path. Both noise scales are known. Accuracy is worst case over this class and both signs.

## Exact single-run ceiling

Even observing the entire continuous Y path without readout noise gives optimal accuracy exactly Φ(sqrt(m³/(3L))/η). Write r(t)=(m+Lt)₊. The admissible drifts ±r are Gaussian shifts in Brownian motion with optimal equal-prior test sign(∫r dY); its accuracy is Φ(||r||₂/η). This can be obtained by Gaussian likelihood on increments and taking the mesh to zero (deterministic square-integrable drift). The very same statistic has at least this accuracy for all positive and negative drifts, by the Lipschitz envelope on the support of r. Thus the pair is least favorable. A finite noisy experiment is a randomization of this full-path experiment. For K truly independent physical runs of the same p with continuous observation, the same proof multiplies the squared signal-to-noise ratio by K. Halving m therefore costs eight independent runs in this ideal drift-limited experiment, versus sixteen independent readouts in the readout-only theorem. The number of runs is rounded up.

## Mixed-channel bound against arbitrary adaptive choices

For n=2k, set λ=4σ²/(nη²) and choose c in (0,m) from (m−c)²/(2L)=λc. Let g=min(r,c), Q=(m−c)²/(2L), and D=||g||₂²+λc². Then every n-read single-run policy has worst-case accuracy at most Φ(U), U=√D/η.

Proof. Construct a positive drift p equal to r on [-h,0], zero before it, except for a negative linear triangular pulse starting at -T and ending at zero after (m−c)/(√2 L) time units, with slope L. The pulse has area −Q/2. The cumulative integral F of p−g therefore stays between −Q/2 and +Q/2, since the late positive residual adds Q. The opposing hypothesis is −p. Give the experimenter an entire Brownian path with drift ±g and n independent scalar readings of ±Q/2 with variance σ². These more informative data simulate the original experiment: at each chosen time add to the drifted path the next scalar multiplied by 2F(t)/Q, plus fresh Gaussian noise to restore variance σ². This works conditional on previous choices. The independent Gaussian channels have total squared separation from the midpoint ||g||₂²/η²+nQ²/(4σ²)=U². Their optimal equal-prior accuracy bounds the original minimax accuracy. The clipping equation minimizes this bound within the displayed split family.

## A schedule with a uniform analytic performance guarantee

Take k early reads at -h+(j−1/2)c/(Lk), j=1,...,k, and k terminal reads at zero. Use the terminal average minus the early average. Define x=c³/(Lk²D). The exact worst-case signal-to-noise ratio achieved is

    Z = U (1+x/24)/sqrt(1+x/6),    0≤x≤3/k².

Proof. Scale the statistic by c. Its stochastic-integral weight K rises from 0 in steps of c/k at the early reads and ends at c. Its independent read variance is η²λc². The minimal signed mean is ∫Kr, achieved by p=r. The shared Brownian variance is η²∫K². For the ideal g, ∫gr=D because the overlap above the clip contributes cQ=λc². Relative to g, rounding the linear rise to nearest levels changes ∫Kr by J/24 and ∫K² by J/6, where J=c³/(Lk²). This follows by integrating on each of the k equal cells: in rescaled cell coordinates g increases linearly from 0 to 1 plus its starting level, while K jumps from 0 to 1 at the midpoint. The constant-level cross terms vanish. Substitution proves the formula. D≥c³/(3L) implies the uniform bound on x. The ratio decreases on [0,3]. With eight reads the achieved Gaussian z-score is at least 99.24% of the upper bound for every positive parameter setting (more precisely 0.9924 rounded down). The bound is finite-sample; the exact mixed-noise optimum is bracketed rather than identified. As n tends to infinity, c tends to m and both bounds converge to the full-path ceiling.

## Frozen test

Run OPENBLAS_NUM_THREADS=1 python exploratory/mixed_resolution.py. The raw covariance computation does not use the rounded-weight identity. Frozen conditions and thresholds are in mixed_protocol.md. All nine n≥8 cells passed; the maximum absolute probability gap there was 0.119 percentage points. The uniform-time comparator uses best nonnegative increment weights on a grid including the start and end. The two-time comparator is the exact readout-only prescription, not a reoptimized mixed-noise two-time schedule. The reported 100000-draw check is one stochastic check of a single extremal drift. No apparatus data, unknown noise calibration, colored departures from Brownian motion, finite-width readout windows, or instrument latency were tested. The ceiling is a one-run conditional guarantee, not a bound on all physical instruments.

## Scientific positioning

Gaussian white-noise optimal recovery and triangular kernels are established. Appropriate comparisons include Armstrong and Kolesár, Simple and honest confidence intervals in nonparametric regression (Quantitative Economics, 2020), and Donoho, Statistical estimation and optimal recovery (Annals of Statistics, 1994). The candidate contribution to check for priority is the finite-read mixed-channel adaptive upper bound together with its constructive near-attainment ratio. No priority claim is established. A prospective physical test should independently separate shared accumulation noise from per-read error, freeze parameters and acquisition limits, and compare the forecast for rereads versus independent resets.
