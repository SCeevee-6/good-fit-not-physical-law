# A sharp read-budget theorem for the sign of instantaneous power

## Instrument and uncertainty class

An integrating meter starts at zero at time -T. A query at time t in [-T,0] returns F(t)+independent N(0,σ²) noise, where F(t)=∫ from -T to t of p(u) du. Power p may be any L-Lipschitz time function, with p(0)≥m or p(0)≤−m. Here L,m,σ>0. Exactly repeated timestamps represent independent readouts of the same stored analog value. The result concerns this instrument; direct gated integration or reuse of an already noisy digital value changes the experiment.

## Exact read-budget theorem

Put h=m/L and Q=m²/(2L). Suppose T≥(1+1/√2)h. For every positive even number n of readouts, the largest worst-case probability of identifying the sign is

    Φ(√n Q/(2σ)) = Φ(√n m²/(4Lσ)).

This holds even if later timestamps depend on earlier observations. An optimal policy takes n/2 independent readings at -h and n/2 at 0 and uses the sign of the difference of the two averages. If timestamps must be strictly distinct the same value is a supremum, approached by tight clusters (bias in differences is controlled for sign-separated p via positive normalized averaging kernels). The two-read result requires no repeated timestamps.

### Achievability

For positive sign, p(t)≥m+Lt on [-h,0]. Hence F(0)-F(-h)≥Q. For negative sign the reverse inequality gives ≤−Q. The averaged difference has noise SD 2σ/√n. Thresholding at zero yields the asserted accuracy.

### Universal converse, including adaptive schedules

Construct one positive-sign history. In the first h/√2 units set p(-T+s)=−L(h/√2−s), then keep p=0 until time -h. On [-h,0] set p(t)=L(t+h). The time assumption makes the pieces disjoint. The first triangle accumulates −Q/2, the last adds Q, and |F(t)|≤Q/2 at every instant. Negating p gives an admissible opposing law. Both paths saturate the decision bound at the selected timestamps.

Even a designer who knows this pair cannot beat n independent measurements with means ±Q/2 and variance σ². To see this, take the ith such measurement X and an independent standard normal Z. At the timestamp selected from prior simulated observations set c=2F(t)/Q. Return cX+σ√(1−c²)Z. Its conditional distribution under either sign is exactly the required next meter reading. This simulation works recursively for every adaptive randomized schedule. Equal-prior Gaussian testing on the more informative experiment has accuracy Φ(√n Q/(2σ)), hence this upper bounds worst-case accuracy. The simple policy attains it.

### Each reading alone can be ambiguous

In the numerical example T=2≥2h. A constant zero integral at either selected timestamp t is achievable with positive p(0)=m, for example p(u)=m+ku with k=2m/(T−t)≤L. Its negative gives the same noiseless result with the opposite sign. These histories can differ between timestamps. The joint result works because a single admissible history must explain both readings.

## Numerical check and baselines

Run OPENBLAS_NUM_THREADS=1 python exploratory/local_power.py. It uses a different route: a feasible bounded piecewise-constant derivative yields a least-favorable candidate and an upper bound; a projected meter statistic and exact integral of its absolute sensitivity kernel give a lower bound over the entire Lipschitz class. Floating-point brackets are numerical checks rather than formal interval arithmetic. Parameters and predictions were committed before evaluation in local_protocol.md.

At m=0.35, L=1, T=2, σ=0.0125 the best two-read spacing in the tested grid was 0.35, independently matching h. Either selected reading alone gives minimax chance. The pair achieves 99.9735%; four equally spaced reads give 99.7661%. When σ=0.05 the corresponding figures are 80.6811% and 76.0250%. A short gate reset at -h and read once has SD σ rather than √2 σ and gives 99.99995% in the low-noise case: a stronger instrument wins. Four equally spaced readings and four reads with two early, low-value timestamps are baseline schedules, not optimal four-read designs.

## Scientific scope

The target is local power sign, allowing unknown even and odd effects and arbitrary histories subject to the slew-rate bound. This extension needs no reflection symmetry. It does need a defensible upper bound on the total rate of change, independently calibrated readout noise and a declared minimum resolvable margin. The least-favorable histories belong to this broad class and may be excluded by a specific oscillator equation. The upper bound need not apply to such a narrower class. Numerical differentiation, optimal recovery and Gaussian convex testing have extensive prior art; a priority claim would require a separate literature review. The physical test is to freeze independently measured m, L and noise assumptions and predict spacing, attainable accuracy and an adversarial waveform before collecting confirmatory traces.
