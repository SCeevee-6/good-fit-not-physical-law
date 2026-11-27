# Can an experiment preserve a target sign when the dynamics are wrong?

Run from the repository root with NumPy and SciPy. The frozen reset protocol is in crossover_protocol.md. Commands:

    OPENBLAS_NUM_THREADS=1 python exploratory/crossover.py
    OPENBLAS_NUM_THREADS=1 python exploratory/phase_control.py

## Frozen test

Training uses a 150-member coefficient grid at assumed coupling 0.3. Evaluation draws 1,500 fresh continuous parameter vectors; 384 positive and 559 negative qualify at the prespecified margin. Readings are integrated using a separate adaptive ODE solver, cross-checked against RK4 for the original form. The metric is the smaller of the two sign-specific mean probabilities of a correct decision, integrating Gaussian instrument noise analytically. Worst individual probabilities are provided separately. These are finite-bank conditional probabilities, not experimental confidence bounds.

The predicted crossover r*=43.45 is outside [0,1]. No interior-switch transfer test was possible. Increasing reset correlation alone left repeats ahead in worst-sign mean throughout both unseen couplings with the original gain form. Adding +0.25 xv to gain reversed the ordering already at zero correlation at both couplings. The extra term is zero at x=0, the target state, but it affects the path followed by an integrating meter. A wrong transfer ranking can thus arise even when the perturbation contributes exactly zero to the target expression.

## Post hoc same-budget control

Use four independent half-unit runs at s=1.5, comparing four identical starts with two replicates at each of two distinct starts. Both designs select starts by maximizing the minimum signed projected development margin over permitted starts. Both use the same centroid-direction classifier. The comparison isolates schedule choice within this class of detectors. Select designs before evaluating a fresh coefficient bank (2,500 draws before margin filtering). The nominal design has no omitted term; a sensitivity arm allows xv coefficients -0.25, 0 and +0.25 at design time. Test both signs of 0.35 at couplings 0.1 and 0.5. All functional-form analyses here are exploratory. Phase indices parameterize (x,v)=sqrt(1.5)(sin p,cos p), p=index*pi/8.

Lowest worst-sign mean over both test couplings and perturbation signs for the nominal design:

| Extra gain term | Four at phase 7 | Two each at phases 0 and 7 |
|---|---:|---:|
| xv | 71.63% | 99.93% |
| x^4 | 99.99% | 99.97% |
| (x²+v²-1.5)² | 99.98% | 99.95% |

The xv benefit is specific. For that perturbation the worst individual pair accuracy across test points was 90.79%, below the 96.34% nominal grid floor. Even this improved schedule does not certify unseen dynamics. The other two perturbations barely challenged the single phase at these coefficients; equal coefficients do not imply equal integrated perturbation sizes. The augmented single detector can perform below chance on individual models because projected classes overlap; raw signed-threshold predictions are reported, without a random-guess fallback.

## Mechanistic question

For fixed trajectories the integrated effect of the xv term is proportional to the integral of x v³. A weighted combination of two starting phases could attenuate this term while preserving sensitivity to the target. With back-action the path itself changes, so this argument motivates a test rather than proving cancellation. Next measure the derivative of the combined readout with respect to each nuisance coefficient, match perturbations by induced meter shift, and test transfer to an independently specified oscillator. Reset correlation models shared readout error, not errors in actual initial state. All follow-up phase comparisons used independent errors.
