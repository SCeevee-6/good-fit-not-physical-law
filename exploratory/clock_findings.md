# Can a clock tell us how small a frequency change it can resolve?

A clock records phase; instantaneous frequency offset is its derivative. This is a physical analogue of an energy meter recording accumulated power. We used measured backgrounds and numerically injected known positive and negative frequency ramps. The question was whether a decision threshold calibrated early in an instrument record transferred to later measurements. The injections were added offline. The original clocks were never commanded to execute the ramps.

## Reproduction and provenance

Run from the repo root:

    OPENBLAS_NUM_THREADS=1 python exploratory/clock_transfer.py
    OPENBLAS_NUM_THREADS=1 python exploratory/clock_replication.py

Both programs fetch pinned AllanTools data and check the Git blob identity before parsing. Alternatively pass local files as indicated by each script. Python, NumPy and SciPy suffice. The raw data are not redistributed here. Source headers describe the [cesium–maser comparison](https://github.com/aewallin/allantools/tree/32ad8fce6be9f51dda4ab407931f71503c50d101/tests/Cs5071A), [counter cable-delay noise-floor measurement](https://github.com/aewallin/allantools/tree/32ad8fce6be9f51dda4ab407931f71503c50d101/tests/Keysight53230A_ti_noise_floor), and [GPS–maser comparison](https://github.com/aewallin/allantools/tree/32ad8fce6be9f51dda4ab407931f71503c50d101/tests/gps). All are related timing-counter measurements. They provide several noise environments within one modality.

## Chronological separation

The first third estimates the background slope and fixes the detector. The second third supplies a tail estimate for a competitor or a resolution threshold for the replication. The last third supplies every reported test outcome, using complete nonoverlapping blocks of 128 samples. Both injection signs are evaluated on each block. Accuracy here means the smaller of the two empirical sign-specific success rates. Complete-block counts in each final third: 1,450 cesium, 145 cable, 628 GPS. Ordinary sampling intervals are one second (nominal 1PPS cadence for GPS).

The first protocol was frozen before scoring cesium. It required all three model forecasts of 90% to land within three percentage points. **Failed:** 91.31%, 90.90%, 82.48% for 4, 16, 64 seconds. At 64 seconds a direct empirical comparator already forecast 82.00%. This identified a forecasting failure before the more elaborate follow-up. The first run assigns both signs but does not estimate an actual clock's physical slew-rate bound.

## A laboratory interpretation of the threshold

Fix nonnegative weights on consecutive phase increments, normalized so the least-favorable unit terminal ramp has weighted response one. With a slew bound L=m/h and terminal magnitude at least m, the integrated ramp with rate m/h is the lower envelope for a positive sign across the last h seconds; the negative envelope is symmetric. Thus the scalar detector response under the least-favorable waveforms is Z±m, where Z is the weighted background. Other admissible waveforms of each sign only improve this detector's decision, provided the background is additive and invariant to the injected signal. When comparing different chosen m at fixed h the implied slew bound scales too. This is not a calibration of L from the apparatus.

**Proposition (fixed-detector resolution).** For a continuous strictly increasing background CDF with quantiles qα and q₁₋α, the smallest m for which any scalar threshold gives at least 1−α success for each sign is (q₁₋α−qα)/2. An attaining threshold is their midpoint. Proof: positive sign requires the threshold to be at most m+qα; negative sign requires it to be at least q₁₋α−m. These requirements overlap exactly at the stated margin. The statement concerns threshold rules on the chosen scalar response; it makes no optimality claim across all nonlinear tests. Empirical ties require care.

The replication uses conservative finite-sample tail order statistics. Under iid blocks with a continuous CDF, each side has at least 97.5% calibration confidence for the target decile (binomial rank construction), hence both bounds hold jointly with at least 95% confidence for a fixed detector by a union bound. This probability concerns calibration, not achieved sample accuracy on an arbitrary future segment. Dependence, drift and digitization limit its interpretation on actual clock data. The chronological tests are therefore empirical.

## Frozen replication outcomes

Before scoring cable and GPS, the gate required >=87% at every one of their six settings for the original filter after tail calibration. **Passed 6/6.** Four-second values fell below the nominal 90% target on both external records. Amplitudes are allowed to change: this asks how large a signal must be to meet a reliability target. It is not a claim of improved detection at fixed amplitude.

| Record | Window (s) | Gaussian forecast amplitude: observed accuracy | Tail-calibrated amplitude: observed accuracy | Tail/Gaussian amplitude |
|---|---:|---:|---:|---:|
| Cable delay | 4 | 87.59% | 88.97% | 1.267 |
| Cable delay | 16 | 92.41% | 95.17% | 1.362 |
| Cable delay | 64 | 88.97% | 93.79% | 1.349 |
| GPS–maser | 4 | 86.31% | 88.38% | 1.124 |
| GPS–maser | 16 | 87.42% | 92.68% | 1.249 |
| GPS–maser | 64 | 93.15% | 92.04% | 0.915 |
| Cesium–maser (development) | 64 | 82.48% | 90.55% | 1.393 |

The same detector overestimated long-window resolution on cesium and underestimated it on GPS. A universal safety factor cannot explain the direction of both corrections. At the same 64-second cesium injection amplitude the mixed filter recovered 82.48%, ordinary endpoint differencing 65.72%, triangular weighting 61.66%. Empirical covariance fitting offered no consistent additional gain at a smaller margin (see all three methods in clock_replication.csv). The simple tail comparator should remain a baseline.

## What would change the scientific claim

These tests concern fixed-window signal resolution in recorded timing noise. They do not test the theorem's optimal time allocation, impose independent physical resets, or establish Brownian dynamics. They also do not measure hidden natural frequency events. A genuinely prospective experiment would lock margins before acquiring new data, command blinded frequency ramps with an independent source, record the reference sign on a second channel, and test both a fresh run and repeated reads of the same underlying run. Include a continuously held-out instrument or lab and compare at equal amplitudes and equal observation budgets. Preserve negative cells and data gaps.
