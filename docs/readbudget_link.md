# Read-budget scaling (exploratory branch) and the fixed-library horizon

Question from the review of `exploratory/`: does the margin scaling in astra's read-budget theorems (m^3 or m^4) say anything about measurement design in this project (0F-4)? Short answer: no. The two results are about different model classes, and in the class this project uses the scaling is m^2. Numbers: `tools/readbudget_scaling.py`, output `results/readbudget_scaling.txt`, figure `figures/fig_readbudget.png`. These are exact formulas, not a simulation.

## The three formulas

To decide a law sign at one point with probability p = 0.9:

| model | measurements needed | slope in the margin | halving the margin costs |
|---|---|---|---|
| fixed library of degree d (Theorem 5, `docs/law_horizon.md`) | n >= sigma^2 Kbar(u) (sqrt q + z_p)^2 / mu0^2 | -2 | x4 |
| slope-bounded class, readout noise (`exploratory/local_theorem.md`) | n >= (4 L sigma z_p / m^2)^2 | -4 | x16 |
| slope-bounded class, accumulating noise (`exploratory/mixed_theorem.md`) | runs >= 3 L eta^2 z_p^2 / m^3 | -3 | x8 |

The fitted slopes in the output file are -2.000, -4.000 and -3.000, as they must be.

## Why they differ

In astra's setting the unknown is any function with |p'| <= L. The worst case puts all the information about p(0) in a window of length h = m / L next to the target and nothing elsewhere. Halving the margin halves the window, so each reading carries less signal (the signal in an integrated reading is about m h, which goes like m^2). That is where the extra powers come from.

In this project the model class is a fixed polynomial library. The law value at a point u is one fixed linear functional of finitely many coefficients, and every data row informs it. Nothing shrinks when the margin shrinks. The only margin dependence is in z = mu0 / sqrt K, so the count goes like mu0^-2 at a fixed point. Corollary 6 then puts the dependence on *where* the point is into K(u), through T_d(g).

The two settings meet only if the library is allowed to grow with the data (a nonparametric fit), so that the degree has to rise as the margin falls to keep the bias below it. That is the optimal-recovery setting astra cites (Donoho 1994; Armstrong and Kolesar 2020). This project does not fit that way, and I do not test it.

## What it means for the project

- For 0F-4 the read-budget theorems predict nothing. 0F-4 asked where to start new trajectories, which changes K(u). The margin scaling is about how much any design costs, and in a fixed library that is mu0^-2 for every design arm alike.
- One sentence for the write-up: in a fixed library, halving a law's margin at a given state costs four times the data. A slope-bounded law class would cost 8 to 16 times, which is one reason to prefer a finite library when one can be justified.
- Nothing here was run as a gate, and no prediction was made.
