# Worked example: perfect terminal data that cannot decide an energy law

This is the opening example for the write-up. It shows the main point of the project in two components and no statistics: data can fit perfectly and still leave a physical law undecided. It comes from the exploratory branch (`exploratory/hidden_storage.md`). Here I recast it in law-zone terms and recheck it with `tools/check_circuit.py` (output `results/check_circuit.txt`, ALL_OK True).

## The circuit

A series R-L branch in parallel with a series R-C branch, both resistors R, with L = tau R and C = tau / R. Start at rest. Let x be the capacitor voltage and y = R i_L. Both obey

    tau z' = v - z,  z(0) = 0,

so y = x for every drive v(t). The current into the port is i_L + i_C = x / R + (v - x) / R = v / R. That holds for every tau and every voltage history, so no measurement at the terminals, however long or precise, can tell one tau from another. The stored energy is E = L i_L^2 / 2 + C x^2 / 2 = tau x^2 / R, and dE/dt = 2 x (v - x) / R. The check script confirms all four identities symbolically (sympy).

## The law and its zones

Take the law "the circuit is storing energy just after the input steps down", dE/dt >= 0 at t = 1+ for the drive v = 2 V on [0, 1), then 1 V. Then x(1) = 2 (1 - exp(-1/tau)), and

| tau | x(1) | dE/dt at 1+ | port current | port power |
|---|---|---|---|---|
| 0.2 | 1.9865 | -3.9195 W | 1 A | 1 W |
| 3 | 0.5669 | +0.4910 W | 1 A | 1 W |

Both circuits take in exactly 1 W through the port at that instant. The first is releasing stored energy faster than it takes energy in, and the second is storing part of what comes in. The sign changes at tau* = 1 / ln 2 = 1.4427. The check script recomputes both rates by integrating the two branches separately, without using y = x, and they agree to 1e-7.

In the language of the main project, the model class is {tau > 0} and the terminal data are consistent with all of it. The data-consistent set is therefore the whole class. The law is certified for tau > tau* and refuted for tau < tau*. Over the set, the law is undecided, and it stays undecided at any amount of terminal data. One internal measurement, the capacitor voltage at t = 1, settles it.

That is the same structure the project finds in learned dynamical systems. In the circuit no data at all decide the law. In the fitted models, data decide it near where they were taken and stop deciding it farther out (`docs/law_horizon.md`). The design question in both is which extra measurement settles the law.

## Figure

`figures/fig_circuit.png`: (a) port current for tau = 0.2 and tau = 3, identical; (b) stored-energy rate for the same two circuits, opposite in sign after the step; (c) the rate at t = 1+ against tau, crossing zero at tau*.

## Limits

Ideal components. With 1 % capacitor mismatch the port behaviour is no longer exactly equal (astra's check: 0.37 % and 0.24 % of peak current on this drive), so a precise enough terminal measurement could in principle tell the circuits apart. The example is about what exact terminal data can and cannot determine, not about a bench experiment.
