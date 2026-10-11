# KdV equation with time-varying delay

MATLAB codes accompanying the numerical simulations in:

> H. Parada, C. Timmimoun, and J. Valein, **Stability results for the KdV equation with time-varying delay**, *Systems & Control Letters* **177** (2023), 105547.  
> DOI: https://doi.org/10.1016/j.sysconle.2023.105547  
> HAL: https://hal.science/hal-03819356v2

The repository contains corrected and cleaned versions of the MATLAB codes recovered from the original 2022 simulation archive.

## Contents

- `boundary-delay/`: nonlinear KdV equation with time-varying or constant delayed boundary feedback.
- `internal-delay/`: nonlinear KdV equation with time-varying or constant delayed internal feedback.
- `smoke_test.m`: short consistency test for both solvers.

## Numerical method

For the KdV part, the space interval is discretized uniformly. The linear operator

```text
D_x + D_x^+ D_x^+ D_x^-
```

is treated implicitly, while the nonlinear term is evaluated explicitly as

```text
y_i^n D_x^+ y_i^n.
```

A time-varying delay is represented by a transport variable. In the boundary case,

```text
z(rho,t) = y_x(0,t-tau(t)*rho),   0 < rho < 1,
```

and in the internal case,

```text
z(x,rho,t) = y(x,t-tau(t)*rho).
```

In both cases the transformed delay equation is

```text
tau(t) z_t + (1-dtau(t)*rho) z_rho = 0.
```

It is discretized implicitly in time and by an upwind finite difference in `rho`. The transport equation and the KdV equation are solved in one coupled linear system at every time step, so the boundary condition at `rho=0` is imposed at the same new time level.

The code computes the discrete counterparts of the energies used in the paper, including the square of the delay state, the quadrature factors, the factor `1/2`, the time-dependent factor `tau(t)`, and, in the internal case, the weight `xi(x)` only on the feedback region.

## Boundary time-varying delay

The main solver is

```matlab
[t,E] = kdv_time_delay(alpha,beta,d);
```

with

```text
tau(t) = d*(1.5+sin(t)).
```

The initial data used in the paper are

```text
y0(x)   = 0.5*(1-cos(2*pi*x)),
z0(rho) = -0.5*sin(2*pi*rho).
```

Main scripts:

- `boundary-delay/figure_delay_size.m`: influence of `d`, with `alpha=beta=0.1` and `d=0.1,0.7,1.3`.
- `boundary-delay/figure_constant_vs_variable.m`: comparison of `tau(t)` with `tau_max=2.5d` and `tau_min=0.5d`.
- `boundary-delay/figure_feedback_amplitude.m`: additional experiment from the recovered archive varying `beta`.

The run with `d=1.3` is intentionally outside the hypothesis `sup_t dtau(t)<1`; the solver issues a warning when the transformed transport speed changes sign.

For a constant boundary delay, use

```matlab
[t,E] = kdv_cte_delay(alpha,beta,h);
```

## Internal time-varying delay

The numerical experiment uses

```text
a(x)=2,  b(x)=1,  xi(x)=2.1   on (0,L/2),
```

and zero feedbacks outside this interval. The initial data are

```text
y0(x)       = 1-cos(2*pi*x),
z0(x,rho)   = (1-cos(2*pi*x))*cos(2*pi*rho).
```

The recovered archive called the three simulations with central delay values `K=1,1.5,2` but labeled the curves `M=1.5,2,2.5`. To remove this ambiguity, the corrected interface uses `M` as the maximum delay:

```matlab
[t,E] = kdvint_timedelay(M);
```

with

```text
tau(t) = (M-0.5)+0.5*sin(t).
```

Thus `M=1.5,2,2.5` gives exactly the three delay functions used by the recovered simulations, while the parameter now agrees with the figure labels.

Main scripts:

- `internal-delay/figure_delay_size.m`: `M=1.5,2,2.5`.
- `internal-delay/figure_constant_vs_variable.m`: comparison of `tau(t)=1+0.5*sin(t)` with the constant delays `1.5` and `0.5`.

For a constant internal delay, use

```matlab
[t,E] = kdvint_ctedelay(h);
```

## Corrections with respect to the recovered archive

The recovered codes contained several inconsistencies that affected the numerical scheme or the plotted energy. The corrected version:

- uses the prescribed transformed initial histories directly, instead of evaluating them again at `-tau(0)*rho` or multiplying `rho` by the delay;
- evaluates the nonlinear term as `y_i D_x^+ y_i`, including the missing `1/dx` factor;
- fixes the first-row coefficient of the boundary transport matrix (the old code accidentally used the final loop index `k=Nrho` instead of `k=1`);
- evaluates `tau` and `dtau` consistently at the new implicit time level;
- uses the correct squared delay norm and the full quadrature factors in the boundary energy;
- uses the same physical energy for constant- and variable-delay comparisons;
- fixes the `d=0.1` boundary experiment, which was previously run with `d=0.2` but labeled `d=0.1`;
- aligns the internal delayed state with the correct spatial grid point, removing the one-cell shift in the old code;
- imposes `z(x,0,t_{n+1})=y(x,t_{n+1})` in the coupled solve instead of injecting the previous KdV state;
- uses `xi=2.1` and integrates the internal memory energy only on the feedback support;
- removes the spurious factor `2` that appeared in the constant internal-delay energy after the first time step;
- resolves the mismatch between the internal parameter values and their figure labels.

## Running the code

Run each figure script from its own directory, or add both code directories to the MATLAB path. For a quick check from the repository root, run

```matlab
smoke_test
```

The grid parameters can be changed through an optional structure. For example,

```matlab
opts = struct('Nx',200,'Nrho',150,'Nt',500,'T',10);
[t,E] = kdv_time_delay(0.1,0.1,0.2,opts);
```
