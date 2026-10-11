# KdV equation with time-varying delay

MATLAB codes accompanying the numerical simulations in:

> H. Parada, C. Timmimoun, and J. Valein, **Stability results for the KdV equation with time-varying delay**, *Systems & Control Letters* **177** (2023), 105547.  
> DOI: https://doi.org/10.1016/j.sysconle.2023.105547  
> HAL: https://hal.science/hal-03819356v2

## Repository layout

### `boundary-delay/`

Boundary feedback with constant or time-varying delay.

Main files:

- `kdv_time_delay.m`: time-varying boundary delay.
- `kdv_constant_delay.m`: constant boundary delay.
- `mat1.m`: implicit KdV matrix and delayed boundary coupling.
- `mat2.m`: transport matrix for the time-varying delay.
- `mat3.m`: transport matrix for the constant delay.
- `Cases1.m`, `Cases2.m`, `Cases3.m`: numerical experiments from the recovered archive.

### `internal-delay/`

Localized internal feedback with constant or time-varying delay.

Main files:

- `kdvint_timedelay.m`: time-varying internal delay.
- `kdvint_constant_delay.m`: constant internal delay.
- `MatA.m`: implicit KdV matrix with instantaneous damping.
- `MatAz.m`: transport matrix for the time-varying memory variable.
- `MatAz_constant.m`: constant-delay transport matrix.
- `Cases1_int.m`, `Cases3int.m`: experiments from the recovered archive.

The names `mat1`, `mat2`, `MatA`, `MatAz`, `Cases1`, etc. are retained because they stay close to the original numerical work; only non-English abbreviations such as `cte` were replaced.

## Boundary delayed feedback

The model is

```text
y_t+y_x+y_xxx+y*y_x=0,
y(0,t)=y(L,t)=0,
y_x(L,t)=alpha*y_x(0,t)+beta*y_x(0,t-tau(t)).
```

The memory variable is

```text
z(rho,t)=y_x(0,t-tau(t)*rho),
```

and satisfies

```text
tau(t) z_t+(1-tau'(t)*rho) z_rho=0,
z(0,t)=y_x(0,t).
```

For `kdv_time_delay.m`,

```text
tau(t)=d*(1.5+sin(t)).
```

The KdV linear part is treated implicitly with

```text
D_x + D_x^+ D_x^+ D_x^-,
```

and the nonlinear term is explicit:

```text
y_i^n D_x^+ y_i^n.
```

`mat1.m` contains the KdV matrix and the elimination of the right-boundary value using

```text
y_{Nx-1}=-alpha*y_1-beta*dx*z_{Nrho}.
```

`mat2.m` contains backward Euler in time and an upwind difference in `rho`:

```text
(1+q_k)z_k^{n+1}-q_k z_{k-1}^{n+1}=z_k^n,
q_k=dt*(1-tau'(t_{n+1})rho_k)/(tau(t_{n+1})drho).
```

At the first memory point, `z_0^{n+1}=y_1^{n+1}/dx` is imposed directly in the coupled linear system.

The discrete energy follows the energy of the paper,

```text
E(t)=1/2 int_0^L |y|^2 dx
    +|beta| tau(t)/2 int_0^1 |z|^2 d rho.
```

## Internal delayed feedback

The internal memory variable is

```text
z(x,rho,t)=y(x,t-tau(t)*rho),
```

with the same transport equation

```text
tau(t) z_t+(1-tau'(t)*rho)z_rho=0,
z(x,0,t)=y(x,t).
```

The numerical experiment uses

```text
a(x)=2, b(x)=1, xi(x)=2.1 on (0,L/2),
```

and zero outside that region.

The recovered code used central delay values `K=1,1.5,2` while the figure labels were `M=1.5,2,2.5`. To make the parameter agree with the figures, `kdvint_timedelay(M)` uses

```text
tau(t)=(M-0.5)+0.5*sin(t),
```

so `M` is the maximum delay.

`MatAz.m` imposes `z(x,0,t_{n+1})=y(x,t_{n+1})` in the same coupled solve. This removes the one-cell spatial shift and old-time inflow that were present in the recovered implementation.

The memory part of the energy is integrated only on the feedback region and uses `xi=2.1`, as in the analytical model.

## Corrections made to the recovered archive

The corrected code keeps the original style but fixes the numerical inconsistencies found in the 2022 files:

1. The nonlinear term contains the missing factor `1/dx` and is implemented as `y_i D_x^+y_i`.
2. In `mat2.m`, the first transport row uses the coefficient at `rho_1=drho`; the recovered code accidentally used the final loop index `k=Nrho`.
3. `tau` and `tau'` are evaluated consistently at the new implicit time level.
4. The transformed initial histories are evaluated on the actual `rho` grid, without applying the delay a second time.
5. Boundary energy uses the square of the memory norm, the quadrature factors, and the factor `1/2`.
6. Constant- and variable-delay simulations use the same physical definition of energy.
7. `Cases2.m` really uses `d=0.1` for the curve labelled `d=0.1`; the recovered script used `d=0.2` while displaying `0.1`.
8. The internal memory is aligned with the same spatial point as the KdV unknown; the recovered code used the previous spatial cell.
9. The internal condition `z(x,0,t)=y(x,t)` is imposed in the coupled linear system rather than by inserting the previous KdV state.
10. The internal memory energy is restricted to the feedback support and uses `xi=2.1`.
11. The spurious factor `2` that appeared after the first time step in the constant internal-delay energy is removed.
12. The parameter used in the internal-delay experiment now agrees with the label `M` in the figure.

## Running the experiments

From `boundary-delay/`, for example:

```matlab
Cases2
Cases3
```

or directly

```matlab
[t,E]=kdv_time_delay(0.1,0.1,0.2);
[t,E]=kdv_constant_delay(0.1,0.1,0.2);
```

From `internal-delay/`:

```matlab
Cases1_int
Cases3int
```

and for a constant internal delay:

```matlab
[t,E]=kdvint_constant_delay(1.5);
```

A short repository-level check is available with

```matlab
smoke_test
```

The cleanup environment did not contain MATLAB/Octave, so the final MATLAB execution should still be run locally before regenerating publication-quality figures.
