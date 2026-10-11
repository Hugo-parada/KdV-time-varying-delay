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

`MatAz.m` imposes `z(x,0,t_{n+1})=y(x,t_{n+1})` in the same coupled solve. 
