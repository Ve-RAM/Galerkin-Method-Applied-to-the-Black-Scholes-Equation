# Solving the Black-Scholes Equation with the Galerkin Method and Hermite Polynomials

Research project (*Práctica Investigativa II*) – Universidad EAFIT

**Authors:** Santiago Vera Ramírez · Fredy Hernán Marín Sánchez

This repository contains a MATLAB implementation of a **spectral Galerkin method based on
Hermite polynomials** for pricing European call options through the Black-Scholes PDE.
The method is compared against:

1. The **analytical Black-Scholes formula** (reference solution).
2. A **finite element method (FEM)** with linear hat functions and a **Crank-Nicolson** time scheme.
3. A **Monte Carlo** simulation with antithetic variates.

## Research hypothesis

> The Galerkin method with Hermite polynomials produces approximations of the Black-Scholes
> solution with a relative error below **0.5 %** with respect to the analytical solution and
> other numerical methods.

This work is also a validation step: the long-term goal is to apply the method to
Black-Scholes variants and other SDE-derived PDEs that have no closed-form solution.

## Methodology

### Problem

Black-Scholes PDE for a European call:

$$\frac{\partial V}{\partial t} + rS\frac{\partial V}{\partial S} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} - rV = 0, \qquad V(S,T)=\max(S-K,0)$$

### 1. Galerkin–Hermite (spectral)
- Change of variables `x = ln(S)`, `τ = T − t`, which maps `(0, ∞)` to ℝ.
- Expansion `u_M(x,τ) = Σ c_k(τ) H_k(x)` using physicists' Hermite polynomials, orthogonal in `L²(ℝ, e^{-x²})`.
- The Galerkin condition yields the ODE system `ċ = −D c`, with `D` a banded matrix (bandwidth 3).
- The system is solved **exactly** with the matrix exponential: `c(T) = expm(−D·T) c(0)`.
- Initial coefficients are computed with adaptive quadrature on `[−10, 10]`.

### 2. Finite elements + Crank-Nicolson
- Truncated domain `[0, 3K]`, `N` uniform elements, linear hat basis functions.
- Tridiagonal mass and stiffness matrices; unconditionally stable, second-order scheme.

### 3. Monte Carlo
- Risk-neutral GBM, exact terminal-price sampling, antithetic variates, and a fixed seed for reproducibility.

### Base parameters

| Parameter | Value |
|-----------|-------|
| r         | 0.05  |
| σ         | 0.20  |
| K         | 100   |
| T         | 1 year |
| S         | 100 (at-the-money) |
| M (Galerkin) | 5, 10, …, 200 |

## Main results

At-the-money comparison (`S = K = 100`):

| Method                              | Price       | Absolute error | Relative error (%) |
|-------------------------------------|-------------|----------------|--------------------|
| Black-Scholes (exact)               | 10.45058357 | 0              | 0                  |
| Galerkin (Hermite, M = 30)          | 10.45058360 | 3.12e-8        | 2.99e-7            |
| Crank-Nicolson (2000 × 2000 nodes)  | 10.45058357 | 1.08e-9        | 1.03e-8            |
| Monte Carlo (500 000 simulations)   | 10.45123456 | 6.51e-4        | 6.23e-3            |

### Key findings

- **Spectral convergence:** the error decays like `e^{-cM}` for `10 ≲ M ≲ 35` and saturates near
  `1e-14 – 1e-15` (round-off limit). With only **31 degrees of freedom**, Galerkin reaches an accuracy
  comparable to Crank-Nicolson with 4 million mesh nodes and far better than Monte Carlo.
- **Robustness to K and r:** relative error stays below 0.31 % for `K ∈ [50, 150]` and below
  0.014 % for `r ∈ [0, 0.20]`.
- **The key quantity is σ²T:** the payoff kink at `x = ln K` is smoothed by diffusion.
  The 0.5 % hypothesis holds for **σ²T ≳ 0.02** (e.g. `σ ≥ 0.15`, `T ≥ 0.5` years).
- **Failure regime:** for `σ²T ≲ 0.01` the method degrades from spectral to algebraic convergence
  and can produce **negative prices** (the Hermite basis is not positive), a behavior not seen in
  FEM or Monte Carlo.
- **Other limitation:** the basis is centered at `x = 0`, so it is less efficient for options far
  in/out-of-the-money (`ln(S/K) ≫ 0`). The matrix exponential costs `O(M³)`.

## Repository structure

> Adapt this section to the actual layout of your code.

```
.
├── src/
│   ├── galerkin_hermite.m      % Galerkin–Hermite solver (expm)
│   ├── fem_crank_nicolson.m    % FEM + Crank-Nicolson solver
│   ├── monte_carlo.m           % Monte Carlo with antithetic variates
│   └── black_scholes_exact.m   % Closed-form reference solution
├── experiments/
│   ├── sensitivity_analysis.m  % Sweeps over K, r, σ, T
│   ├── convergence_M.m         % Error vs. number of polynomials
│   └── compare_methods.m       % Comparison table
├── figures/
├── paper/                      % PDF report
└── README.md
```

## Requirements

- MATLAB (R2019b or later recommended)
- No additional toolboxes expected beyond core MATLAB (`expm`, `integral`, `normcdf` equivalents may be needed).

## Usage

```matlab
% Example: price an ATM call with the Galerkin–Hermite method
r = 0.05; sigma = 0.20; K = 100; T = 1; S = 100; M = 30;
V = galerkin_hermite(S, K, r, sigma, T, M);
```

## Future work

- Rigorous numerical analysis of the method (convergence, stability, consistency, monotonicity).
- Applying the method to Black-Scholes variants and other PDEs derived from SDEs with no analytical solution.
- Exploring shifted/centered Hermite bases to handle deep in/out-of-the-money options.

## References

- F. Black and M. Scholes, *The Pricing of Options and Corporate Liabilities*, J. Political Economy, 81(3), 1973.
- F. Baustian, K. Filipová, J. Pospíšil, *Solution of option pricing equations using orthogonal polynomial expansion*, Applications of Mathematics, 66(4), 2021.
- P. Wilmott, S. Howison, J. Dewynne, *The Mathematics of Financial Derivatives*, Cambridge University Press, 1995.
- P. Glasserman, *Monte Carlo Methods in Financial Engineering*, Springer, 2004.
- G. Strang and G. Fix, *An Analysis of the Finite Element Method*, Prentice-Hall, 1973.

## License

Add your preferred license (e.g., MIT) here.

## Contact

- Santiago Vera Ramírez – sverar1@eafit.edu.co
- Fredy Hernán Marín Sánchez – fmarinsa@eafit.edu.co
