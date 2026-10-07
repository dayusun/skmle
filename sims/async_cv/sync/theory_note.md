# Bandwidth choice for the bias-corrected kernel GEE (working note, 2026-10-05)

Status: derivations below are exact where marked; the selection result (Section 4) is a
proof sketch with one gap stated at the end. Simulation numbers live in `summary_*.csv`.

## 1. Setup

Subjects i.i.d. Response $Y(T)$ at times $T_{ij}$, covariate $x(S)$ at times $S_{ik}$,
$X(t)=(1,x(t))^\top$, $E[Y(T)\mid x]=\mu\{X(T)^\top\beta_0\}$. Pairs $(T,S)$ within a subject,
lag $u=T-S$, weight $K_h(u)$.

- (A1) $x$ second-order stationary, $\gamma(u)=E\{x(t)-x(t-u)\}^2$, $\gamma(0)=0$.
- (A2) Observation times independent of $(x,\varepsilon)$.
- (A3) $Y(T)\perp x(S)\mid x(T)$.
- (A5, logistic only) $x$ Gaussian.

## 2. Bias of the three estimators

**Plain kernel GEE (identity).** Exact: $E[X(S)X(T)^\top]=\Sigma_X-\tfrac12\gamma(u)e_2e_2^\top$, so

$$E\,U_h(\beta)=\bar\Sigma_h(\beta_0-\beta)-\tfrac12\beta_{0x}\,\bar\gamma_h\,e_2,\qquad
\bar\gamma_h=\frac{E\,K_h(u)\gamma(|u|)}{E\,K_h(u)}.$$

Bias of $\hat\beta_x$: $-\tfrac12\beta_{0x}\bar\gamma_h/\operatorname{Var}x$, order $h^\alpha$ when
$\gamma(u)\approx c|u|^\alpha$.

**LVCF.** Same identity with the weight replaced by the distribution of the gap to the
last covariate time, $U_L=T-\max\{S\le T\}$:

$$\operatorname{bias}_L=-\tfrac12\beta_{0x}\,E\gamma(U_L)/\operatorname{Var}x .$$

$U_L$ has a fixed distribution (its scale is set by the visit rate, not by $n$), so
$\operatorname{bias}_L$ is a constant and $\operatorname{bias}_L/\operatorname{SE}\asymp\sqrt n$:
LVCF coverage $\to0$.

**Corrected GEE.** Adding $\tfrac12\beta_x\gamma(u)e_2$ to each pair's contribution
(identity), or replacing $\mu(\eta_S)$ by $E[\mu\{X(T)^\top\beta\}\mid x(S)]$ (logistic, under
A5 with $\rho(u)=1-\gamma(u)/2s^2$), gives $E\,U^c_h(\beta_0)=0$ **for every $h$**. Exact under
(A1)-(A3) [+A5].

So with $\gamma$ known: plain bias is $O(h^\alpha)$, LVCF bias is $O(1)$, corrected bias is $0$.

## 3. What "a proper $h$" means

Call $h$ *proper* for an estimator if, along the sequence $h=h_n$,

1. $|\operatorname{bias}(h)|<|\operatorname{bias}_L|$ (better than LVCF), and
2. $\operatorname{bias}(h)/\operatorname{SE}(h)\to0$ (Wald CI has asymptotic level $1-a$).

- Plain GEE: (1) holds iff $\bar\gamma_h<E\gamma(U_L)$, i.e. $h$ below a constant of order
  the mean LVCF gap. (2) holds iff $h=o(n^{-1/(2\alpha+1)})$. The MSE-optimal rate
  $n^{-1/(2\alpha+1)}$ is **not** proper: bias/SE $\to(2\alpha)^{-1/2}$. This is why plain CV
  coverage stalls at 88-91%.
- Corrected GEE: (1) and (2) hold for **every** $h$ in $(0,\max\text{lag}]$ with $nh\to\infty$
  (for fixed $h$, $\hat\beta^c(h)$ is $\sqrt n$-consistent). The bias constraint disappears;
  $h$ only trades variance.

The MSE-optimal $h$ is one member of the proper set, not the target. Any data-driven
$\hat h$ that stays in the proper set with probability $\to1$ gives a consistent,
asymptotically unbiased estimator.

## 4. CV for the corrected fit

Sync-loss CV (intercept of held-out pair loss on $(1,\hat\gamma(u))$) applied to corrected
fold fits estimates, up to an $h$-free constant,

$$L(h)=\delta_h^\top\Sigma_X\delta_h,\qquad \delta_h=\hat\beta^{c,(-f)}_h-\beta_0 ,$$

(Lemma 1 of the earlier note, unchanged: it holds for any fold fit). With no bias,
$E\,L(h)=\operatorname{tr}\{\Sigma_XV(h)\}/n'$, minimised at a fixed $h^*=\arg\min V(h)$ that
does not depend on $n$.

Consequences:

- $\hat h_{cvc}$ is always in the proper set (every $h$ is), so $\hat\beta^c(\hat h)$ is
  consistent and asymptotically unbiased whatever $\hat h$ is.
- The CV differences $L(h)-L(h')$ are $O_p(n^{-1})$, the same order as their noise, so
  $\hat h$ does **not** concentrate at $h^*$; it stays random as $n\to\infty$. This matches the
  simulation ($\hat h$ 10th-90th percentile 0.005-0.997 at every $n$).
- $\sqrt n(\hat\beta^c(h)-\beta_0)\Rightarrow\mathbb G(h)$, a Gaussian process in $h$. The CI
  $\hat\beta^c(\hat h)\pm z\,\widehat{\operatorname{SE}}(\hat h)$ is centred correctly but, because
  $\hat h$ depends on the same data, its level equals
  $P\{|\mathbb G(\hat h)|\le z\,\sigma(\hat h)\}$, which need not be exactly $1-a$.
  A level guaranteed for **any** selection rule: replace $z$ by the $1-a$ quantile of
  $\sup_{h\in\text{grid}}|\mathbb G(h)|/\sigma(h)$, obtainable by multiplier bootstrap of the
  per-subject scores already computed for the sandwich. Not implemented.

## 5. Gaps

1. The sandwich treats $\hat\gamma$ as known. At fixed $h$ the term
   $A^{-1}e_2\tfrac12\beta_x\,n^{-1}\sum K_h(u)\{\hat\gamma(u)-\gamma(u)\}$ is a smooth functional
   of the isotonic estimator, so it is $O_p(n^{-1/2})$, the same order as the SE. Simulation
   SE/ESE is 0.95-1.01 at $\hat h_{cvc}$, so the omitted part is small in this design, but a
   proof needs a stacked variance or a subject bootstrap that re-estimates $\hat\gamma$.
2. With bias removed by modelling rather than by localisation, validity rests entirely on
   (A1), (A3), and (A5) for logistic; Section 2 shows misspecification gives an $O(1)$ bias
   at fixed $h$, like LVCF. Robustness run (job 8449854) measures how large.

## 6. Simulation results (job 8449854, 1000 reps per cell, 0 failed fits)

Producers: `cov3.R`, `run3.R`, `agg3.R`; tables `summary3.csv`, `oracle_cmp.csv`.

**Oracle curves (correct model).** Plain bias at n=8000: -0.003 (h=0.005), -0.021 (0.038),
-0.060 (0.114), -0.343 (1.0) identity; LVCF bias -0.200 identity, -0.177 logistic, LVCF CP 0%
for n>=2000. Plain bias is below LVCF for h<=0.386, but plain CP>=93% only for h<=0.01-0.02
at n=8000. Corrected |bias|<=0.002 at every h for n>=2000; CP 94.7-96.1% (identity) and
92.1-94.8% (logistic) over the 40-point grid at n=8000. MSE ratio of h-hat_cvc to the best
fixed corrected h: 1.08-1.12 at every n.

**Misspecification, CP% at n=8000 (bias).**

| scenario | cv+corr | cvc+corr | LVCF |
|---|---|---|---|
| nonstat, identity | 97.8 (+0.001) | 97.3 (-0.001) | 0.0 (-0.086) |
| nonstat, logistic | 97.7 (-0.002) | 92.8 (-0.011) | 0.0 (-0.110) |
| merr, identity | 96.4 (+0.063) | 96.0 (-0.005) | 0.0 (-0.460) |
| merr, logistic | 87.5 (+0.173) | 55.0 (+0.051) | 0.0 (-0.361) |
| skew, identity | 96.8 (-0.001) | 95.6 (-0.007) | 0.0 (-0.355) |
| skew, logistic | 96.6 (-0.009) | 82.6 (-0.026) | 0.0 (-0.370) |

Under misspecification the corrected bias grows with h (skew logistic, fixed h, n=8000:
-0.004 at h=0.01, -0.038 at 0.1, -0.159 at 1.0), so the bias-free property of Section 3
fails and h again controls bias. cvc picks larger h than cv and loses coverage in the
logistic cells; cv+corr picks small h and keeps CP>=96.6% except under measurement error
with the logistic link.

**Bootstrap at the cv+corr bandwidth (job 8456435, correct model, 500 reps, B=200; `cov4.R`,
`agg4.R`, `summary4_*.csv`).** Re-estimating gamma-hat, mean and variance with h fixed at h-hat
gives a bootstrap SE equal to the sandwich SE (median ratio 1.00 in all four cells), so the
gamma-hat term of Gap 1 is negligible in this design. Sandwich SE/ESE is 0.98-1.09; sandwich
CP 96.2-97.2%, none significant after BH. Re-selecting h inside each resample is biased: a
duplicated subject lands in different CV folds, the resample picks a larger h (median 0.024 vs
0.020 identity, 0.047 vs 0.040 logistic at n=8000) and the SE is 4-14% too small (CP 92.0% at
logistic n=8000).

**Other links, correct model (job 8462346, 1000 reps per cell, 0 errors; `cov5.R`, `agg5.R`,
`summary5.csv`).** Correction E[mu(eta_T) | x_S] by Gauss-Hermite for every link; plain fit equals
the package fit to 1e-15 (log). n=8000, CP% (RMSE): probit cvc+corr 95.9 (0.0113), cv+corr 96.9
(0.0144), plain cv 90.0, LVCF 0.0; cloglog cvc+corr 96.5 (0.0123), cv+corr 97.8 (0.0152), plain 90.1,
LVCF 0.0; log cvc+corr 95.6, cv+corr 96.9, plain 87.6, LVCF 0.0. cvc+corr SE/ESE 0.96-1.01 for
probit and cloglog. Log link: ESE 0.10-0.19 does not shrink with n because 1.3-1.9% of reps select
an h at the bottom of the grid, where a training fold has 4-7 pairs, the fold fit is wild, and the
OLS intercept of the held-out loss on gamma-hat extrapolates to a large negative value. The same
degenerate pick occurs in 3-7% of identity reps and 0.2-0.9% of logistic reps in the earlier run.
Lemma 2 requires delta_h -> 0 on the whole grid, which fails at these h; a grid lower end set by a
diverging number of subjects per fold is needed (not yet implemented).

**Grid floor (job 8466882, five links, 1000 reps per cell, 0 errors; `cov6.R`, `agg6.R`,
`summary6.csv`).** Grid lower end = smallest h at which every training fold has at least
ceil(sqrt(n_train)) subjects with a pair inside the kernel (about 0.0013, 0.00056, 0.00027 at
n = 400, 2000, 8000). Same seeds as job 8462346. Log link ESE of cvc+corr is now 0.027, 0.012,
0.006 at n = 400, 2000, 8000 (was 0.135, 0.186, 0.102); |error| > 0.3 in at most 0.2% of reps at
n = 400 and none at n >= 2000; h-hat at the floor in at most 0.7% of reps (identity, n = 400) and
none at n = 8000. cvc+corr at n = 8000, CP% (SE/ESE, RMSE): identity 96.8 (1.03, 0.0146), logistic
97.5 (1.00, 0.0190), log 96.1 (0.99, 0.0061), probit 96.0 (0.96, 0.0113), cloglog 96.3 (0.97,
0.0123). cv+corr 96.5-97.8 with SE/ESE 1.07-1.12; plain cv 88.1-91.0; LVCF 0.0.
