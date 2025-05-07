# Simulation studies for the Area Under the Mean Cumulative Function (AUMCF)

## Description
This simulation study is used to evaluate the finite sample properties of the difference in the AUMCF across two treatment groups, including the bias, empirical standard error (ESE), asymptotic standard error (ASE), coverage probability of the 95% normal-based confidence interval (CP), and relative efficiency (RE) of the covariate augmentation approach. All analyses were conducted using the [Mean Cumulative Curve (MCC)](https://github.com/zrmacc/MCC.git) package.

The simulation compares AUMCF-based methods with the following commonly used methods:
- Cox Proportional Hazards Model (First Event Only)
- Lin Wei Yang Ying (LWYY) Model
- Negative Binomial Model Rate Ratio
- Joint Frailty Model (without terminal event)
- Win Ratio

The implementation of these methods can be found in the files `jacc_helper.R` and `JACC_methods.R`. The main simulation script is `sim_comparison.R`.

## Simulation inputs
- `n`: integer, arm size $n$.
- `time`: double, truncation time $\tau$.
- `censor`: double, censoring rate $\lambda_C$.
- `frailtyVar`: double, variance of the gamma frailty $\sigma_0^2$.
- `BaseDeath0`: double, baseline arrival rate for the terminal event in the reference arm $\lambda_{D1}$.
- `BaseDeath1`: double, baseline arrival rate for the terminal event in the control arm $\lambda_{D2}$.
- `BetaDeath`: double, log rate ratio for the death rate $\beta_D$.
- `BaseEvent0`: double, baseline arrival rate for the recurrent events in the reference arm $\lambda_1$.
- `BaseEvent1`: double, baseline arrival rate for the recurrent events in the control arm $\lambda_2$.
- `BetaEvent`: double, log rate ratio for the event rate $\beta_E$.
- `reps`: integer, simulation replicates.
- `adjusted`: integer, indicator of covariate adjustment. Adjusted if `adjusted = 1`; unadjusted if `adjusted = 0`.
- `tv`: double, true value, needs to be calculated before running the simulation.
- `experiment`: integer, index of the simulation experiment, from 1 to 3. For more detailed information about the experiment, please see below.
- `out`: character, where to store the outputs.


## True value
When the true value is not 0, it is the empirical average of 2,000 (`reps=2000`) simulated realizations without any censoring for each setting. This is achieved by setting the censoring rate (`censor`) to 0 and using very large values for the truncation time (`time`) and arm size (`n`), such as `time` = $10^{35}$ and `n` = $10^4$, during the data generation step (`Gen_data()`).

## Simulation Design
Each scenario is based on a fixed number of replicated datasets.

### Validity (E1)
Evaluate bias, coverage probability of 95% CI, and type I error under the null hypothesis (no treatment effect).
- Event times: Poisson(λ_E = 1)
- Censoring: Exponential(λ_C = 0.2)
- Death: Exponential(λ_D = 0.2)
- Sample sizes: n ∈ {50, 100, 200, 400}
- Truncation times: τ ∈ {1, 2, 3, 4}
- Replicates: 10,000

### Power (E2)
Assess estimation and inference when a treatment effect is present.
- Event rates: $λ_{E1}$ = 1 (treatment), $λ_{E2}$ = 2 (reference)
- Other settings same as Validity
- Replicates: 1,000

### Covariate Adjustment (E3)
Evaluate the effect of covariate augmentation under the null (no treatment effect).
- Covariate: X ∼ N(0, 1)
- Two settings:
  1. Uninformative: no effect on events or death
  2. Informative:
     - Death rate: scaled by $exp(X β_D), β_D = log(0.5)$
     - Event rate: scaled by $exp(X β_E), β_E = log(2.0)$
  - Note: in the code, set BetaEvent = 0 and experiment = 3 to run case (1), and other value of BetaEvent will run case (2) automatically
- Truncation time: τ = 2
- Sample sizes: n ∈ {50, 100, 200, 400}

## Output Files
Each run generates two `.rds` files:

### 1. Simulation-level estimates (contains per-replicate results)
File name pattern:
```
simN<n>_T<time>_l0<BaseEvent0>_l1<BaseEvent1>_adj<adjusted>.rds
```
Includes the following columns:
- `value`: point estimate of AUMCF difference
- `se`: standard error
- `lower`, `upper`: confidence interval bounds
- `p_value`: p-value from hypothesis test
- `type`: method name
- `true_value`: true value for the setting

### 2. Summary statistics
File name pattern:
```
N<n>_T<time>_l0<BaseEvent0>_l1<BaseEvent1>_adj<adjusted>.rds
```
Includes the following columns:
- `type`: method name
- `bias`: average bias
- `lower`, `upper`: average CI bounds
- `cov_p`: empirical coverage
- `ase`: average ASE
- `ese`: empirical standard deviation of estimates
- `p_value`: average p-value
- `true_value`: setting-specific true value
- `n`: sample size
- `time`: truncation time
- `rep`: number of replicates

## Running the Code

### Local Execution
Run `sim_comparison.R` after setting parameters in the `params` list. Output will be saved in the specified format.

### Cluster Execution
Use the provided `Sim.sh` file as a template (configured for Niagara cluster). Adjust it to meet your cluster environment's requirements.
