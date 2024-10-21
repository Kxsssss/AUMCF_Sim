# Extensive simulation studies for the Area Under the Mean Cumulative Function (AUMCF)

## Description
This simulation study is used to evaluate the finite sample properties of AUMCF, by considering the bias, empirical standard error (ESE), asymptotic standard error (ASE), coverage probability of the 95% normal-based confidence interval (CP), and relative efficiency (RE) of the difference in AUCMF. All analyses were conducted using the [Mean Cumulative Curve (MCC)](https://github.com/zrmacc/MCC}) package.

## Controllable variables
- `n`: integer, arm size.
- `time`: double, truncation time $\tau$.
- `censor`: double, censoring rate $\lambda_C$.
- `frailtyVar`: double, variance of the gamma frailty $\sigma_0^2$.
- `BaseDeath0`: double, baseline arrival rate for the terminal event in the reference arm $\lambda_{D1}$.
- `BaseDeath1`: double, baseline arrival rate for the terminal event in the control arm $\lambda_{D2}$.
- `BetaDeath`: double, log rate ratios for the death rate $\beta_D$.
- `BaseEvent0`: double, baseline arrival rate for the recurrent events in the reference arm $\lambda_1$.
- `BaseEvent1`: double, baseline arrival rate for the recurrent events in the control arm $\lambda_2$.
- `BetaEvent`: double, log rate ratios for the event rate $\beta_E$.
- `reps`: integer, simulation replicates.
- `boot`: integer, number of bootstrap samples per replicate.
- `adjusted`: integer, indicator of adjustment. Adjusted if `adjusted = 1`; unadjusted if `adjusted = 0`.
- `tv`: numeric, true value, need to be calculated before run the simulation.
- `ei`: integer, index of experience, from 1 to 4.
- `out`: character, where to store the outputs.


## True value
The true value for each setting is the empirical average of 2,000 simulated realizations without any censoring. This is achieved by setting the censoring rate (censor) to 0 and using a very large value for the truncation time (tau), such as $10^{35}$, during the data generation step. 
 
## Simulation Design

### Experiment 1: Null Case (No Difference Between Groups)
- Objective: Evaluate bias, ESE, ASE, and CP under the null hypothesis of no difference in event rates between groups.
- Parameters:
  - $\lambda_1 = \lambda_2 =1, \lambda_D = \lambda_D = 0.2, \lambda_C = 0.2$, indicating patients are expected to live for 5 years and experience 1 recurrent event per year. 
  - No frailty or covariates
  - Vary n in {50, 100, 200, 400}.
  - Vary $\tau$ in {1, 2, 3, 4}.

### Experiment 2: Non-null Case (Difference in Event Rates)
- Objective: Introduce a difference in event rates between the groups from E1, $\lambda_1=1$ and $\lambda_2 =2$, and evaluate asymptotic unbiasedness.

### Experiment 3: Covariate Effect
- Objective: Demonstrate the effect of covariate augmentation on estimator performance.
- Parameters:
  - Fixed $\tau = 2$ under the null.
  - Two settings:
    (i) No effect of covariates, beta event rate set to be (log(1), log(1)), 
    (ii) Strong effect of covariates, beta event rate set to be (log(0.5), log(2)).
    - Note: in the code, set `BetaEvent` = 0 and `ei` = 3 to run case i, and other value of  `BetaEvent` will run case ii automatically.
  - Vary n in {50, 100, 200, 400}.
  
### Experiment 4: Frailty Effect
- Objective: Assess the impact of frailty (unobserved heterogeneity) on the AUMCF estimator's unbiasedness and coverage probability.
- Parameters:
  - Fixed $\tau = 2$ under the null.
  - Vary frailty variance in {0, 2, 4, 8}.
  - Vary n in {50, 100, 200, 400}.
  
  
  