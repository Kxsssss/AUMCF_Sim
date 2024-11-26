# Simulation studies for the Area Under the Mean Cumulative Function (AUMCF)

## Description
This simulation study is used to evaluate the finite sample properties of the difference in the AUMCF across two treatment groups, including the bias, empirical standard error (ESE), asymptotic standard error (ASE), coverage probability of the 95% normal-based confidence interval (CP), and relative efficiency (RE) of the covariate augmentation approach. All analyses were conducted using the [Mean Cumulative Curve (MCC)](https://github.com/zrmacc/MCC}) package.

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
- `tv`: numeric, true value, needs to be calculated before running the simulation.
- `ei`: integer, index of the simulation experiment, from 1 to 4.
- `out`: character, where to store the outputs.


## True value
The true value for each setting is the empirical average of 2,000 simulated realizations without any censoring. This is achieved by setting the censoring rate (`censor`) to 0 and using a very large value for the truncation time (`time`), such as $10^{35}$, during the data generation step (`Gen_data()`). 
 
## Simulation Design
Each simulation study is based on 1000 replicated datasets (`reps=1000`).

### Experiment 1: Null Case (No Difference Between Groups)
- Objective: Evaluate bias, ESE (empirical standard error), ASE (asymptotic standard error), CP (coverage probability), p-value, and MSE (mean square error) under the null hypothesis of no difference in event rates between groups.
- Parameters:
  - $\lambda_1 = \lambda_2 =1, \lambda_D = \lambda_D = 0.2, \lambda_C = 0.2$, indicating patients are expected to live for 5 years and experience 1 recurrent event per year
  - No frailty or covariates
  - Vary n in {50, 100, 200, 400}
  - Vary $\tau$ in {1, 2, 3, 4}

### Experiment 2: Non-null Case (Difference in Event Rates)
- Objective: Introduce a difference in event rates between the groups from E1, with $\lambda_1=1$ and $\lambda_2=2$.  The rest of the settings remain the same.

### Experiment 3: Covariate Effect
- Objective: Demonstrate the effect of covariate augmentation on estimator performance.
- Parameters:
  - Fixed $\tau = 2$ under the null
  - Two settings:
    - (i) No effect of covariates: beta event rate set to be (log(1), log(1)) 
    - (ii) Strong effect of covariates: beta event rate set to be (log(0.5), log(2))
    - Note: in the code, set `BetaEvent` = 0 and `ei` = 3 to run case (i), and other value of  `BetaEvent` will run case (ii) automatically
  - Vary n in {50, 100, 200, 400}
  
### Experiment 4: Frailty Effect (Additional simulation not included in the text)
- Objective: Assess the impact of frailty (unobserved heterogeneity) on estimator performance. 
- Parameters:
  - Fixed $\tau = 2$ under the null
  - Vary frailty variance in {0, 2, 4, 8}
  - Vary n in {50, 100, 200, 400}

## How to run the code

### On your computer
If you select all and run the `Simulation.R` file, it will run using the default values in the `params`.  The current setup is under the null. 

The output includes:

- two RDS files containing:
  - A summary table, including...
  - All the simulation data and the running time information.  The name of the file starts with "sim". 

If you change the default values of `params` (in the **Command line arguments** section) and "select all" to run each time, the simulation setting will change.

### On a cluster
Since different clusters may require different `.sh` formats, `Sim.sh` is an example that works on the Niagara cluster. Please adapt the format according to the requirements of your specific cluster.

  
  
