# Extensive simulation studies for the Area Under the Mean Cumulative Function (AUMCF)

## Description
This simulation study is used to evaluate the finite sample properties of AUMCF, by considering the bias, empirical standard error (ESE), asymptotic standard error (ASE), coverage probability of the 95% normal-based confidence interval (CP), and relative efficiency (RE) of the difference in AUCMF. All analyses were conducted using the [Mean Cumulative Curve (MCC)](https://github.com/zrmacc/MCC}) package.

## Documents
- `Simulation.R`: R file includes the main code for the simulation.
- `Simulation.sh`: A shell script that can test several settings sequentially with the help of `SimConfig.txt`.

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

## Experiences

Will update later

