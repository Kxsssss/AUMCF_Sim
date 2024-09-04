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
- `adj`: integer, indicator of adjustment. Adjusted if `adj = 1`; unadjusted if `adj = 0`.
- `out`: character, where to store the outputs.

## True value
The true value of the parameter can be calculated when specifying $\lambda_C=0$ and `reps = 10,000,000`, and replacing the `Loop()` function with `Loop_true()` in the $\textit{Simulation}$ section of `Simulation.R` as shown below:

```r
Loop_true <- function(i) {
  data <- Gen_data(params)
  boot <- try(
    MCC::CompareAUCs(
      tau = params$time,
      # no censoring
      data %>% dplyr::filter(status != 0),
      cens_after_last = FALSE
    )
  )
  if (class(boot) != "try-error") {
    return(boot@CIs$observed[1])
  }
}
```
Once the 10 million estimates are obtained, their average value will be used to represent the true value. Additionally, please replace the value of `true_value` (in the $\textit{Summarize} section of `Simulation.R`) under different parameter settings.

