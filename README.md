# Extensive simulation studies for the Area Under the Mean Cumulative Function (AUMCF)

## Description
This simulation study is used to evaluate the finite sample properties of AUMCF, by considering the bias, empirical standard error (ESE), asymptotic standard error (ASE), coverage probability of the 95% normal-based confidence interval (CP), and relative efficiency (RE) of the difference in AUCMF. All analyses were conducted using the [Mean Cumulative Curve (MCC)](https://github.com/zrmacc/MCC}) package.

## Documents
- `Simulation.R`: R file includes the main code for the simulation with unadjusted estimators.
- `Simulation_adj.R`: R file includes the main code for the simulation with adjusted estimators.
- `Simulation.sh`: A shell script that can test several settings sequentially with the help of `SimConfig.txt`.