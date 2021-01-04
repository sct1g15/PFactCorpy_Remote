from math import exp

#SM8 Model
###Fitting Params###
Beta_SASA = 11.6408
Gamma_SASA = 1.2054
Delta_SASA = 10.0026
Beta_DIST = 9.5150
Gamma_DIST = 1.1609
Delta_DIST = 6.8153

#Model Run
def run_SM8(SASAi, DISTi):
    lnPF = max(Beta_SASA/(1+exp(Gamma_SASA*(SASAi - Delta_SASA))), Beta_DIST/(1+exp(Gamma_DIST*(DISTi - Delta_DIST))))
    if Beta_SASA/(1+exp(Gamma_SASA*(SASAi - Delta_SASA))) > Beta_DIST/(1+exp(Gamma_DIST*(DISTi - Delta_DIST))):
        print('SASA contrib')
    else:
        print('DIST contrib')
    return lnPF