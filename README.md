# SSPMultiStageTwoDerivativeMethods
This is a collection of codes which we used to find Optimal SSP multiderivative multistage time stepping methods. 

opt_mdrk:         Optimization Driver file, uses Fmincon to find MDRK methods with largest SSP coeficient

Butcher2ShuOsher: Converts a scheme from Butcher representation to an optimal Shu-Osher decomposition 

mdrk_am_obj	:     Defines the objective function for Fmincon

oc_mdrk :         Creates our equality constraints for our optimization routine based off MDRK order conditions

nlc_mdrk:         Creates our equality constraints for our optimization routine based off SSP conditions

unpackMSMDRK :    Unpacks coeficients from optimization vector into Butcher Matrices

packingMSMDRK:    Transfer methods coeficients into optimization vector to feed into Fmincon

revealR:          Find largest allowable r for a given method with violating SSP conditions





