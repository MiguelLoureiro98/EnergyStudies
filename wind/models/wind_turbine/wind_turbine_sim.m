%
%
%

clear
clc

%% Run necessary definition scripts.

fit_var_name = "Cp";
run data_preparation.m
run surface_fitting.m

Cp_coef = coeffvalues(surf_fit);

fit_var_name = "Ct";
run data_preparation.m
run surface_fitting.m
Ct_coef = coeffvalues(surf_fit);

run parameter_definition.m

%% Define simulation parameters.

sim_time = 10;