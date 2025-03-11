%% Create variables.

run initialise_vars.m;

%% Surface fitting.

COP = COP / max(COP) * 4.9;
[fit, gof] = fit_COP_data_Hitachi(COP, Tin, Tout);

%% Goodness of fit.

fprintf(1, "Goodness of fit: %f\n", gof.adjrsquare);

%% Plots.

figure;
plot(fit, [Tin, Tout], COP);

figure;
plot(fit, [Tin, Tout], COP, "Style", "Residuals");