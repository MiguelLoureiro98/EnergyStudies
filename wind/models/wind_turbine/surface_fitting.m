%
%
%

%% Select variable to fit.

%if(strcmp(fit_var_name, "Cp") == 1)

%    fit_var = Cp;

%elseif(strcmp(fit_var_name, "Ct") == 1)

%    fit_var = Ct;

%else

%    fprintf(2, "A valid variable should be selected.\n");

%end

%fit_var_name = "Cp";
fit_var = Cp;

%% Select fit type.

fit_type = "poly";
deg_TSR = 5;
deg_beta = 5;

%% Fit surface to the data.

R = 125.88 / 2;

if(strcmp(fit_type, "exp") == 1)

    ft = fittype('exponential_approximation(TSR, beta, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13)', ...
        "independent", {'TSR', 'beta'});
    options = fitoptions(ft);
    options.StartPoint = [1, 0.01, 0, 1, 0.01, 0, 0, 0, 0, 0.01, 0.01, 0, 0];
    options.Lower = [-2, -2, 0, 1, -2, 0, 0, 0, 0, -2, -2, 0, 0];
    options.Upper = [2, 2, 0, 1, 2, 0, 0, 0, 0, 2, 2, 0, 0];
    options.Algorithm = "Trust-Region";
    options.Display = "Iter";
    surf_fit = fit([TSR, beta], fit_var, ft, options);

elseif(strcmp(fit_type, "poly") == 1)

    surf_fit = fit([TSR, beta], fit_var, "poly" + num2str(deg_TSR) + num2str(deg_beta));

elseif(strcmp(fit_type, "trig") == 1)

    ft = fittype("trigonometric_approximation(TSR, beta, c1, c2, c3, c4, c5, c6)", ...
         "independent", {'TSR', 'beta'});
    options = fitoptions(ft);
    options.StartPoint = [0.01, 0.01, 0.01, 1, 0.01, 0.01];
    options.Lower = [-2, -2, -2, 1, -2, -2];
    options.Upper = [2, 2, 2, 1, 2, 2];
    options.Algorithm = "Trust-Region";
    options.Display = "Iter";
    surf_fit = fit([TSR, beta], fit_var, ft, options);

else

    fprintf(2, "There are only three options for the surface fit: " + ...
               "'exp', 'poly' and 'trig'.");

end

%% Plotting.

figure
plot(surf_fit, [TSR, beta], fit_var), xlabel("TSR"), ylabel("beta"), ...
    zlabel(fit_var_name), grid on, ...
    title("Wind turbine aerodynamic coefficients from data");