% Performs a surface fit on the heat pump's COP data with Tin and Tout
%   (the evaporator and condenser inlet temperatures) as the dependent
%   variables.

function [surf, good] = fit_COP_data_Hitachi(COP, Tin, Tout)

    [surf, good] = fit([Tin, Tout], COP, 'poly22');

end