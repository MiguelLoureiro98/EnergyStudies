%
%
%

function power_coef = trigonometric_approximation(TSR, beta, c1, c2, c3, c4, c5, c6)

    sin_arg = (pi .* (TSR - c3)) ./ (c4 - c5 .* beta);

    power_coef = (c1 - c2 .* beta) .* sin(sin_arg) - c6 .* beta .* (TSR - c3);

end