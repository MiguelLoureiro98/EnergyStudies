%
%
%

function power_coef = exponential_approximation(TSR, beta, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13)

    lambda_i_inv = 1 ./ (c10 .* TSR + c11 .* beta) - c12 ./ (1 + (c13 .* beta).^3);
    lambda_i = 1 ./ lambda_i_inv;

    exponent = -c7 .* (lambda_i_inv - c8);

    power_coef = c1 .* (c2 ./ lambda_i - c3 .* beta - c4 .* beta.^c5 - c6) .* exp(exponent) + c9 .* TSR;

end