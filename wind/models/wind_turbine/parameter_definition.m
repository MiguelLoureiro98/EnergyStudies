%
%
%

%% Wind generation parameters.

Vm = 6;

Lm = 150;
sigma = 0.15;
Tv = Lm / Vm; 

a1 = 0.4;
a2 = 0.25;
Kv = sqrt(2 * Tv * (1 - a2^2) / (a1^2/a2 - a2 + 1 - a1^2));

Noise_power = 1;
Noise_Ts = 0.001;

%% Aerodynamic properties.

