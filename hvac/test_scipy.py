import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def HVAC_system(t, x, mdot_sa, d, mdot_cool, mdot_dehum, Q_electric, T_out, w_out, CO2_out):

    dT_zone = (mdot_sa * 1005 * (x[10] - x[0]) + 4 * 2 * 12 * (x[2] - x[0]) + 1 * 9 * (x[1] - x[0])) / 47100;
    dT_roof = 1 * 9 * (x[0] - 2 * x[1] + T_out) / 80000;
    dT_walls = 2 * 12 * (x[0] - 2 * x[2] + T_out) / 65000;
    dw_zone = mdot_sa * (x[8] - x[3]) / (1.225 * 1.44);
    dCO2_zone = mdot_sa * (CO2_out - x[9]) / (1.225 * 1.44);
    dT_cool = 0;
    dT_dehum = 0;
    dw_dehum = 0;
    dT_s_dehum = 0;
    dT_s_cooler = 0;
    dT_heat = (mdot_sa * 1005 * (x[6] - x[10]) + 0.8 * Q_electric) / 4500;

    return [dT_zone, dT_roof, dT_walls, dw_zone, dCO2_zone, dT_cool, dT_dehum, dw_dehum, dT_s_dehum, dT_s_cooler, dT_heat];

def first_order(t, x, u):

    # x[0] = y
    dydt = u - x[0];

    return [dydt];

def controller(x):

    u = 1;

    return u;

if __name__ == "__main__":

    # Initial conditions
    x0 = [0];

    # Time span
    t_span = [0, 10];
    Ts = 0.1;
    n = int((t_span[1] - t_span[0]) / Ts + 1);

    if(int((t_span[1] - t_span[0]) % Ts) != 0):
    
        n += 1;
    
    true_time = np.array([np.fmin(t_span[1], t_span[0] + Ts * i) for i in range(n)]);

    results = np.array([[0.0], [0.0]]);

    for t in range(n-1):

        # Control input
        u = controller(x0);

        # Solve ODE
        sol = solve_ivp(first_order, [true_time[t], true_time[t+1]], x0, args=(u,));

        # Store results
        new_data = np.vstack((sol.t, sol.y[0]));
        results = np.hstack((results, new_data));

        # Update initial conditions
        x0 = sol.y[:, -1];

print(sol.t);

# Plot results
plt.plot(results[0, :], results[1, :]);
plt.show();