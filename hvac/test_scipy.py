import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def eps(NTU: float, Cr: float) -> float:

    if(NTU == 0.0 and Cr == 1.0):

        epsilon = 0.0;
    
    else:

        epsilon = (1 - np.exp(-NTU * (1 - Cr))) / (1 - Cr * np.exp(-NTU * (1 - Cr)));

    return epsilon;

def Cr(Cmin: float, Cmax: float) -> float:

    if(Cmax == 0.0):

        Cr = 1.0;
    
    else:

        Cr = Cmin / Cmax;

    return Cr;

def Cminmax(mdot_1: float, cp_1: float, mdot_2: float, cp_2: float) -> float:

    Cmin = np.fmin(mdot_1 * cp_1, mdot_2 * cp_2);
    Cmax = np.fmax(mdot_1 * cp_1, mdot_2 * cp_2);

    return Cmin, Cmax;

def efficiency(mdot_1: float, 
               cp_1: float, 
               mdot_2: float, 
               cp_2: float,
               component: str) -> float:

    Cmin, Cmax = Cminmax(mdot_1, cp_1, mdot_2, cp_2);
    Cratio = Cr(Cmin, Cmax);

    correlations = {
        "cooling_coil": (4.6879 * mdot_1 ** (0.472)) / (1 + 15.6315 * (mdot_1 / (mdot_2 + 1e-12)) ** (0.472)), # 1 - air; 2 - water
        "dehumidifier": (1.993 * mdot_1 ** (1.1324)) / (1 + 1.7540 * (mdot_1 / (mdot_2 + 1e-12)) ** (1.1324)), # 1 - air; # 2 - solution
        "cooler": (12.8299 * mdot_1 ** (0.8505)) / (1 + 4.888 * (mdot_1 / (mdot_2 + 1e-12)) ** (0.8505)) # 1 - solution; 2 - water
    };

    if(Cmin == 0.0):

        NTU = 0.0;
    
    else:
        
        NTU = correlations[component] / Cmin;

    eff = eps(NTU, Cratio);

    return eff;

def fan_power(mdot: float | np.ndarray) -> float | np.ndarray:

    power = 0.0;

    if(mdot != 0.0):

        power = -40803.019 * mdot**2 + 6737.594 * mdot + 157.166;

    return power;

def cooling_coil_power(mdot: float | np.ndarray, 
                       Tin: float | np.ndarray, 
                       Tout: float | np.ndarray) -> float | np.ndarray:

    power = 0.0;

    return power;

def pump_power(mdot: float | np.ndarray) -> float | np.ndarray:

    power = mdot**3 - 30927273.097 * mdot**2 + 137809.092 * mdot + 157.5;

    return power;

def dehumidifier_power(mdot: float | np.ndarray, 
                       Tsin: float | np.ndarray, 
                       Tsout: float | np.ndarray) -> float | np.ndarray:

    COP_chiller = 3;

    power = mdot * (Tsin - Tsout) / COP_chiller;

    return power;

def HVAC_system(t: np.ndarray, 
                x: list, 
                mdot_sa: float,
                d: float,
                mdot_cooling_coil: float, 
                mdot_dehum: float, 
                mdot_cooler: float,
                Q_electric: float,
                w_equ: float, 
                T_out: float, 
                w_out: float, 
                CO2_out: float,
                N: int) -> float:

    cp_air = 1005;
    cp_vapour = 1820;
    cp_water = 4186.8;
    cp_solution = 4027;
    rho_air = 1.225;

    w_new = d * w_out + (1 - d) * x[3];

    cp_air_new = cp_air + w_new * cp_vapour;
    cp_air_sa = cp_air + x[7] * cp_vapour;
    cp_air_zone = cp_air + x[3] * cp_vapour;

    cp_air_out = cp_air + w_out * cp_vapour;
    T_new = (d * cp_air_out * T_out + (1 - d) * cp_air_zone * x[0]) / cp_air_new;

    eff_cc = efficiency(mdot_sa, cp_air_new, mdot_cooling_coil, cp_water, "cooling_coil");
    Cmin_cc, _ = Cminmax(mdot_sa, cp_air_new, mdot_cooling_coil, cp_water);
    eff_deh = efficiency(mdot_sa, cp_air_new, mdot_dehum, cp_solution, "dehumidifier");
    Cmin_deh, _ = Cminmax(mdot_sa, cp_air_new, mdot_dehum, cp_solution);
    eff_cooler = efficiency(mdot_dehum, cp_solution, mdot_cooler, cp_water, "cooler");
    Cmin_cooler, _ = Cminmax(mdot_dehum, cp_solution, mdot_cooler, cp_water);

    dT_zone = (mdot_sa * cp_air_zone * (x[10] - x[0]) + 4 * 2 * 12 * (x[2] - x[0]) + 1 * 9 * (x[1] - x[0]) + N * 185) / 47100;
    dT_roof = 1 * 9 * (x[0] - 2 * x[1] + T_out) / 80000;
    dT_walls = 2 * 12 * (x[0] - 2 * x[2] + T_out) / 65000;
    dw_zone = mdot_sa * (x[7] - x[3]) / (rho_air * 1.44) + N * (0.025 / 3600);
    dCO2_zone = mdot_sa * (CO2_out - x[4]) / (rho_air * 1.44) + N * 0.554;
    dT_cool = (eff_cc * Cmin_cc * (10 - T_new) + mdot_sa * cp_air_new * (T_new - x[5])) / (7.9941 * cp_air_new);
    dT_dehum = (eff_deh * Cmin_deh * (x[9] - x[5]) + mdot_sa * cp_air_new * (x[5] - x[6])) / (10.798 * cp_air_new);
    dw_dehum = (0.8 * (w_equ - w_out) + mdot_sa * (w_out - x[7])) / (10.798);
    dT_s_dehum = (eff_deh * Cmin_deh * (x[5] - x[9]) + mdot_dehum * cp_solution * (x[9] - x[8]) + mdot_sa * 2501300 * (w_out - x[7])) / (cp_solution * 15.6571);
    dT_s_cooler = (eff_cooler * Cmin_cooler * (1 - x[9]) + mdot_dehum * cp_solution * (x[8] - x[9])) / (cp_solution * 15.6571);
    dT_heat = (mdot_sa * cp_air_sa * (x[6] - x[10]) + 0.8 * Q_electric) / 4500;

    return [dT_zone, dT_roof, dT_walls, dw_zone, dCO2_zone, dT_cool, dT_dehum, dw_dehum, dT_s_dehum, dT_s_cooler, dT_heat];

def HVAC_controller(x):

    return (20, 1, 0, 0, 0, 0, 0.3);

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

    #print(sol.t);

    # Plot results
    #plt.plot(results[0, :], results[1, :]);
    #plt.show();

    # HVAC system

    # Initial conditions
    x2 = [25, 20, 20, 0.2, 0, 15, 15, 0.2, 15, 15, 25];

    # Time span
    t_span = [0, 100];
    Ts = 0.1;
    n = int((t_span[1] - t_span[0]) / Ts + 1);

    if(int((t_span[1] - t_span[0]) % Ts) != 0):
    
        n += 1;
    
    true_time = np.array([np.fmin(t_span[1], t_span[0] + Ts * i) for i in range(n)]);

    results = np.array([[0.0], [25.0], [20.0], [20.0], [0.2], [0.0], [15.0], [15.0], [0.2], [15.0], [15.0], [25.0]]);

    for t in range(n-1):

        # Control input
        mdot_sa, d, mdot_cooling_coil, mdot_dehum, mdot_cooler, Q, w_equ = HVAC_controller(x2);

        # Solve ODE
        sol = solve_ivp(HVAC_system, [true_time[t], true_time[t+1]], x2, args=(mdot_sa, d, mdot_cooling_coil, mdot_dehum, mdot_cooler, Q, w_equ, 15, 0.3, 0, 0));

        # Store results
        new_data = np.vstack((sol.t, sol.y));
        results = np.hstack((results, new_data));

        # Update initial conditions
        x2 = sol.y[:, -1];

    # Plot results
    plt.plot(results[0, :], 15 * np.ones(results.shape[1]), label="Outside Temperature");
    plt.plot(results[0, :], results[1, :], label="Zone Temperature");
    plt.legend();
    plt.show();