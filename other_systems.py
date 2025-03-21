import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def pitch_actuator(t: float,
                   x: list,
                   u: float,
                   eps: float=0.6,
                   wn: float=11.11) -> list:

    dbetadot = wn**2 * u - 2 * eps * wn * x[0] - wn**2 * x[1];
    dbeta = x[0];

    return [dbetadot, dbeta];

def pitch_controller(x):

    return 1.0;

def yaw_actuator(t: float,
                 x: list,
                 T_yaw: float,
                 J_yaw: float=2607890,
                 D_yaw: float=19160000,
                 K_yaw: float=9028320000) -> list:
    
    dgamadot = (K_yaw * T_yaw - D_yaw * x[0] - K_yaw * x[1]) / J_yaw;
    dgama = x[0];

    return [dgamadot, dgama];

def heat_pump(t: float,
              x: list,
              N_people: int,
              mdot_HP: float,
              mdot_rad: float,
              Q: float,
              T_outside: float) -> list:
    
    dT_HP = (mdot_HP * 4186 * (x[1] - x[0]) + 4.9 * Q) / 0.3153;
    dT_BT = (mdot_HP * 4186 * (x[0] - x[1]) + mdot_rad * 4186 * (x[2] - x[1])) / (0.2 * 4186.8 * 998);
    dT_rad = (mdot_rad * 4186 * (x[1] - x[2]) + 500 * (x[2] - x[1]) / np.log((x[3] - x[2]) / (x[3] - x[1]))) / (46.7 * 447);
    dT_zone = (4 * 2 * 12 * (x[5] - x[3]) + 1 * 9 * (x[4] - x[3]) - 500 * (x[2] - x[1]) / np.log((x[3] - x[2]) / (x[3] - x[1]))) / 47100;
    dT_roof = 1 * 9 * (x[0] - 2 * x[1] + T_outside) / 80000;
    dT_walls = 2 * 12 * (x[0] - 2 * x[2] + T_outside) / 65000;

    return [dT_HP, dT_BT, dT_rad, dT_zone, dT_roof, dT_walls];

if __name__ == "__main__":

    # Pitch Actuator     
    # Initial conditions
    x0 = np.zeros(2).tolist();

    # Time span
    t_span = [0, 10];
    Ts = 0.01;
    n = int((t_span[1] - t_span[0]) / Ts + 1);

    if(int((t_span[1] - t_span[0]) % Ts) != 0):
    
        n += 1;
    
    true_time = np.array([np.fmin(t_span[1], t_span[0] + Ts * i) for i in range(n)]);

    results = np.zeros((3, 1));

    for t in range(n-1):

        # Control input
        u = pitch_controller(x0);

        # Solve ODE
        sol = solve_ivp(pitch_actuator, [true_time[t], true_time[t+1]], x0, args=(u, 0.6, 11.11));

        # Store results
        new_data = np.vstack((sol.t, sol.y));
        results = np.hstack((results, new_data));

        # Update initial conditions
        x0 = sol.y[:, -1];

    # Plot results
    fig = plt.figure(figsize=(5, 5));

    ax = fig.subplots(2, 1, sharex=True);

    

    ax[0].plot(results[0, :], results[1, :]);
    ax[0].set_title("Beta_dot");
    ax[0].grid();

    ax[1].plot(results[0, :], results[2, :]);
    ax[1].set_title("Beta");
    ax[1].grid();
    ax[1].set_xlabel("t (s)");

    plt.show();

    # Yaw Actuator     
    # Initial conditions
    x0 = np.zeros(2).tolist();

    # Time span
    t_span = [0, 10];
    Ts = 0.01;
    n = int((t_span[1] - t_span[0]) / Ts + 1);

    if(int((t_span[1] - t_span[0]) % Ts) != 0):
    
        n += 1;
    
    true_time = np.array([np.fmin(t_span[1], t_span[0] + Ts * i) for i in range(n)]);

    results = np.zeros((3, 1));

    for t in range(n-1):

        # Control input
        u = pitch_controller(x0);

        # Solve ODE
        sol = solve_ivp(yaw_actuator, [true_time[t], true_time[t+1]], x0, args=(u,));

        # Store results
        new_data = np.vstack((sol.t, sol.y));
        results = np.hstack((results, new_data));

        # Update initial conditions
        x0 = sol.y[:, -1];

    # Plot results
    fig = plt.figure(figsize=(5, 5));

    ax = fig.subplots(2, 1, sharex=True);

    

    ax[0].plot(results[0, :], results[1, :]);
    ax[0].set_title("Gama_dot");
    ax[0].grid();

    ax[1].plot(results[0, :], results[2, :]);
    ax[1].set_title("Gama");
    ax[1].grid();
    ax[1].set_xlabel("t (s)");

    plt.show();

    # Heat pump     
    # Initial conditions
    x0 = np.zeros(6).tolist();

    # Time span
    t_span = [0, 10];
    Ts = 0.01;
    n = int((t_span[1] - t_span[0]) / Ts + 1);

    if(int((t_span[1] - t_span[0]) % Ts) != 0):
    
        n += 1;
    
    true_time = np.array([np.fmin(t_span[1], t_span[0] + Ts * i) for i in range(n)]);

    results = np.zeros((7, 1));

    for t in range(n-1):

        # Solve ODE
        sol = solve_ivp(yaw_actuator, [true_time[t], true_time[t+1]], x0, args=(u,));

        # Store results
        new_data = np.vstack((sol.t, sol.y));
        results = np.hstack((results, new_data));

        # Update initial conditions
        x0 = sol.y[:, -1];

    # Plot results
    fig = plt.figure(figsize=(5, 5));

    ax = fig.subplots(6, 1, sharex=True);

    ax[0].plot(results[0, :], results[1, :]);
    ax[0].set_title("Gama_dot");
    ax[0].grid();

    ax[1].plot(results[0, :], results[2, :]);
    ax[1].set_title("Gama");
    ax[1].grid();
    
    ax[1].set_xlabel("t (s)");

    plt.show();