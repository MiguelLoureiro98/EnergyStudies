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

def COP_model(T_outside: float, T_BT: float) -> float:

    COP = np.fmin(4.9, 5.71 + 0.1615 * T_outside - 0.1236 * T_BT + 0.001217 * T_outside**2 - 0.002318 * T_BT**2 + 0.0008342 * T_outside * T_BT);

    return COP;

def log_mean_T(Tzone: float, Tin: float, Tout: float) -> float:

    if(np.sign(Tin) == np.sign(Tout)):

        T = 500 * (Tin - Tout) / (np.log((Tzone - Tout + 1e-7) / (Tzone - Tin + 1e-7)));
    
    else:

        T = 0.0;

    return T;

def heat_pump(t: float,
              x: list,
              N_people: int,
              mdot_HP: float,
              mdot_rad: float,
              Q: float,
              T_outside: float) -> list:
    
    dT_HP = (mdot_HP * 4186 * (x[1] - x[0]) + COP_model(T_outside, x[1]) * Q) / 0.3153;
    dT_BT = (mdot_HP * 4186 * (x[0] - x[1]) + mdot_rad * 4186 * (x[2] - x[1])) / (0.2 * 4186.8 * 998);
    dT_rad = (mdot_rad * 4186 * (x[1] - x[2]) + log_mean_T(x[3], x[1], x[2])) / (46.7 * 447);
    dT_zone = (4 * 2 * 12 * (x[5] - x[3]) + 1 * 9 * (x[4] - x[3]) - log_mean_T(x[3], x[1], x[2])) / 47100;
    dT_roof = 1 * 9 * (x[3] - 2 * x[4] + T_outside) / 80000;
    dT_walls = 2 * 12 * (x[3] - 2 * x[5] + T_outside) / 65000;

    return [dT_HP, dT_BT, dT_rad, dT_zone, dT_roof, dT_walls];

def XY_table(t: float,
             x: list,
             Tm_x: float,
             Tm_y: float) -> list:
    
    ddpsi_x = (Tm_x + 2.6682 * (x[3] - x[1]) + 5e-4 * (x[2] - x[0]) - 1.39e-4 * x[0]) / 2.3e-5;
    dpsi_x = x[0];
    ddphi_x = (2.6682 * (x[1] - x[3]) + 5e-4 * (x[0] - x[2]) - 0.001 * x[2]) / 2.1e-5;
    dphi_x = x[2];
    ddpsi_y = (Tm_y + 2.3017 * (x[7] - x[5]) + 0.0011 * (x[6] - x[4]) - 1.17e-4 * x[4]) / 2.3e-5;
    dpsi_y = x[4];
    ddphi_y = (2.3017 * (x[5] - x[7]) + 0.0011 * (x[4] - x[6]) - 7.25e-4 * x[6]) / 2.33e-5;
    dphi_y = x[6];
    X = 0.7958 * x[3];
    Y = 0.7958 * x[7];

    return [ddpsi_x, dpsi_x, ddphi_x, dphi_x, ddpsi_y, dpsi_y, ddphi_y, dphi_y, X, Y];

def tower_crane():

    pass

def CSTR():

    pass

def bioreactor():
    
    pass

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
    x0 = [50.0, 40.0, 30.0, 25.0, 20.0, 20.0];

    # Time span
    t_span = [0, 20000];
    Ts = 10.0;
    n = int((t_span[1] - t_span[0]) / Ts + 1);

    if(int((t_span[1] - t_span[0]) % Ts) != 0):
    
        n += 1;
    
    true_time = np.array([np.fmin(t_span[1], t_span[0] + Ts * i) for i in range(n)]);

    results = np.zeros((7, 1));

    #for t in range(n-1):

        # Solve ODE
    #    sol = solve_ivp(heat_pump, [true_time[t], true_time[t+1]], x0, args=(0, 1, 1, 1000, 10.0));

        # Store results
    #    new_data = np.vstack((sol.t, sol.y));
    #    results = np.hstack((results, new_data));

        # Update initial conditions
    #    x0 = sol.y[:, -1];

    # Plot results
    #fig = plt.figure(figsize=(10, 5));

    #ax = fig.subplots(6, 1, sharex=True);

    #ax[0].plot(results[0, :], results[1, :]);
    #ax[0].set_title("T Heat Pump");
    #ax[0].grid();

    #ax[1].plot(results[0, :], results[2, :]);
    #ax[1].set_title("T Buffer Tank");
    #ax[1].grid();

    #ax[2].plot(results[0, :], results[3, :]);
    #ax[2].set_title("T Radiator");
    #ax[2].grid();

    #ax[3].plot(results[0, :], results[4, :]);
    #ax[3].set_title("T Zone");
    #ax[3].grid();

    #ax[4].plot(results[0, :], results[5, :]);
    #ax[4].set_title("T Roof");
    #ax[4].grid();

    #ax[5].plot(results[0, :], results[6, :]);
    #ax[5].set_title("T Walls");
    #ax[5].grid();
    #ax[5].set_xlabel("t (s)");

    #plt.show();

    # Pitch Actuator     
    # Initial conditions
    x0 = np.zeros(10).tolist();

    # Time span
    t_span = [0, 10];
    Ts = 0.01;
    n = int((t_span[1] - t_span[0]) / Ts + 1);

    if(int((t_span[1] - t_span[0]) % Ts) != 0):
    
        n += 1;
    
    true_time = np.array([np.fmin(t_span[1], t_span[0] + Ts * i) for i in range(n)]);

    results = np.zeros((11, 1));

    for t in range(n-1):

        # Solve ODE
        sol = solve_ivp(XY_table, [true_time[t], true_time[t+1]], x0, args=(0.0, 0.0));

        # Store results
        new_data = np.vstack((sol.t, sol.y));
        results = np.hstack((results, new_data));

        # Update initial conditions
        x0 = sol.y[:, -1];

    # Plot results
    fig = plt.figure(figsize=(5, 5));

    ax = fig.subplots(2, 1, sharex=True);

    ax[0].plot(results[0, :], results[9, :]);
    ax[0].set_title("X");
    ax[0].grid();

    ax[1].plot(results[0, :], results[10, :]);
    ax[1].set_title("Y");
    ax[1].grid();
    ax[1].set_xlabel("t (s)");

    plt.show();