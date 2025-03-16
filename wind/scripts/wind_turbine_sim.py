import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def Ct_TSR(TSR: float| np.ndarray, beta: float | np.ndarray) -> float | np.ndarray:

   Ft = -0.08958 + 0.1408 * TSR + 0.01699 * beta \
        - 0.003703 * TSR**2 - 0.007462 * TSR * beta + 0.001383 * beta**2 \
        - 0.0001281 * TSR**3 + 0.0003408 * TSR**2 * beta - 0.0004531 * TSR * beta**2 - 0.0003628 * beta**3 \
        + 1.005e-05 * TSR**4 - 5.927e-05 * TSR**3 * beta + 3.485e-05 * TSR**2 * beta**2 + 7.359e-05 * TSR * beta**3 - 2.542e-05 * beta**4 \
        - 1.403e-07 * TSR**5 + 7.613e-07 * TSR**4 * beta + 8.261e-07 * TSR**3 * beta**2 - 7.007e-06 * TSR**2 * beta**3 - 3.734e-06 * TSR * beta**4 + 5.575e-06 * beta**5;

   return Ft;

def Cp_TSR(TSR: float | np.ndarray, beta: float | np.ndarray) -> float | np.ndarray:

   T = -0.09119 + 0.03993 * TSR + 0.005038 * beta \
       + 0.02204 * TSR**2 - 0.0001244 * TSR * beta + 0.003486 * beta**2 \
       - 0.003337 * TSR**3 - 0.0002893 * TSR**2 * beta - 0.001403 * TSR * beta**2 - 0.0002117 * beta**3 \
       + 0.0001546 * TSR**4 + 5.636e-05 * TSR**3 * beta + 9.427e-05 * TSR**2 * beta**2 + 0.0001068 * TSR * beta**3 - 4.918e-05 * beta**4 \
       - 2.423e-06 * TSR**5 - 1.882e-06 * TSR**4 * beta - 3.938e-06 * TSR**3 * beta**2 - 1.416e-05 * TSR**2 * beta**3 + 4.316e-06 * TSR * beta**4 + 3.648e-06 * beta**5;

   return T;

def Ct(beta: float | np.ndarray, omega_r: float | np.ndarray, wind_speed: float) -> float | np.ndarray:
    
   TSR = omega_r * 62.94 / wind_speed;
   Ft = - 0.08958 + 0.1408 * TSR + 0.01699 * beta \
        - 0.003703 * TSR**2 - 0.007462 * TSR * beta + 0.001383 * beta**2 \
        - 0.0001281 * TSR**3 + 0.0003408 * TSR**2 * beta - 0.0004531 * TSR * beta**2 - 0.0003628 * beta**3 \
        + 1.005e-05 * TSR**4 - 5.927e-05 * TSR**3 * beta + 3.485e-05 * TSR**2 * beta**2 + 7.359e-05 * TSR * beta**3 - 2.542e-05 * beta**4 \
        - 1.403e-07 * TSR**5 + 7.613e-07 * TSR**4 * beta + 8.261e-07 * TSR**3 * beta**2 - 7.007e-06 * TSR**2 * beta**3 - 3.734e-06 * TSR * beta**4 + 5.575e-06 * beta**5;

   return Ft;

def Cp(beta: float | np.ndarray, omega_r: float | np.ndarray, wind_speed: float) -> float | np.ndarray:

   TSR = omega_r * 62.94 / wind_speed;
   T = - 0.09119 + 0.03993 * TSR + 0.005038 * beta \
       + 0.02204 * TSR**2 - 0.0001244 * TSR * beta + 0.003486 * beta**2 \
       - 0.003337 * TSR**3 - 0.0002893 * TSR**2 * beta - 0.001403 * TSR * beta**2 - 0.0002117 * beta**3 \
       + 0.0001546 * TSR**4 + 5.636e-05 * TSR**3 * beta + 9.427e-05 * TSR**2 * beta**2 + 0.0001068 * TSR * beta**3 - 4.918e-05 * beta**4 \
       - 2.423e-06 * TSR**5 - 1.882e-06 * TSR**4 * beta - 3.938e-06 * TSR**3 * beta**2 - 1.416e-05 * TSR**2 * beta**3 + 4.316e-06 * TSR * beta**4 + 3.648e-06 * beta**5;

   return T;

def WT(t: float, 
       x: list,
       wind_speed: float,
       beta_ref: float,
       T_ref: float,
       HSSB_ref: float, 
       T_yaw: float) -> list:
   
   R = 62.94;
   Jr = 55 * 10**6;
   Kdt = 2.7 * 10**9;
   Bdt = 775.49;
   Br = 7.11;
   gear_ratio = 95.0;
   dt_eff = 0.97;
   Jg = 390;
   Bg = 45.6;
   tau_g = 1 / 50;
   tau_HSSB = 0.6;
   eps = 0.6;
   wn = 11.11;
   J_yaw = 2607890;
   D_yaw = 19160000;
   K_yaw = 9028320000;
   M_tower = 347460;
   K_tower = M_tower * (2 * np.pi * 0.3240)**2;
   D_tower = 2 * 0.001 * np.sqrt(K_tower / M_tower);

   eff_wind_speed = wind_speed #- np.sqrt(x[3] + x[5]);
   TSR = x[0] * R / eff_wind_speed;
   Cp = Cp_TSR(TSR, x[10]);
   Ct = Ct_TSR(TSR, x[10]);
   rho_air = 1.25;
   T_aero = np.pi / 2 * rho_air * R**3 * Cp / (TSR + 0.00001) * eff_wind_speed**2 * np.cos(x[8]);
   Ft = np.pi / 2 * rho_air * R**2 * Ct * eff_wind_speed**2;

   domega_r = (T_aero - Kdt * x[2] - (Bdt + Br) * x[0] + Bdt / gear_ratio * x[1] - x[12]) / Jr;
   domega_g = (dt_eff * Kdt / gear_ratio * x[2] + dt_eff * Bdt / gear_ratio * x[0] - (dt_eff * Bdt / gear_ratio**2 + Bg) * x[1] - x[11]) / Jg;
   dtheta = x[0] - x[1] / gear_ratio;
   dxdot = (Ft * np.cos(x[8]) - D_tower * x[3] - K_tower * x[4]) / M_tower;
   dx = x[3];
   dydot = (Ft * np.sin(x[8]) - D_tower * x[4] - K_tower * x[5]) / M_tower;
   dy = x[5];
   dgamadot = (T_yaw - D_yaw * x[7] - K_yaw * x[8]) / J_yaw;
   dgama = x[7];
   dbetadot = beta_ref - 2 * eps * wn * x[9] - wn**2 * x[10];
   dbeta = x[9];
   dTg = (T_ref - x[11]) / tau_g;
   dTHSSB = (HSSB_ref - x[12]) / tau_HSSB;

   return [domega_r, domega_g, dtheta, dxdot, dx, dydot, dy, dgamadot, dgama, dbetadot, dbeta, dTg, dTHSSB];

def controller(x):

   return [0.0, 0.0, 0.0, 0.0];

if __name__ == "__main__":

   # Plot aerodynamic power and thrust force surfaces.
    
   tsr = np.arange(0, 25, 0.1);
   beta = np.arange(-5, 9, 0.1);
   tsr, beta = np.meshgrid(tsr, beta);
   Cp_surf = Cp_TSR(tsr, beta);
   Ct_surf = Ct_TSR(tsr, beta);

   fig, ax = plt.subplots(2, 1, subplot_kw={"projection": "3d"});

   ax[0].plot_surface(tsr, beta, Cp_surf);
   ax[0].set_title("Cp (TSR, beta)");
   ax[0].set_xlabel("TSR");
   ax[0].set_ylabel("beta");
   ax[0].set_zlabel("Cp");
   ax[0].view_init(elev=25, azim=45, roll=0);

   ax[1].plot_surface(tsr, beta, Ct_surf);
   ax[1].set_title("Ct (TSR, beta)");
   ax[1].set_xlabel("TSR");
   ax[1].set_ylabel("beta");
   ax[1].set_zlabel("Ct");
   ax[1].view_init(elev=25, azim=45, roll=0);

   plt.show();

   # Initial conditions
   x0 = np.zeros(13).tolist();

   # Time span
   t_span = [0, 10];
   Ts = 0.1;
   n = int((t_span[1] - t_span[0]) / Ts + 1);

   if(int((t_span[1] - t_span[0]) % Ts) != 0):
    
     n += 1;
    
   true_time = np.array([np.fmin(t_span[1], t_span[0] + Ts * i) for i in range(n)]);

   results = np.zeros((14, 1));

   for t in range(n-1):

        # Control input
        beta_ref, T_ref, T_yaw, T_HSSB = controller(x0);

        # Solve ODE
        sol = solve_ivp(WT, [true_time[t], true_time[t+1]], x0, args=(10, beta_ref, T_ref, T_yaw, T_HSSB));

        # Store results
        new_data = np.vstack((sol.t, sol.y));
        results = np.hstack((results, new_data));

        # Update initial conditions
        x0 = sol.y[:, -1];

   # Plot results
   fig = plt.figure(figsize=(5, 5));

   ax = fig.subplots(6, 1, sharex=True);

   ax[0].set_title("Wind Turbine response");
   ax[0].set_xlabel("t (s)");
   ax[0].set_ylabel("Gen. speed (rad/s)");
   ax[0].plot(results[0, :], results[2, :] * 60 / (2 * np.pi));
   ax[0].grid();

   ax[1].set_title("Wind Turbine response");
   ax[1].set_xlabel("t (s)");
   ax[1].set_ylabel("Rotor speed (rad/s)");
   ax[1].plot(results[0, :], results[1, :] * 60 / (2 * np.pi));
   ax[1].grid();

   ax[2].set_title("Wind Turbine response");
   ax[2].set_xlabel("t (s)");
   ax[2].set_ylabel("Torsion angle (rad)");
   ax[2].plot(results[0, :], results[3, :] / (2 * np.pi));
   ax[2].grid();

   ax[3].set_title("Wind Turbine response");
   ax[3].set_xlabel("t (s)");
   ax[3].set_ylabel("Beta (deg)");
   ax[3].plot(results[0, :], results[11, :]);
   ax[3].grid();

   ax[4].set_title("Wind Turbine response");
   ax[4].set_xlabel("t (s)");
   ax[4].set_ylabel("Gama (deg)");
   ax[4].plot(results[0, :], results[9, :]);
   ax[4].grid();

   ax[5].set_title("Wind Turbine response");
   ax[5].set_xlabel("t (s)");
   ax[5].set_ylabel("x_dot (m/s)");
   ax[5].plot(results[0, :], results[5, :]);
   ax[5].grid();

   plt.show();