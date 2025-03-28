import numpy as np
import matplotlib.pyplot as plt

def fan_power_curve(x: float, a: float, b: float, c: float) -> float:

    return a * x**2 + b * x + c;

def pump_power_curve(x: float, a: float, b: float, c: float, d: float) -> float:

    return a * x**3 + b * x**2 + c * x + d;

def print_params(opt_param: np.ndarray) -> None:

    print("Optimal curve");
    print("-" * 13);
    print(f"{opt_param[0]:.3f} * x**2 + {opt_param[1]:.3f} * x + {opt_param[2]:.3f}");

    return;

def print_pump_params(opt_param: np.ndarray) -> None:

    print("Optimal curve");
    print("-" * 13);
    print(f"{opt_param[0]:.3f} * x**3 + {opt_param[1]:.3f} * x**2 + {opt_param[2]:.3f} * x + {opt_param[3]:.3f}");

    return;

def plot_fit(xdata: np.ndarray, ydata: np.ndarray, func: callable, opt_param: np.ndarray) -> None:

    fig = plt.figure(figsize=(5, 5));

    ax = fig.add_subplot(1, 1, 1);

    ax.set_title("Curve fit");
    ax.scatter(xdata, ydata, label="Data", c="b");
    ax.plot(xdata, func(xdata, *opt_param), label="Fit", c="r");
    ax.grid();
    ax.legend();
    
    plt.show();

    return;

if __name__ == "__main__":

    from scipy.optimize import curve_fit

    xdata_fan = np.array([2, 42, 85, 127, 170, 212, 255, 297, 340]) / 3600;
    ydata_fan = np.array([170, 220, 290, 340, 390, 420, 430, 430, 430]);

    opt_param, _ = curve_fit(fan_power_curve, xdata_fan, ydata_fan);
    print_params(opt_param);
    plot_fit(xdata_fan, ydata_fan, fan_power_curve, opt_param);

    xdata_pump = np.arange(1, 10, 1) / 3600;
    ydata_pump = np.array([200, 220, 245, 270, 290, 305, 315, 310, 305]);

    opt_param_pump, _ = curve_fit(pump_power_curve, xdata_pump, ydata_pump);
    print_pump_params(opt_param_pump);
    plot_fit(xdata_pump, ydata_pump, pump_power_curve, opt_param_pump);