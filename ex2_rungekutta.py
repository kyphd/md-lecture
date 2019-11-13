# ex2_rungekutta
import numpy as np
import matplotlib.pyplot as plt

# initial
t = 0
y_0 = 1
dt = 0.1        # change dt to compare

t_ary = [t]
y_ary = [y_0]

while t < 5:

    k1 = y_0
    k2 = y_0 + dt / 2.0 * k1
    k3 = y_0 + dt / 2.0 * k2
    k4 = y_0 + dt * k3

    y_1 = y_0 + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
    t += dt

    t_ary.append(t)
    y_ary.append(y_1)

    y_0 = y_1


# 解析解
t_ary_analytic = np.arange(0, 5.1, 0.1)
y_ary_analytic = np.exp(t_ary_analytic)

plt.plot(t_ary, y_ary, label="RungeKutta")
plt.plot(t_ary_analytic, y_ary_analytic, label="Analytic")
plt.legend()
plt.show()
