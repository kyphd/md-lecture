# ex1_euler
import numpy as np
import matplotlib.pyplot as plt

# initial
t = 0
y_0 = 1
dt = 0.1        # change dt to compare

# euler
t_ary = [t]
y_ary = [y_0]

while t < 5:
    y_1 = y_0 + y_0 * dt
    t += dt

    t_ary.append(t)
    y_ary.append(y_1)

    y_0 = y_1


# analytic
t_ary_analytic = np.arange(0, 5.1, 0.1)
y_ary_analytic = np.exp(t_ary_analytic)

# plot
plt.plot(t_ary, y_ary, label="euler")
plt.plot(t_ary_analytic, y_ary_analytic, label="analytic")
plt.legend()
plt.show()
