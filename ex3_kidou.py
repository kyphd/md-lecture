# ex3_kidou
import numpy as np
import matplotlib.pyplot as plt

def euler(t_start, t_end, dt, x_0, u_0, y_0, v_0, ax):

    t = t_start

    tary = [t]
    xary = [x_0]
    uary = [u_0]
    yary = [y_0]
    vary = [v_0]

    while t < t_end:
        r3 = np.sqrt(x_0 ** 2 + y_0 ** 2) ** 3

        x_1 = x_0 + u_0 * dt
        u_1 = u_0 - x_0 / r3 * dt

        y_1 = y_0 + v_0 * dt
        v_1 = v_0 - y_0 / r3 * dt

        t += dt

        tary.append(t)
        xary.append(x_1)
        uary.append(u_1)
        yary.append(y_1)
        vary.append(v_1)

        x_0 = x_1
        u_0 = u_1
        y_0 = y_1
        v_0 = v_1

    ax.scatter(xary, yary, label="euler dt=" + str(dt), s=1)
    ax.legend(loc="upper right", fontsize=8)


def rungekutta(t_start, t_end, dt, x_0, u_0, y_0, v_0, ax):

    t = t_start

    tary = [t]
    xary = [x_0]
    uary = [u_0]
    yary = [y_0]
    vary = [v_0]

    while t < t_end:

        r3 = np.sqrt(x_0 ** 2 + y_0 ** 2) ** 3

        x_k1 = u_0
        u_k1 = -x_0 / r3
        y_k1 = v_0
        v_k1 = -y_0 / r3

        r3_k1 = np.sqrt((x_0 + x_k1 * dt/2.0) ** 2 + (y_0 + y_k1 * dt/2.0) ** 2) ** 3

        x_k2 = u_0 + dt / 2.0 * u_k1
        u_k2 = - (x_0 + x_k1 * dt / 2.0) / r3_k1
        y_k2 = v_0 + dt / 2.0 * v_k1
        v_k2 = - (y_0 + y_k1 * dt / 2.0) / r3_k1

        r3_k2 = np.sqrt((x_0 + x_k2 * dt/2.0) ** 2 + (y_0 + y_k2 * dt/2.0) ** 2) ** 3

        x_k3 = u_0 + dt / 2.0 * u_k2
        u_k3 = - (x_0 + x_k2 * dt / 2.0) / r3_k2
        y_k3 = v_0 + dt / 2.0 * v_k2
        v_k3 = - (y_0 + y_k2 * dt / 2.0) / r3_k2

        r3_k3 = np.sqrt((x_0 + x_k3 * dt) ** 2 + (y_0 + y_k3 * dt) ** 2) ** 3

        x_k4 = u_0 + dt * u_k3
        u_k4 = - (x_0 + x_k3 * dt) / r3_k3
        y_k4 = v_0 + dt * v_k3
        v_k4 = - (y_0 + y_k3 * dt) / r3_k3

        x_1 = x_0 + dt / 6.0 * (x_k1 + 2 * x_k2 + 2 * x_k3 + x_k4)
        u_1 = u_0 + dt / 6.0 * (u_k1 + 2 * u_k2 + 2 * u_k3 + u_k4)
        y_1 = y_0 + dt / 6.0 * (y_k1 + 2 * y_k2 + 2 * y_k3 + y_k4)
        v_1 = v_0 + dt / 6.0 * (v_k1 + 2 * v_k2 + 2 * v_k3 + v_k4)

        t += dt

        tary.append(t)
        xary.append(x_1)
        uary.append(u_1)
        yary.append(y_1)
        vary.append(v_1)

        x_0 = x_1
        u_0 = u_1
        y_0 = y_1
        v_0 = v_1

    ax.scatter(xary, yary, label="RungeKutta dt=" + str(dt), s=1)
    ax.legend(loc="upper right", fontsize=8)


# main
if __name__ == "__main__":

    # initial condition
    t = 0.0
    t_end = 10.0
    x_0 = 1.0
    u_0 = 0.0
    y_0 = 0.0
    v_0 = 1.0

    # for plot
    fig = plt.figure(figsize=(8, 4))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    ax1.set_aspect('equal')
    ax2.set_aspect('equal')

    ax1.set_xlim(-2, 2)
    ax1.set_ylim(-2, 2)
    ax2.set_xlim(-2, 2)
    ax2.set_ylim(-2, 2)

    # calculation
    euler(t, t_end, 0.1, x_0, u_0, y_0, v_0, ax1)
    euler(t, t_end, 0.01, x_0, u_0, y_0, v_0, ax1)
    euler(t, t_end, 0.001, x_0, u_0, y_0, v_0, ax1)
    rungekutta(t, t_end * 10, 0.5, x_0, u_0, y_0, v_0, ax2)
    rungekutta(t, t_end * 10, 0.4, x_0, u_0, y_0, v_0, ax2)
    rungekutta(t, t_end * 10, 0.1, x_0, u_0, y_0, v_0, ax2)

    # plot
    plt.show()
