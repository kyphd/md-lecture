# ex4_kidou_longtime
import numpy as np

def rungekutta(t_start, t_end, dt, x_0, u_0, y_0, v_0):

    t = t_start
    i = 0

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
        i += 1

        x_0 = x_1
        u_0 = u_1
        y_0 = y_1
        v_0 = v_1

        energy = 1.0 / 2.0 * (u_0**2 + v_0**2) - 1.0 / np.sqrt(x_0**2 + y_0**2)

        if i % 1000 == 0:
            print("{0:10f} {1:5.4f} {2:5.4f} {3:5.4f}".format(t, x_0, y_0, energy))


# main
if __name__ == "__main__":
    t = 0.0
    t_end = 300000.0
    x_0 = 1.0
    u_0 = 0.0
    y_0 = 0.0
    v_0 = 1.0

    rungekutta(t, t_end, 0.1, x_0, u_0, y_0, v_0)