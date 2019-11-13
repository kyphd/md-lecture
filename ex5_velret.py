# ex5_velret
import numpy as np

def velret(t_start, t_end, dt, x_0, u_0, y_0, v_0):

    t = t_start
    i = 0

    while t < t_end:

        t += dt
        i += 1

        r3 = np.sqrt(x_0 ** 2 + y_0 ** 2) ** 3
        u_tmp = u_0 + dt / 2.0 * -x_0 / r3
        v_tmp = v_0 + dt / 2.0 * -y_0 / r3

        x_1 = x_0 + u_tmp * dt
        y_1 = y_0 + v_tmp * dt

        r3_1 = np.sqrt(x_1 ** 2 + y_1 ** 2) ** 3
        u_1 = u_tmp + dt / 2.0 * -x_1 / r3_1
        v_1 = v_tmp + dt / 2.0 * -y_1 / r3_1

        x_0 = x_1
        u_0 = u_1
        y_0 = y_1
        v_0 = v_1

        energy = 1.0 / 2.0 * (u_0 ** 2 + v_0 ** 2) - 1.0 / np.sqrt(x_0 ** 2 + y_0 ** 2)

        if i % 1000 == 0:
            print("{0:10f} {1:5.4f} {2:5.4f} {3:10.8f}".format(t, x_0, y_0, energy))


# main
if __name__ == "__main__":
    t = 0.0
    t_end = 300000.0
    x_0 = 1.0
    u_0 = 0.0
    y_0 = 0.0
    v_0 = 1.0

    velret(t, t_end, 0.1, x_0, u_0, y_0, v_0)
