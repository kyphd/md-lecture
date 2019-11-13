# ex6_nve
import numpy as np
import random
import matplotlib.pyplot as plt
import matplotlib.animation as animation


# Constant
Kb = 8.61e-5        # Boltzmann constant (eV/K)
Na = 6.02e23        # Avogadro constant (/mol)

# Nondimensionalization Constant
unit_length = 3.4                       # Length (A)
unit_mass = 39.95 / Na * 1e-3           # Mass (kg)    (for Ar)
unit_temp = 120.0                       # Temperature (K)     (for Ar)
unit_eng = unit_temp * Kb               # Energy (eV)
unit_time = np.sqrt(unit_mass * unit_length ** 2 / unit_eng)  # Time (s)
unit_pressure = unit_eng * 1.602e-19 / (unit_length * 1e-10) ** 2 / 1e-10    # Pressure (Pa) for 2D

# cell
lat_const = 3.4     # Lattice Constant (A)
Lx = 5              # Number of x lattice
Ly = 5              # Number of y lattice
Natom = Lx * Ly     # Number of Atoms
dt = 2e-15          # Time step (s)
Max_step = 2000     # Calc step
Tinitial = 100      # Initial Temperature [K]
rcut_ = 3.0         # Cut off (Nondimension)

# Nondimensionalization
m_ = 1.0                                # Mass
Tinit_ = Tinitial / unit_temp           # Temperature
dt_ = dt / unit_time                    # Time step
lat_ = lat_const / unit_length          # Lattice Constant
cellsize_x_ = lat_ * Lx                 # Cell size : x
cellsize_y_ = lat_ * Ly                 # Cell size : y
volume_ = cellsize_x_ * cellsize_y_     # Volume : 2D

# for plot
fig = plt.figure(figsize=(4, 4))
ims = []


# Velocity of atoms: Maxwell distribution (Box-Muller)
def boxmuller(vx_, vy_, T_):
    for i in range(Natom):
        R1 = random.random()
        R2 = random.random()
        R3 = random.random()
        R4 = random.random()
        vx_[i] = (np.sqrt(-2 * (T_ / m_) * np.log(R1))) * np.cos(2 * np.pi * R2)
        vy_[i] = (np.sqrt(-2 * (T_ / m_) * np.log(R3))) * np.cos(2 * np.pi * R4)


# Initial position and velocity of atoms
def initialize_pos_vel(x_, y_, vx_, vy_, T_):
    i = -1
    for ix in range(Lx):
        for iy in range(Ly):
            i += 1
            x_[i] = lat_ * ix
            y_[i] = lat_ * iy

    boxmuller(vx_, vy_, T_)


# Calculate Force
def forces(x_, y_, fx_, fy_):

    PE_ = 0
    virial_ = 0     # Virial

    rcut_2 = rcut_ ** 2

    fx_[:] = 0
    fy_[:] = 0
    for i in range(Natom - 1):
        for j in range(i + 1, Natom):
            dx_ = x_[i] - x_[j]
            dy_ = y_[i] - y_[j]

            # periodic image of particle?
            if abs(dx_) > 0.5 * cellsize_x_:
                dx_ = dx_ - np.sign(dx_) * cellsize_x_
            if abs(dy_) > 0.5 * cellsize_y_:
                dy_ = dy_ - np.sign(dy_) * cellsize_y_

            r_2 = dx_ ** 2 + dy_ ** 2
            if r_2 < rcut_2:
                if r_2 == 0:
                    r_2 = 0.0001
                invr_2 = 1.0 / r_2
                invr_6 = invr_2**3
                invr_12 = invr_6**2

                virialij_ = 48 * (invr_12 - 0.5 * invr_6)
                fijx_ = virialij_ * invr_2 * dx_
                fijy_ = virialij_ * invr_2 * dy_

                fx_[i] += fijx_
                fy_[i] += fijy_
                fx_[j] -= fijx_
                fy_[j] -= fijy_

                PE_ += 4 * (invr_12 - invr_6)
                virial_ += virialij_

    return PE_, virial_


# Calculate pressure
def pressure(n, volume, T, virial):
    return n / volume * T + virial / (2 * volume)


# main: time evolution
def nve():

    # Position
    x_ = np.zeros([Natom], float)
    y_ = np.zeros([Natom], float)

    # Velocity
    vx_ = np.zeros([Natom], float)
    vy_ = np.zeros([Natom], float)

    # Force
    fx_ = np.zeros([Natom], float)
    fy_ = np.zeros([Natom], float)

    n_step = 0
    while n_step <= Max_step:

        if n_step == 0:

            # Initial position and velocity
            initialize_pos_vel(x_, y_, vx_, vy_, Tinit_)
            PE_, virial_ = forces(x_, y_, fx_, fy_)  # Calculate f_0

        else:
            # Velocity Velret
            # update velocity'
            vx_ += dt_ * fx_ / 2
            vy_ += dt_ * fy_ / 2

            # update position for all atoms
            for i in range(Natom):

                # update position
                x_[i] += dt_ * vx_[i]
                y_[i] += dt_ * vy_[i]

                # periodic boundary
                if x_[i] < 0:
                    x_[i] += cellsize_x_
                if x_[i] >= cellsize_x_:
                    x_[i] -= cellsize_x_
                if y_[i] < 0:
                    y_[i] += cellsize_y_
                if y_[i] >= cellsize_y_:
                    y_[i] -= cellsize_y_

            # Calculate f_n+1
            PE_, virial_ = forces(x_, y_, fx_, fy_)

            # Calculate v_n+1
            vx_ += dt_ * fx_ / 2
            vy_ += dt_ * fy_ / 2

        # Calculate Energy, Temperature, Pressure
        KE_ = (np.sum(vx_ ** 2) + np.sum(vy_ ** 2)) / 2
        T_ = KE_ / Natom
        P_ = pressure(Natom, volume_, T_, virial_)

        # 5ステップ毎にプロット
        if n_step % 5 == 0:
            im = plt.scatter(x_, y_, s=120, c="red", alpha=0.9, linewidths="1.5", edgecolors="black")
            ims.append([im])

            # Convert values
            T = T_ * unit_temp  # K
            P = P_ * unit_pressure  # Pa
            KE = KE_ * unit_eng / Natom  # eV/atom
            PE = PE_ * unit_eng / Natom  # eV/atom
            Etot = KE + PE  # eV/atom

            print("{0:5d} {1:10.5f} {2:10.5f} {3:10.5f} {4:10.5f} {5:10.5f}".format(n_step, T, P, KE, PE, Etot))

        n_step += 1


# メイン
if __name__ == "__main__":
    nve()

    # プロット+動画の保存
    plt.axis([0, cellsize_x_, 0, cellsize_y_])
    ani = animation.ArtistAnimation(fig, ims, interval=1)
    plt.show()
    ani.save('anim.mp4', writer="ffmpeg", fps=30)
