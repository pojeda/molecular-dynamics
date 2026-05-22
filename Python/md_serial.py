from numba import njit
import math
import numpy as np
import time

# Parameters
NUM_RES = 1000
NSTEPS = 10000
WRITE_FREQ = 1000

DT = 0.001
DT2 = DT * DT
DTP5 = 0.5 * DT

EPS_CONST = 40.0
SIGMA_CONST = 6.5
R_CERO = 3.8
A_CONST = 50.0

OLD = 0
NEW = 1

# Arrays
coor = np.zeros((NUM_RES, 3), dtype=np.float64)
vel = np.zeros((NUM_RES, 3), dtype=np.float64)
force = np.zeros((2, NUM_RES, 3), dtype=np.float64)
mass = np.ones(NUM_RES, dtype=np.float64)

pot_ener = 0.0
kin_ener = 0.0


def initialize():
    global pot_ener, kin_ener

    for p in range(NUM_RES):
        coor[p] = [(p + 1 - 15) * 4.0, 0.0, 0.0]
        vel[p] = np.random.normal(0.0, 1.0, 3)

    force[:, :, :] = 0.0
    mass[:] = 1.0

    pot_ener = 0.0
    kin_ener = 0.0


@njit
def compute_forces_numba(coor, force_new):
    pot = 0.0
    force_new[:, :] = 0.0

    # Spring forces
    for i in range(NUM_RES - 1):
        dx = coor[i, 0] - coor[i + 1, 0]
        dy = coor[i, 1] - coor[i + 1, 1]
        dz = coor[i, 2] - coor[i + 1, 2]

        r = math.sqrt(dx*dx + dy*dy + dz*dz)
        dr = r - R_CERO

        pot += A_CONST * dr * dr

        f = -2.0 * A_CONST * dr / r

        fx = f * dx
        fy = f * dy
        fz = f * dz

        force_new[i, 0] += fx
        force_new[i, 1] += fy
        force_new[i, 2] += fz

        force_new[i + 1, 0] -= fx
        force_new[i + 1, 1] -= fy
        force_new[i + 1, 2] -= fz

    sigma2 = SIGMA_CONST * SIGMA_CONST

    # Lennard-Jones forces
    for i in range(NUM_RES - 2):
        for j in range(i + 2, NUM_RES):
            dx = coor[i, 0] - coor[j, 0]
            dy = coor[i, 1] - coor[j, 1]
            dz = coor[i, 2] - coor[j, 2]

            r2 = dx*dx + dy*dy + dz*dz

            sr2 = sigma2 / r2
            sr4 = sr2 * sr2
            sr6 = sr4 * sr2
            sr12 = sr6 * sr6

            pot += EPS_CONST * (sr12 - sr6)

            fterm = EPS_CONST * (12.0 * sr12 - 6.0 * sr6) / r2

            fx = fterm * dx
            fy = fterm * dy
            fz = fterm * dz

            force_new[i, 0] += fx
            force_new[i, 1] += fy
            force_new[i, 2] += fz

            force_new[j, 0] -= fx
            force_new[j, 1] -= fy
            force_new[j, 2] -= fz

    return pot

def compute_forces():
    global pot_ener, kin_ener

    pot_ener = compute_forces_numba(coor, force[NEW])
    kin_ener = 0.0


def update_coordinates():
    coor[:] += DT * vel + 0.5 * DT2 * force[OLD] / mass[:, None]


def update_velocities():
    vel[:] += DTP5 * (force[NEW] + force[OLD]) / mass[:, None]


def kinetic_energy():
    global kin_ener

    kin_ener = 0.5 * np.sum(mass[:, None] * vel * vel)


def write_xyz(file):
    file.write(f"{NUM_RES}\n")
    file.write("\n")

    for p in range(NUM_RES):
        x, y, z = coor[p]
        file.write(f"C {x:.8f} {y:.8f} {z:.8f}\n")


def main():
    global OLD, NEW

    np.random.seed()

    initialize()
    compute_forces()

    start = time.time()

    with open("KOORDINATEN_1T.xyz", "w") as f:
        for nstep in range(1, NSTEPS + 1):
            update_coordinates()
            compute_forces()
            update_velocities()
            kinetic_energy()

            if nstep % WRITE_FREQ == 0:
                write_xyz(f)
                print(pot_ener, kin_ener, pot_ener + kin_ener)

            OLD, NEW = NEW, OLD

    end = time.time()
    print(f"SERIAL TIME = {end - start:.6f} SECONDS")


if __name__ == "__main__":
    main()
