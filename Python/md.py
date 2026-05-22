import numpy as np
import time

# Parameters
NUM_RES = 1000
NSTEPS = 1000
WRITE_FREQ = 100

PI = np.pi
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


def spring_force():
    global pot_ener

    for i in range(NUM_RES - 1):
        rij = coor[i] - coor[i + 1]
        r = np.linalg.norm(rij)
        dr = r - R_CERO

        pot_ener += A_CONST * dr * dr

        fij = -2.0 * A_CONST * dr * rij / r

        force[NEW, i] += fij
        force[NEW, i + 1] -= fij


def lennard_jones_force():
    global pot_ener

    sigma2 = SIGMA_CONST * SIGMA_CONST

    for i in range(NUM_RES - 2):
        for j in range(i + 2, NUM_RES):
            rij = coor[i] - coor[j]
            r2 = np.dot(rij, rij)

            sr2 = sigma2 / r2
            sr4 = sr2 * sr2
            sr6 = sr4 * sr2
            sr12 = sr6 * sr6

            eij = EPS_CONST * (sr12 - sr6)
            fterm = EPS_CONST * (12.0 * sr12 - 6.0 * sr6) / r2

            pot_ener += eij

            force[NEW, i] += fterm * rij
            force[NEW, j] -= fterm * rij


def compute_forces():
    global pot_ener, kin_ener

    force[NEW, :, :] = 0.0
    pot_ener = 0.0
    kin_ener = 0.0

    spring_force()
    lennard_jones_force()


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
