import numpy as np
import time
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

NUM_RES = 1000
NSTEPS = 1000
WRITE_FREQ = 100

DT = 0.001
DT2 = DT * DT
DTP5 = 0.5 * DT

EPS_CONST = 40.0
SIGMA_CONST = 6.5
R_CERO = 3.8
A_CONST = 50.0

OLD = 0
NEW = 1

coor = np.zeros((NUM_RES, 3), dtype=np.float64)
vel = np.zeros((NUM_RES, 3), dtype=np.float64)
force = np.zeros((2, NUM_RES, 3), dtype=np.float64)
mass = np.ones(NUM_RES, dtype=np.float64)

pot_ener = 0.0
kin_ener = 0.0


def get_range(n, rank, size):
    base = n // size
    rem = n % size

    if rank < rem:
        start = rank * (base + 1)
        end = start + base + 1
    else:
        start = rem * (base + 1) + (rank - rem) * base
        end = start + base

    return start, end


istart, iend = get_range(NUM_RES, rank, size)
nlocal = iend - istart


counts = np.array(
    [3 * (get_range(NUM_RES, r, size)[1] - get_range(NUM_RES, r, size)[0])
     for r in range(size)],
    dtype=np.int32
)

displs = np.array(
    [3 * get_range(NUM_RES, r, size)[0]
     for r in range(size)],
    dtype=np.int32
)


def initialize():
    global pot_ener, kin_ener

    rng = np.random.default_rng(seed=int(time.time()) + 104729 * rank)

    coor[:, :] = 0.0
    vel[:, :] = 0.0
    force[:, :, :] = 0.0
    mass[:] = 1.0

    for p in range(istart, iend):
        coor[p] = [(p + 1 - 15) * 4.0, 0.0, 0.0]
        vel[p] = rng.normal(0.0, 1.0, 3)

    pot_ener = 0.0
    kin_ener = 0.0


def exchange_coordinates_collective():
    sendbuf = np.ascontiguousarray(coor[istart:iend].ravel())

    comm.Allgatherv(
        sendbuf,
        [coor.ravel(), counts, displs, MPI.DOUBLE]
    )


def spring_force(slot):
    global pot_ener

    for i in range(istart, iend):

        if i < NUM_RES - 1:
            rij = coor[i] - coor[i + 1]
            r = np.linalg.norm(rij)
            dr = r - R_CERO

            fij = -2.0 * A_CONST * dr * rij / r
            force[slot, i] += fij

            pot_ener += A_CONST * dr * dr

        if i > 0:
            rij = coor[i] - coor[i - 1]
            r = np.linalg.norm(rij)
            dr = r - R_CERO

            fij = -2.0 * A_CONST * dr * rij / r
            force[slot, i] += fij


def lennard_jones_force(slot):
    global pot_ener

    sigma2 = SIGMA_CONST * SIGMA_CONST

    for i in range(istart, iend):
        for j in range(NUM_RES):

            if j == i:
                continue
            if abs(j - i) <= 1:
                continue

            rij = coor[i] - coor[j]
            r2 = np.dot(rij, rij)

            sr2 = sigma2 / r2
            sr4 = sr2 * sr2
            sr6 = sr4 * sr2
            sr12 = sr6 * sr6

            eij = EPS_CONST * (sr12 - sr6)
            fterm = EPS_CONST * (12.0 * sr12 - 6.0 * sr6) / r2

            force[slot, i] += fterm * rij

            if j > i + 1:
                pot_ener += eij


def compute_forces(slot):
    global pot_ener, kin_ener

    force[slot, istart:iend, :] = 0.0
    pot_ener = 0.0
    kin_ener = 0.0

    spring_force(slot)
    lennard_jones_force(slot)


def update_coordinates():
    coor[istart:iend] += (
        DT * vel[istart:iend]
        + 0.5 * DT2 * force[OLD, istart:iend] / mass[istart:iend, None]
    )


def update_velocities():
    vel[istart:iend] += (
        DTP5
        * (force[NEW, istart:iend] + force[OLD, istart:iend])
        / mass[istart:iend, None]
    )


def kinetic_energy():
    global kin_ener

    local_vel = vel[istart:iend]
    local_mass = mass[istart:iend]

    kin_ener = 0.5 * np.sum(local_mass[:, None] * local_vel * local_vel)


def sum_energies_collective():
    local = np.array([pot_ener, kin_ener], dtype=np.float64)
    total = np.zeros(2, dtype=np.float64)

    comm.Reduce(local, total, op=MPI.SUM, root=0)

    if rank == 0:
        return total[0], total[1]

    return None, None


def write_xyz(file):
    file.write(f"{NUM_RES}\n")
    file.write("\n")

    for p in range(NUM_RES):
        x, y, z = coor[p]
        file.write(f"C {x:.8f} {y:.8f} {z:.8f}\n")


def main():
    global OLD, NEW

    initialize()

    exchange_coordinates_collective()
    compute_forces(OLD)

    if rank == 0:
        start = time.time()
        f = open("KOORDINATEN_1T.xyz", "w")
    else:
        start = None
        f = None

    for nstep in range(1, NSTEPS + 1):

        update_coordinates()
        exchange_coordinates_collective()

        compute_forces(NEW)
        update_velocities()
        kinetic_energy()

        if nstep % WRITE_FREQ == 0:
            pot_total, kin_total = sum_energies_collective()

            if rank == 0:
                write_xyz(f)
                print(pot_total, kin_total, pot_total + kin_total)

        OLD, NEW = NEW, OLD

    if rank == 0:
        end = time.time()
        f.close()
        print(f"PARALLEL COLLECTIVE TIME = {end - start:.6f} SECONDS")


if __name__ == "__main__":
    main()
