from numba import njit
import math
import numpy as np
import time
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

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

@njit
def compute_local_forces_numba(coor, force_slot, istart, iend):
    pot = 0.0
    force_slot[istart:iend, :] = 0.0

    # Spring force
    for i in range(istart, iend):

        if i < NUM_RES - 1:
            dx = coor[i, 0] - coor[i + 1, 0]
            dy = coor[i, 1] - coor[i + 1, 1]
            dz = coor[i, 2] - coor[i + 1, 2]

            r = math.sqrt(dx*dx + dy*dy + dz*dz)
            dr = r - R_CERO

            f = -2.0 * A_CONST * dr / r

            force_slot[i, 0] += f * dx
            force_slot[i, 1] += f * dy
            force_slot[i, 2] += f * dz

            pot += A_CONST * dr * dr

        if i > 0:
            dx = coor[i, 0] - coor[i - 1, 0]
            dy = coor[i, 1] - coor[i - 1, 1]
            dz = coor[i, 2] - coor[i - 1, 2]

            r = math.sqrt(dx*dx + dy*dy + dz*dz)
            dr = r - R_CERO

            f = -2.0 * A_CONST * dr / r

            force_slot[i, 0] += f * dx
            force_slot[i, 1] += f * dy
            force_slot[i, 2] += f * dz

    sigma2 = SIGMA_CONST * SIGMA_CONST

    # Lennard-Jones force
    for i in range(istart, iend):
        for j in range(NUM_RES):

            if j == i:
                continue
            if abs(j - i) <= 1:
                continue

            dx = coor[i, 0] - coor[j, 0]
            dy = coor[i, 1] - coor[j, 1]
            dz = coor[i, 2] - coor[j, 2]

            r2 = dx*dx + dy*dy + dz*dz

            sr2 = sigma2 / r2
            sr4 = sr2 * sr2
            sr6 = sr4 * sr2
            sr12 = sr6 * sr6

            fterm = EPS_CONST * (12.0 * sr12 - 6.0 * sr6) / r2

            force_slot[i, 0] += fterm * dx
            force_slot[i, 1] += fterm * dy
            force_slot[i, 2] += fterm * dz

            if j > i + 1:
                pot += EPS_CONST * (sr12 - sr6)

    return pot

def compute_forces(slot):
    global pot_ener, kin_ener

    pot_ener = compute_local_forces_numba(coor, force[slot], istart, iend)
    kin_ener = 0.0

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

@njit
def kinetic_energy_numba(vel, mass, istart, iend):
    kin = 0.0

    for i in range(istart, iend):
        kin += mass[i] * (
            vel[i, 0] * vel[i, 0]
            + vel[i, 1] * vel[i, 1]
            + vel[i, 2] * vel[i, 2]
        )

    return 0.5 * kin

def kinetic_energy():
    global kin_ener
    kin_ener = kinetic_energy_numba(vel, mass, istart, iend)

def exchange_coordinates_collective():
    sendbuf = np.ascontiguousarray(coor[istart:iend].ravel())

    comm.Allgatherv(
        sendbuf,
        [coor.ravel(), counts, displs, MPI.DOUBLE]
    )






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
