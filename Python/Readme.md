ml GCCcore/11.2.0 Python/3.9.6
python -m venv vpyenv-python-course
source vpyenv-python-course/bin/activate  
ml GCC/11.2.0 OpenMPI/4.1.1
pip install mpi4py 
pip install numpy 
pip install numba

python md_serial.py
mpirun -np 10 md_mpi_p2p.py
mpirun -np 10 md_mpi_collectives.py
