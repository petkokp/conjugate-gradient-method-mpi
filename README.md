# Conjugate gradient method with MPI

N - size of matrix

## Generate matrix:

python generate_matrix.py N

#

## Serial:

g++ serial.cpp shared.cpp -o serial

./serial N

#

## Parallel:

mpic++ parallel.cpp shared.cpp -o parallel

mpiexec -n 2 ./parallel N