#!/bin/bash

N=$1

python generate_matrix.py $N

g++ serial.cpp shared.cpp -o serial
serial_time=$(./serial $N)

echo "Serial time: $serial_time"

mpic++ parallel.cpp shared.cpp -o parallel

declare -a parallel_times

for ((num_procs=2; num_procs<=8; num_procs++))
do
    echo "Running with $num_procs processes:"
    time_output=$(mpiexec -n $num_procs ./parallel $N)
    echo "Time: $time_output"
    parallel_times+=("$time_output")
done

python plot_speedup.py $serial_time "${parallel_times[@]}"
