import sys
import matplotlib.pyplot as plt

def calculate_speedup(serial_time, parallel_times):
    return [serial_time / t for t in parallel_times]

def plot_speedup(num_processes, speedup_values):
    plt.figure(figsize=(8, 6))
    plt.plot(num_processes, speedup_values, marker='o')
    plt.title('Parallel Speedup vs. Number of Processes')
    plt.xlabel('Number of Processes')
    plt.ylabel('Speedup')
    plt.grid(True)
    plt.xticks(num_processes)
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 9:
        print("Usage: python calculate_speedup.py serial_time (parallel_times)")
        sys.exit(1)

    serial_time = float(sys.argv[1])
    parallel_times = [float(time) for time in sys.argv[2:]]

    num_processes = list(range(2, 9))
    speedup_values = calculate_speedup(serial_time, parallel_times)
    plot_speedup(num_processes, speedup_values)
