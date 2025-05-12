import matplotlib.pyplot as plt

# Process counts and measured speedups
processes = [1, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40]
speedup =   [1.00, 3.54, 6.19, 8.04, 9.59, 9.80, 9.91, 10.38, 11.33, 11.36, 11.84]
ideal =     processes  # Ideal linear speedup

plt.figure(figsize=(6, 4))
plt.plot(processes, speedup, marker='o', label='Observed Speedup')
plt.plot(processes, ideal, linestyle='--', label='Ideal Speedup')

plt.title('Speedup vs. Process Count')
plt.xlabel('Number of Processes')
plt.ylabel('Speedup')
plt.xticks(processes)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("speedup_plot.png", dpi=300)
plt.show()
