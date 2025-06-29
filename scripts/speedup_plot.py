import matplotlib.pyplot as plt
import os

# Process counts and measured speedups
processes = [1, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40]
speedup =   [1.00, 3.54, 6.19, 8.04, 9.59, 9.80, 9.91, 10.38, 11.33, 11.36, 11.84]
ideal =     processes  # Ideal linear speedup

# Calculate efficiency
efficiency = [s/p * 100 for s, p in zip(speedup, processes)]

# Create output directory if it doesn't exist
os.makedirs("../results", exist_ok=True)

# Create speedup plot
plt.figure(figsize=(10, 6))

# Subplot 1: Speedup
plt.subplot(1, 2, 1)
plt.plot(processes, speedup, marker='o', label='Observed Speedup', linewidth=2, markersize=6)
plt.plot(processes, ideal, linestyle='--', label='Ideal Speedup', linewidth=2)
plt.title('Speedup vs. Process Count')
plt.xlabel('Number of Processes')
plt.ylabel('Speedup')
plt.xticks(processes[::2])  # Show every other tick for readability
plt.grid(True, alpha=0.3)
plt.legend()

# Subplot 2: Efficiency
plt.subplot(1, 2, 2)
plt.plot(processes, efficiency, marker='s', color='red', linewidth=2, markersize=6)
plt.title('Parallel Efficiency vs. Process Count')
plt.xlabel('Number of Processes')
plt.ylabel('Efficiency (%)')
plt.xticks(processes[::2])
plt.grid(True, alpha=0.3)
plt.ylim(0, 105)

plt.tight_layout()
plt.savefig("../results/speedup_analysis.png", dpi=300, bbox_inches='tight')
print(f"Speedup analysis saved to ../results/speedup_analysis.png")

# Print performance summary
print("\nPerformance Summary:")
print("===================")
print("Processes | Speedup | Efficiency")
print("----------|---------|----------")
for p, s, e in zip(processes, speedup, efficiency):
    print(f"{p:8d} | {s:7.2f} | {e:7.1f}%")

plt.show()
