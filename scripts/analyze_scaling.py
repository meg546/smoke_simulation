import numpy as np
import matplotlib.pyplot as plt

# Load data
data = np.loadtxt('scaling_results.txt')
procs = data[:,0]
times = data[:,1]
speedup = data[:,2]
efficiency = data[:,3]
serial_fraction = data[:,4]

# Create subplots
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))

# Speedup plot
ax1.plot(procs, speedup, 'bo-', label='Actual')
ax1.plot(procs, procs, 'r--', label='Linear')
ax1.set_xlabel('Number of Processors (p)')
ax1.set_ylabel('Speedup')
ax1.legend()
ax1.grid(True)
ax1.set_title('Speedup vs. Processors')

# Efficiency plot
ax2.plot(procs, efficiency*100, 'go-')
ax2.set_xlabel('Number of Processors (p)')
ax2.set_ylabel('Efficiency (%)')
ax2.grid(True)
ax2.set_title('Parallel Efficiency')

# Execution time plot
ax3.plot(procs, times, 'mo-')
ax3.set_xlabel('Number of Processors (p)')
ax3.set_ylabel('Execution Time (s)')
ax3.grid(True)
ax3.set_title('Execution Time vs. Processors')

# Serial fraction plot
ax4.plot(procs, serial_fraction, 'ko-')
ax4.set_xlabel('Number of Processors (p)')
ax4.set_ylabel('Serial Fraction')
ax4.grid(True)
ax4.set_title('Estimated Serial Fraction')

plt.tight_layout()
plt.savefig('scaling_analysis.png')
plt.show() 