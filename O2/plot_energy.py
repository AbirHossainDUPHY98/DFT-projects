import matplotlib.pyplot as plt
import numpy as np

# Load data
data = np.loadtxt('energy.dat', comments='#')
bond_lengths = data[:, 0]
energies = data[:, 1]

# Plot
plt.figure(figsize=(8, 6))
plt.plot(bond_lengths, energies, 'o-', color='blue')
plt.xlabel('Bond Length (Å)')
plt.ylabel('Total Energy (Ry)')
plt.title('O₂ Potential Energy Curve')
plt.grid(True)
plt.tight_layout()
plt.savefig('o2_energy_curve.png')
plt.show()
