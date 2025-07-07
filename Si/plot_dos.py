import numpy as np
import matplotlib.pyplot as plt

# Load the DOS data
# Format: energy (eV), total DOS, integrated DOS (3 columns)
data = np.loadtxt('si.dos')

energy = data[:, 0]
dos = data[:, 1]

plt.figure(figsize=(8,6))
plt.plot(energy, dos, color='blue')
plt.axvline(x=0, color='red', linestyle='--', label='Fermi Level')
plt.xlabel('Energy (eV)')
plt.ylabel('Density of States (states/eV)')
plt.title('Density of States for Silicon')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('si_dos.png', dpi=300)
plt.show()
