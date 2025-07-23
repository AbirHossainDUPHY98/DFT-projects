import matplotlib.pyplot as plt

data_file = '/home/abir-hossain/H2/energy.dat'
bond_lengths = []
energies = []

with open(data_file, 'r') as file:
    for line in file:
        values = line.split()
        bond_lengths.append(float(values[0]))
        energies.append(float(values[1]))

plt.figure(figsize=(8, 6))
plt.plot(bond_lengths, energies, marker='o', linestyle='-', color='b')
plt.xlabel('Bond Length (Å)')
plt.ylabel('Energy (Ry)')
plt.title('Potential Energy Curve')
plt.grid()
plt.savefig('potential_energy_curve.png')
plt.show()
