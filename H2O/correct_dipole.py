import numpy as np

def read_dipole_data(file_path):
    steps = []
    electronic_dipoles = []
    ionic_dipoles = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('|')
            if len(parts) != 5:
                continue
            step = int(parts[0].strip())
            e_dip = [float(x) for x in parts[1].strip().split(',')]
            i_dip = [float(x) for x in parts[2].strip().split(',')]
            steps.append(step)
            electronic_dipoles.append(e_dip)
            ionic_dipoles.append(i_dip)
    return steps, electronic_dipoles, ionic_dipoles

def calculate_total_dipole(electronic_dipoles, ionic_dipoles):
    total_dipoles_au = []
    total_dipoles_debye = []
    for e_dip, i_dip in zip(electronic_dipoles, ionic_dipoles):
        total_dip = np.array(e_dip) + np.array(i_dip)
        magnitude_au = np.sqrt(np.sum(total_dip**2))
        magnitude_debye = magnitude_au * 2.54175
        total_dipoles_au.append(magnitude_au)
        total_dipoles_debye.append(magnitude_debye)
    return total_dipoles_au, total_dipoles_debye

def save_corrected_data(file_path, steps, electronic_dipoles, ionic_dipoles, total_dipoles_au, total_dipoles_debye):
    with open(file_path, 'w') as f:
        f.write("# Step | Electronic Dipole (x,y,z) | Ionic Dipole (x,y,z) | Total Dipole (a.u.) | Total Dipole (Debye)\n")
        for step, e_dip, i_dip, t_au, t_debye in zip(steps, electronic_dipoles, ionic_dipoles, total_dipoles_au, total_dipoles_debye):
            f.write(f"{step} | {e_dip[0]:.6f},{e_dip[1]:.6f},{e_dip[2]:.6f} | {i_dip[0]:.6f},{i_dip[1]:.6f},{i_dip[2]:.6f} | {t_au:.6f} | {t_debye:.6f}\n")

def main():
    input_file = 'dipole_data.txt'
    output_file = 'dipole_data_corrected.txt'
    try:
        steps, electronic_dipoles, ionic_dipoles = read_dipole_data(input_file)
        if not steps:
            print("Error: No data read from file. Check 'dipole_data.txt'.")
            return
        total_dipoles_au, total_dipoles_debye = calculate_total_dipole(electronic_dipoles, ionic_dipoles)
        save_corrected_data(output_file, steps, electronic_dipoles, ionic_dipoles, total_dipoles_au, total_dipoles_debye)
        print(f"Corrected dipole data saved to '{output_file}'")
    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found in the H2O directory")
    except ValueError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()
