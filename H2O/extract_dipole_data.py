import re
import numpy as np

def extract_dipoles(file_path):
    steps = []
    electronic_dipoles = []
    ionic_dipoles = []
    step_count = 0
    electronic_dip = [0.0, 0.0, 0.0]
    ionic_dip = [0.0, 0.0, 0.0]
    in_electronic = False
    in_ionic = False
    with open(file_path, 'r') as file:
        for line in file:
            if 'iteration #' in line:
                step_count += 1
                steps.append(step_count)
            elif 'Electronic Dipole on Cartesian axes' in line:
                in_electronic = True
                electronic_dip = [0.0, 0.0, 0.0]
                continue
            elif 'Ionic Dipole on Cartesian axes' in line:
                in_electronic = False
                in_ionic = True
                ionic_dip = [0.0, 0.0, 0.0]
                continue
            elif in_electronic:
                match = re.search(r'^\s*(\d+)\s+([-]?\d*\.\d*(?:E[-+]\d+)?)', line)
                if match:
                    axis = int(match.group(1)) - 1
                    electronic_dip[axis] = float(match.group(2))
                    if axis == 2:  # Last axis (z)
                        electronic_dipoles.append(electronic_dip.copy())
            elif in_ionic:
                match = re.search(r'^\s*(\d+)\s+([-]?\d*\.\d*(?:E[-+]\d+)?)', line)
                if match:
                    axis = int(match.group(1)) - 1
                    ionic_dip[axis] = float(match.group(2))
                    if axis == 2:  # Last axis (z)
                        ionic_dipoles.append(ionic_dip.copy())
                        in_ionic = False
    if not steps or len(steps) != len(electronic_dipoles) or len(steps) != len(ionic_dipoles):
        raise ValueError("Mismatch in number of steps or dipole data")
    return steps, electronic_dipoles, ionic_dipoles

def calculate_total_dipole(electronic_dipoles, ionic_dipoles):
    total_dipoles_au = []
    total_dipoles_debye = []
    for e_dip, i_dip in zip(electronic_dipoles, ionic_dipoles):
        total_dip = np.array(e_dip) + np.array(i_dip)
        magnitude_au = np.sqrt(np.sum(total_dip**2))
        magnitude_debye = magnitude_au * 2.54175  # Convert Ry a.u. to Debye
        total_dipoles_au.append(magnitude_au)
        total_dipoles_debye.append(magnitude_debye)
    return total_dipoles_au, total_dipoles_debye

def save_dipole_data(file_path, steps, electronic_dipoles, ionic_dipoles, total_dipoles_au, total_dipoles_debye):
    with open(file_path, 'w') as f:
        f.write("# Step | Electronic Dipole (x,y,z) | Ionic Dipole (x,y,z) | Total Dipole (a.u.) | Total Dipole (Debye)\n")
        for step, e_dip, i_dip, t_au, t_debye in zip(steps, electronic_dipoles, ionic_dipoles, total_dipoles_au, total_dipoles_debye):
            f.write(f"{step} | {e_dip[0]:.6f},{e_dip[1]:.6f},{e_dip[2]:.6f} | {i_dip[0]:.6f},{i_dip[1]:.6f},{i_dip[2]:.6f} | {t_au:.6f} | {t_debye:.6f}\n")

def main():
    input_file = 'h2o_dipole.out'
    output_file = 'dipole_data.txt'
    try:
        steps, electronic_dipoles, ionic_dipoles = extract_dipoles(input_file)
        if not steps:
            print("Error: No dipole data extracted. Check file format.")
            return
        total_dipoles_au, total_dipoles_debye = calculate_total_dipole(electronic_dipoles, ionic_dipoles)
        save_dipole_data(output_file, steps, electronic_dipoles, ionic_dipoles, total_dipoles_au, total_dipoles_debye)
        print(f"Dipole data saved to '{output_file}'")
    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found in the H2O directory")
    except ValueError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()
