import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt

def read_dipole_data(file_path):
    steps = []
    total_dipoles_au = []
    total_dipoles_debye = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('|')
            if len(parts) != 5:
                continue
            step = int(parts[0].strip())
            t_au = float(parts[3].strip())
            t_debye = float(parts[4].strip())
            steps.append(step)
            total_dipoles_au.append(t_au)
            total_dipoles_debye.append(t_debye)
    return steps, total_dipoles_au, total_dipoles_debye

def plot_dipole_trend(steps, total_dipoles_au, total_dipoles_debye):
    plt.figure(figsize=(10, 6))
    plt.plot(steps, total_dipoles_au, label='Total Dipole Moment (a.u.)', marker='o', color='#1f77b4')
    plt.plot(steps, total_dipoles_debye, label='Total Dipole Moment (Debye)', marker='s', color='#ff7f0e')
    plt.xlabel('SCF Iteration')
    plt.ylabel('Dipole Moment')
    plt.title('Total Dipole Moment Trend Over SCF Iterations')
    plt.legend()
    plt.grid(True)
    plt.savefig('dipole_trend.png')
    plt.close()

def main():
    file_path = 'dipole_data_corrected.txt'
    try:
        steps, total_dipoles_au, total_dipoles_debye = read_dipole_data(file_path)
        if not steps:
            print("Error: No data read from file. Check 'dipole_data_corrected.txt'.")
            return
        print("Steps:", steps)
        print("Total Dipole (a.u.):", total_dipoles_au)
        print("Total Dipole (Debye):", total_dipoles_debye)
        plot_dipole_trend(steps, total_dipoles_au, total_dipoles_debye)
        print("Plot saved as 'dipole_trend.png' in the H2O directory")
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found in the H2O directory")
    except ValueError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()
