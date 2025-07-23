import os

# Directory containing QE output files
output_dir = os.path.expanduser("~/H2")
output_files = [f"H2_{bond_length}.out" for bond_length in [0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9]]

# File to save energy data
output_data_file = os.path.join(output_dir, "energy.dat")

# Initialize results
results = []

# Parse each output file
for filename in output_files:
    filepath = os.path.join(output_dir, filename)
    try:
        with open(filepath, "r") as file:
            for line in file:
                if "!    total energy" in line:
                    energy = float(line.split()[-2])  # Extract energy value
                    bond_length = filename.split("_")[1].replace(".out", "")
                    results.append(f"{bond_length} {energy}")
                    break
    except FileNotFoundError:
        print(f"Warning: {filename} not found.")
    except Exception as e:
        print(f"Error reading {filename}: {e}")

# Write results to energy.dat
if results:
    with open(output_data_file, "w") as file:
        file.write("\n".join(results))
    print(f"Data written to {output_data_file}")
else:
    print("No data to write. Check your QE output files.")
