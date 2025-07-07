# DFT Projects

This repository contains Density Functional Theory (DFT) simulations performed using Quantum ESPRESSO (QE). The `Si` subfolder includes simulations for silicon.

## Si Subfolder
Simulations for silicon in the diamond cubic structure.

### Files
- **Input Files**:
  - `si_scf.in`: SCF calculation input.
  - `si_nscf.in`: NSCF calculation input.
  - `si_dos.in`: DOS calculation input.
  - `si_vc_relax.in`: Variable-cell relaxation input.
- **Output Files**:
  - `si_scf.out`: SCF output.
  - `si_nscf.out`: NSCF output (Fermi energy: 6.2851 eV, bandgap: 0.7527 eV).
  - `si_dos.out`: DOS output.
  - `si_vc_relax.out`: Variable-cell relaxation output.
- **Data Files**:
  - `si.dos`: DOS data.
- **Scripts**:
  - `plot_dos.py`: Python script to plot DOS.
  - `Task`: Job script or documentation.
- **Excluded Files**:
  - `out/`: Large files (e.g., `silicon.save`) not uploaded.
  - `Si.pbe-n-kjpaw_psl.1.0.0.UPF`: Pseudopotential (download from [QE website](https://www.quantum-espresso.org/pseudopotentials)).

### Requirements
- Quantum ESPRESSO v7.4.1
- Python with Matplotlib

### How to Reproduce
1. Download `Si.pbe-n-kjpaw_psl.1.0.0.UPF` from QE website.
2. Run: `pw.x < si_scf.in > si_scf.out`
3. Run: `pw.x < si_nscf.in > si_nscf.out`
4. Run: `dos.x < si_dos.in > si_dos.out`
5. Plot: `python plot_dos.py`

### Notes
- DOS bandgap (~0.260 eV) is underestimated due to smearing (0.0200 Ry). NSCF bandgap is 0.7527 eV.
