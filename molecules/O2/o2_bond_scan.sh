#!/bin/bash

# Number of MPI processes
NPROCS=6

for bond in $(seq 1.0 0.1 2.0); do
    alat=10.0
    bond_au=$(echo "$bond * 1.889725989" | bc -l)
    frac=$(echo "scale=6; $bond_au / $alat / 2.0" | bc -l)
    
    # Format fraction to 6 decimals for neatness
    frac_plus=$(printf "%.6f" $(echo "0.5 + $frac" | bc -l))
    frac_minus=$(printf "%.6f" $(echo "0.5 - $frac" | bc -l))
    
    cat > o2_${bond}.in << EOF
&CONTROL
  calculation = 'scf',
  prefix = 'o2',
  outdir = './tmp_${bond}',
  pseudo_dir = './',
  verbosity = 'high',
/

&SYSTEM
  ibrav = 1,
  celldm(1) = ${alat},
  nat = 2,
  ntyp = 1,
  ecutwfc = 60.0,
  ecutrho = 240.0,
  occupations = 'smearing',
  smearing = 'gaussian',
  degauss = 0.01,
  nspin = 2,
  starting_magnetization(1) = 1.0,
/

&ELECTRONS
  conv_thr = 1.0d-8,
  mixing_beta = 0.7,
/

ATOMIC_SPECIES
  O  15.999  O.pbe-n-kjpaw_psl.1.0.0.UPF

ATOMIC_POSITIONS crystal
  O  0.0  0.0  ${frac_plus}
  O  0.0  0.0  ${frac_minus}

K_POINTS automatic
  1 1 1 0 0 0
EOF

    echo "Running for bond length ${bond} Å..."
    mpirun -np $NPROCS pw.x -i o2_${bond}.in > o2_${bond}.out

done
