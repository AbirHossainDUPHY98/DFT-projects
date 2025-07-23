#!/bin/bash

# Create input files for bond lengths 0.5 to 1.9
for R_BOND in $(seq 0.5 0.2 1.9); do
    cat > H2_${R_BOND}.in << EOF
&CONTROL
    calculation = 'scf',
    outdir = './',
    prefix = 'h2',
    pseudo_dir = './',
    verbosity = 'high',
/

&SYSTEM
    ibrav = 1,
    celldm(1) = 10.0,
    nat = 2,
    ntyp = 1,
    ecutwfc = 40.0,
    ecutrho = 320.0,
/

&ELECTRONS
    conv_thr = 1.0d-8,
    mixing_beta = 0.7,
/

ATOMIC_SPECIES
    H  1.00794  H.pbe-kjpaw_psl.1.0.0.UPF

ATOMIC_POSITIONS {angstrom}
    H  0.0 0.0 0.0
    H  0.0 0.0 ${R_BOND}

K_POINTS {automatic}
    1 1 1 0 0 0
EOF
done
