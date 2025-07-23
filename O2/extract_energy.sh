#!/bin/bash

echo "# BondLength(Å)    TotalEnergy(Ry)" > energy.dat

for bond in $(seq 1.0 0.1 2.0); do
    energy=$(grep "!    total energy" o2_${bond}.out | awk '{print $5}')
    printf "%4.1f           %s\n" "$bond" "$energy" >> energy.dat
done
