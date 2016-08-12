#!/usr/bin/env python
###############################################
# calc.py
# Template for calculating energies using an ASE compliant calculator
# Required vars: calc (ase.calculators)

if 'energies' not in locals():
    energies = {}

for d in db_rows:
    atoms = db.get_atoms(d.id)
    atoms.set_calculator(calc)
    energy = atoms.get_potential_energy()
    energies[d.id: energy]
