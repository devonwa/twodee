#!/usr/bin/env python
###############################################
# ase_db_get.py
# Template for selecting a list of atoms objects from an ASE DB

db_path = "{{ db_path }}"
select = "{{ select }}"

from ase.db import connect
db = connect(db_path)

images = []
for d in db.select(select):
    atoms = db.get_atoms(d.id)
    del atoms.constraints
    images += [atoms]
