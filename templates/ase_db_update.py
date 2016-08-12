#!/usr/bin/env python
###############################################
# ase_db_update.py
# Template for updating the rows of an ASE database
# Required vars: db (ase.db)
#                db_rows (Generator[ase.db row instance])

key = {{ key }}

for d in zip(db_rows, energies):
    db.update(d.id, key)
