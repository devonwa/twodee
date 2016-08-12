#!/usr/bin/env python
###############################################
# ase_db_get.py
# Template for getting a generator from an ASE database

db_path = "{{ db_path }}"
select = {{ select }}

from ase.db import connect
db = connect(db_path)

db_rows = db.select(select)
