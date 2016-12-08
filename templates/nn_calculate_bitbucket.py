#!/usr/bin/env python
from amp import Amp
from ase.db import connect
import os
import shutil

import twodee as td

# Required arguments from template
iteration = {{ iteration }}
framework = {{ framework }}
work_dir = "{{ work_dir }}"
load_path = "{{ load_path }}"
db_path = "{{ db_path }}"
selection = {{ selection }}
{% if overwrite %}
overwrite = {{ overwrite }}
{% else %}
overwrite = False
{% endif %}
cores = {{ cores }}

# Set working directory
os.chdir(work_dir)

# Create DB update strings
update_str = "_".join([str(f) for f in framework])
update_key = "nn_" + update_str + "_iter{}".format(iteration)

# Load Amp object
calc = Amp.load(load_path, cores=cores)

# Get atoms from the database
images = []
new_db = os.path.join(os.path.dirname(load_path), "DB.db")
shutil.copyfile(db_path, new_db)
db = connect(new_db)
for d in db.select(selection):
    if not overwrite and update_key in d.key_value_pairs.keys():
        continue

    atoms = db.get_atoms(d.id)
    atoms.set_calculator(calc)
    energy = atoms.get_potential_energy()

    update = {update_key: energy}
    db.update(d.id, **update)
