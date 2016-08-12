#!/usr/bin/env python
###############################################
# amp_train.py
# Template for training an Amp object
# Required vars: db_rows (Generator[ase.db row instance])
#                calc (Amp)
#                images (Atoms)

train_args = {{ train_args }}

images = []
for d in db_rows:
    atoms = db.get_atoms(d.id)
    del atoms.constraints
    images += [atoms]

{% if gs_type %}
# Define global search object
from amp import {{ gs_type }}
gs = {{ gs_type }}(**{{ gs_args }})
train_args['global_search'] = gs
{% endif %}

# Run training calculations
calc.train(images=images,
           **train_args)
