#!/usr/bin/env python
from amp import Amp
from ase.db import connect
import numpy as np
import os
import shutil

import twodee as td

# Required arguments from template
from ase.optimize import {{ optimizer }}
work_dir = "{{ work_dir }}"
db_path = "{{ db_path }}"
nn_path = "{{ nn_path }}"
output_path = "{{ output_path }}"
selection = ["{{ selection }}"]
cores = {{ cores }}
fmax = {{ fmax }} # force goal
steps = {{ steps }} # max steps


# Set working directory
os.chdir(work_dir)


# Get atoms from the database
db = connect(db_path)
d = next(db.select(selection))
atoms = db.get_atoms(d.id)


# Load Amp object
converged = "trained-parameters.json"
checkpoints = "checkpoint-parameters.json"
nn_files = os.listdir(nn_path)

if converged in nn_files:
    params = converged
elif checkpoints in nn_files:
    params = checkpoints

load_path = os.path.join(nn_path, params)
calc = Amp(load=load_path, cores=cores)
atoms.set_calculator(calc)


# Run MD simulation
traj_path = os.path.join(output_path, "out.traj")
log_path = os.path.join(output_path, "log.txt")
dyn = {{ optimizer }}(atoms, trajectory=traj_path, logfile=log_path)
dyn.run(fmax=fmax, steps=steps)
