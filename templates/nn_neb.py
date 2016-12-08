#!/usr/bin/env python
from amp import Amp
from ase.db import connect
from ase.io.trajectory import Trajectory
from ase.neb import NEB
import numpy as np
import os
import shutil

import twodee as td

# Required arguments from template
from ase.optimize import {{ optimizer }}
work_dir = "{{ work_dir }}"
nn_path = "{{ nn_path }}"
output_path = "{{ output_path }}"
cores = {{ cores }}
fmax = {{ fmax }} # force goal
steps = {{ steps }} # max steps


# Set working directory
os.chdir(work_dir)


# Get images
traj = Trajectory(os.path.join(output_path, "in.traj"))
images = [i for i in traj]
neb = NEB(images)
neb.interpolate()


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

for image in images[1:-1]:
    image.set_calculator(calc)


# Optimize
traj_path = os.path.join(output_path, "out.traj")
log_path = os.path.join(output_path, "log.txt")
optimizer = {{ optimizer }}(neb, trajectory=traj_path, logfile=log_path)
optimizer.run(fmax=fmax, steps=steps)
