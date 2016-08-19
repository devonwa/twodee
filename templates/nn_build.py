#!/usr/bin/env python
from amp import Amp
from amp.descriptor.gaussian import Gaussian
from amp.model import LossFunction
from amp.model.neuralnetwork import NeuralNetwork
from amp.utilities import Annealer
from ase.db import connect
import os
import twodee as td

# Required arguments from template
iteration = {{ iteration }}
framework = {{ framework }}
work_dir = "{{ work_dir }}"
db_path = "{{ db_path }}"
label = "{{ label }}"
dblabel = "{{ dblabel }}"
cutoff = {{ cutoff }}
energy_rmse = {{ energy_rmse }}
force_rmse = {{ force_rmse }}
cores = {{ cores }}

# Set working directory
os.chdir(work_dir)

# Get atoms from the database
images = []
db = connect(db_path)
for d in db.select(['iteration<={}'.format(iteration, 'train_set=True')]):
    atoms = db.get_atoms(d.id)
    del atoms.constraints
    images.append(atoms)

# Build Amp object
framework_str = "-".join([str(f) for f in framework])
label = label
dblabel = dblabel
desc = Gaussian(cutoff=cutoff)
model = NeuralNetwork(hiddenlayers=framework)
calc = Amp(label=label,
           dblabel=dblabel,
           descriptor=desc,
           model=model,
	   cores=cores)
loss = LossFunction(convergence={'energy_rmse': energy_rmse,
                                 'force_rmse': force_rmse})
calc.model.lossfunction = loss
           
# Perform simulated annealing for global search
Annealer(calc=calc, images=images)

# Train the network
calc.train(images=images)
