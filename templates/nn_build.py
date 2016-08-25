#!/usr/bin/env python
from amp import Amp
from amp import SimulatedAnnealing
from amp.descriptor import Gaussian
from amp.regression import NeuralNetwork
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
for d in db.select(['iteration<={}'.format(iteration), 'train_set=True']):
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
           regression=model)
           
# Train the network
calc.train(images=images,
           data_format='db',
           cores=cores,
           energy_goal=energy_rmse,
           force_goal=force_rmse,
	   global_search=SimulatedAnnealing(temperature=70,
					    steps=50),
           extend_variables=False)
