#!/usr/bin/env python
###############################################
# amp_remove_params.py
# Template for removing training data from an Amp directory

wd = "{{ workdir }}"

import os
os.unlink(os.path.join(wd, 'train-log.txt'))
os.unlink(os.path.join(wd, 'trained-parameters.json'))
