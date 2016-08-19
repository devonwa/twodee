#!/usr/bin/env python

import os
import twodee as td


work_dir = "{{ work_dir }}"
json_path = "{{ json_path }}"
size = {{ size }}
pore_size = {{ pore_size }}

# Set working directory
os.chdir(work_dir)

cans = td.candidates(mat='graphene', layers=1, size=size, pores=pore_size,
                     json_path=json_path, silent=False,
                     write=True, overwrite_file=False)
