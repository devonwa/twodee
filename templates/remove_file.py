#!/usr/bin/env python
###############################################
# remove_file.py
# Template for removing a file
# TODO devon: Deprecate? Consider a more generic setup.

file_path = "{{ file_path }}"

import os
os.unlink(file_path)
