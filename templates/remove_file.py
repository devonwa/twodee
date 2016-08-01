"""Remove a file"""
file_path = "{{ file_path }}"

import os
os.unlink(file_path)
