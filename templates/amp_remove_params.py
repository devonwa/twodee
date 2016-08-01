import os
wd = "{{ workdir }}"
os.unlink(os.path.join(wd, 'train-log.txt'))
os.unlink(os.path.join(wd, 'trained-parameters.json'))
