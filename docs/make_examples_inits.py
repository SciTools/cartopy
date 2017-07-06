from glob import glob
import os
import shutil as sh
import sys

kind = sys.argv[1]

# First delete `__init__.py` in all example folders
path_examples = os.path.join('..', 'lib', 'cartopy', 'examples')
path_folders = os.path.join(path_examples, '*')
init_folders = glob(path_folders)
init_folders += [path_examples]

for folder in init_folders:
    if not os.path.isdir(folder):
        continue

    # If specified, remove all init files
    if kind == 'remove':
        if os.path.exists(os.path.join(folder, '__init__.py')):
            os.remove(os.path.join(folder, '__init__.py'))
    # Otherwise add the init files back
    elif kind == 'add':
        with open(os.path.join(folder, '__init__.py'), 'w'):
            pass
