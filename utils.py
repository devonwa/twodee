import os
import re
import shutil
import subprocess
import sys
import time
from types import ModuleType

from ase import Atoms
from ase.io import write as ase_write
from ase.visualize import view
import numpy as np
from vasp import Vasp
from vasp.vasprc import VASPRC


def bp(info=None):
    """A breakpoint to view something and stop the rest of the script."""
    if isinstance(info, Atoms):
        view(info)    
    elif isinstance(info, list) and info:
        if all(isinstance(i, Atoms) for i in info):
            view(info)
        else:
            print(info)
    elif info is not None:
        print(info)

    sys.exit()


def print_code(code):
    """Print an orgmode code block containing the given code."""
    print("#+BEGIN_SRC python :results output org drawer")
    print(code)
    print("#+END_SRC")



def broken_calcs(dft_path, delete=False, silent=False):
    """Return a list of broken calculation paths."""

    VASPRC['mode'] = None

    # Don't print to stdout for awhile. Vasp() prints error traces.
    stdout = sys.stdout
    f = open(os.devnull, 'w')
    sys.stdout = f

    broken_dirs = []
    for path in calc_paths(dft_path):
        try:
            Vasp(path)
        except:
            broken_dirs.append(path)

    sys.stdout = stdout
    if not silent:
        if broken_dirs:
            print("Broken directories:")
            for b in broken_dirs:
                print("    {}".format(b))
        else:
            print("No broken directories found.")

    # Recursively delete the broken directories
    if delete:
        for b in broken_dirs:
            shutil.rmtree(b)
    
    return broken_dirs


def calc_paths(path):
    """Walk a path and return directories with Vasp() calculations."""
    calc_paths = []
    for pwd in os.walk(path):
        path = pwd[0]
        if calc_output_files(path):
            calc_paths.append(pwd[0])

    return calc_paths


def calc_output_files(path):
    """Return a list of the output files in the directory."""
    files = os.listdir(path)
    regex = '[0-9]*.gilgamesh.cheme.cmu.edu.OU'
    output_files = [f for f in files if re.search(regex, f)]
    return output_files


def center(atoms):
    """Return the position (x,y,z) of the center of the cell."""
    cell = np.array(atoms.get_cell())
    center = (cell[0] + cell[1]) / 2
    center += cell[2] / 2
    return center


def is_the_same(x, fun, *args):
    """True if the object is unchanged during the function call."""
    import copy
    y = copy.deepcopy(x)
    fun(*args)
    return x == y


def paint_atoms(atoms, indices, sym=None, layers=None):
    """Update the chemical symbol of atoms in the list of indices."""
    if sym is not None:
        symbols = sym
    else:
        symbols = ["N", "O", "B", "F"]

    if layers is not None:
        for i in indices:
            for j, l in enumerate(layers):
                if i in l:
                    atoms[i].symbol = symbols[j % len(symbols)]
    else:
        for i in indices:
            atoms[i].symbol = symbols[0]


def qsub(file_path, walltime='168:00:00', mem="2GB", nodes=1, ppn=1):
    """Submit a script to the Torque queue on Gilgamesh."""
    wd = os.path.dirname(os.path.realpath(file_path))
    os.chmod(file_path, 0777)

    cmd = '''#!/bin/bash
#PBS -N {0}
#PBS -o {0}
#PBS -e {0}
#PBS -l nodes={2}:ppn={3}
#PBS -l walltime={4}
#PBS -l mem={5}
#PBS -joe
cd $PBS_O_WORKDIR
{1}
#end'''.format(wd, file_path, nodes, ppn, walltime, mem)

    submit_path = os.path.join(wd, 'submit.sh')
    with open(submit_path, 'w') as f:
        f.write(cmd)

    subprocess.call(['qsub', submit_path])
    time.sleep(1)
    os.remove(submit_path)
    print("Job submitted.")


def set_vacuum(atoms, vacuum):
    """Center atoms in the z-direction in a cell of size vacuum.

    Centers atoms in a unitcell with space above and below of 1/2 * vacuum. Assumes the current unitcell is centered and cell length changes only in the z-direction.
    
    Args:
        atoms (Atoms): Unitcell of atoms
        vacuum (float): Height of new unitcell

    Returns:
        An Atoms object with the new cell height.
    """
    cell = atoms.get_cell()
    center_old = cell[2][2] / 2.
    center_new = vacuum / 2.
    cell[2][2] = vacuum
    atoms.set_cell(cell)

    for atom in atoms:
        atom.position[2] = center_new - (center_old - atom.position[2])


def result(name, calc, fu=None, per_atom=False):
    """Print a brief calculation report."""
    atoms = calc.get_atoms()
    energy = atoms.get_potential_energy()

    if energy is None:
        stat = "Inprogress."
    else:
        time = calc.get_elapsed_time()
        if per_atom:
            stat = "Energy/atom = {:0.3f}. Calc time: {:.0f} min.".format(energy / len(atoms), time/60.)
        else:
            stat = "Energy = {:0.3f}. Calc time: {:.0f} min.".format(energy, time/60.)

    print(name + ": " + stat)


def status_converged(energy, time):
    print("Final structure calculation: Energy/f.u. = {:0.3f}. Calculation time: {:.0f} min.".format(energy, time/60.))


def status_inprogress():
    print("Final structure calculation: In progress.")


def status_unconverged(i):
    print("Distance: {:5.2f}. Did not converge.".format(i))


def spline(x, y, points=200):
    """Return x and y spline values over the same range as x."""
    from scipy.interpolate import interp1d
    spline = interp1d(x, y, kind='cubic')
    x_lin = np.linspace(x[0], x[-1], points)
    y_interp = spline(x_lin)

    return [x_lin, y_interp]


def print_image(path, data, fig_name=None, caption=None):
    if caption is not None:
        print('#+CAPTION: {}'.format(caption))
    if fig_name is not None:
        print('#+NAME: fig:{}'.format(fig_name))
    print(write_image(path, data))


def write_image(path, data, options=None):
    #file_path = './img/' + path
    file_path = path
    directory = file_path[:file_path.rfind('/')]
    if not os.path.exists(directory):
        os.makedirs(directory)

    if isinstance(data, ModuleType):
        if data.__name__ == "matplotlib.pyplot":
            data.savefig(file_path)
    elif isinstance(data, Atoms):
        #file_path += '.png'
        ase_write(file_path, data)

    else:
        print("No functionality for type = {}".format(type(data)))

    return '[[' + file_path + ']]'
