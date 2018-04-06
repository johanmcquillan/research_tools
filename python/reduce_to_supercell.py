#!/usr/bin/env python

"""Use periodic boundary conditions to reduce atomic coordinates to a single supercell"""

import os
from argparse import ArgumentParser
from packages.trajectories import read_trajectory, write_trajectory, reduce_to_supercell

__author__ = "Johan G. McQuillan"
__email__ = "johan.mcquillan.13@ucl.ac.uk"

parser = ArgumentParser()
parser.add_argument('traj', type=str, help='Trajectory files')
parser.add_argument('suffix', type=str, nargs='?', default='reduced', help='Suffix for new file')
args = parser.parse_args()

# Generate new file name
file_name, ext = os.path.splitext(os.path.basename(args.traj))
new_name = file_name + '_' + args.suffix + ext

atoms, frames, cell, species_list = read_trajectory(args.traj, start=0, end=1, cell=True,
                                                    spec_list=True, species=('O', 'H'),
                                                    axes=(0, 1, 2))
atoms = reduce_to_supercell(atoms, cell, molecules=True)
write_trajectory(new_name, atoms, cell, species_list)
