#!/usr/bin/env python2.7

"""Convert an MB-pol water PDB trajectory into an XYZ file."""

import numpy as np
import os
from argparse import ArgumentParser
from packages.ipi_builder import mbpol_config, mbpol_control, mbpol_field
from packages.trajectories import read_pdb, write_xyz

parser = ArgumentParser()
parser.add_argument('pdb', type=str,
                    help='PDB file')
parser.add_argument('-z', action='store_true',
                    help='Offset z values of first frame such that the average z of O is half of cell height')
parser.add_argument('--dlpoly', action='store_true',
                    help='Generate DLPOLY2 FIELD and CONFIG.01 files')
args = parser.parse_args()

in_file_name = os.path.splitext(args.pdb)[0]

coords, frames, cells, species = read_pdb(args.pdb, species=('O', 'H', 'H1', 'H2'), cell=True,
                                          spec_list=True, axes=(0, 1, 2))

h = cells[0]
c = cells[0, 2]
# Get indices of oxygen atoms
o_mask = np.array([i for i, s in enumerate(species) if s == 'O'])

N = coords.shape[1]
assert N % 3 == 0
n = N / 3

if args.z:
    z_offset = c/2. - np.mean([coords[0, o_mask, 2] for coord in coords])
    coords[:, :, 2] += z_offset

write_xyz(in_file_name+'.xyz', coords, cells, species, frames=frames)

if args.dlpoly:
    with open('FIELD', 'w') as field_file:
        field_file.write(mbpol_field(n))
    with open('CONFIG.01', 'w') as config_file:
        config_file.write(mbpol_config(n, h))

# if args.control:
#     control = 'CONTROL'
#     with open(control, 'w') as control_file:
#         control_file.write(mbpol_control())
#
