#!/usr/bin/env python2.7

"""Generate water layers.
 
Places water molecules in a hexagonal close packed (HCP) with random
orientations for a given volume density and number of layers. The volume
density is calculated according to
Kumar et al., Phys. Rev. E, 72, 5, 051503 (2005).
"""

import numpy as np
from argparse import ArgumentParser
from ase import Atoms, io
from ase.calculators.tip3p import rOH, angleHOH
from potentials import help_text, get_potential
from random import random

__author__ = "Johan G. McQuillan"
__email__ = "johan.mcquillan.13@ucl.ac.uk"

L = 6.0221E23   # Avogadro's number
M = 18          # Molar mass of H2O in g/mol
m = M / L       # Molecular mass of H2O in g

parser = ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
parser.add_argument('N', type=int,
                    help='Number of molecules')
parser.add_argument('rho', type=float,
                    help='Volume density in g/cm^3')
parser.add_argument('c', type=float,
                    help='Super cell z dimension in Angstroms')
parser.add_argument('-L', type=int, default=1,
                    help='Number of layers')
group.add_argument('-w', type=float,
                   help='Effective confinement width in Angstroms')
group.add_argument('-V', type=str, nargs=2,
                   help=help_text())
args = parser.parse_args()

N = args.N
layers = args.L
n = N / layers
c = args.c

# Get effective width of confinement
if args.w is not None:
    w = args.w
else:
    # Get effective confinement width from the PotentialEnergySurface
    PES = get_potential(args.V[0])
    W = float(args.V[1])
    PES = get_potential(args.V[0])(W, c, au=False)
    w = PES.w_eff_su

# Calculate side length in Agstroms
a = np.sqrt(m * N / (w * args.rho))*1E12

# Get atomic coordinates of water molecule and create Atoms object
angle = angleHOH * np.pi / 180 / 2
pos =  [[0.,                0.,                 0.],
        [0., rOH*np.cos(angle),  rOH*np.sin(angle)],
        [0., rOH*np.cos(angle), -rOH*np.sin(angle)]]
h2o = Atoms('OH2', positions=pos)

# x and y length of cuboid around each H2O
h2o_xy = a / np.sqrt(n)
# z height of cuboid
h2o_z = w / layers
origin = np.array([h2o_xy, h2o_xy, h2o_z]) / 2

h2o.set_cell([h2o_xy, h2o_xy, h2o_z])
h2o.center()

# Initialise holder for all molecules
atoms = Atoms()

# Unit vectors
x = np.array([1, 0, 0])
y = np.array([0, 1, 0])
z = np.array([0, 0, 1])

# Iterate over all layers and grid points
# Fill up HCP lattice until N z_coords
count = 0
for l in range(layers):
    for i in range(int(np.ceil(np.sqrt(n)))):
        for j in range(int(np.ceil(np.sqrt(n)))):
            if count >= N:
                break
            
            # Randomise angles
            phi = 180.*random()
            theta = 180.*random()
            psi = 180.*random()
            
            # Arrange in HCP arrangement with minor random x and y displacement
            x = h2o_xy * (i - 0.01*(random() - 0.5)*0 + (1+(-1)**j+(-1)**l)/4.)
            y = h2o_xy * (j - 0.01*(random() - 0.5)*0 + (1+(-1)**l)/4.)
            z = c / 2. + ((1. - layers) / 2. + l) * h2o_z
            offset = np.array([x, y, z])
            
            h2o = Atoms('OH2', positions=pos+offset)
            Atoms.euler_rotate(h2o, phi, theta, psi,
                               center=h2o.get_center_of_mass())
            atoms += h2o
            count += 1

atoms.set_cell([a, a, c])
atoms.set_pbc(1)

io.write('test.xyz', atoms)
io.write('test.pdb', atoms)
