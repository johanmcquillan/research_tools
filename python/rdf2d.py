#!/usr/bin/env python

"""Calculate the two-dimensional radial distribution function.

Given an i-PI trajectory in PDB or XYZ format, plot the 2D RDF. Note that this is the
partial RDF for the same species (g_{OO}(r) or g_{HH}(r)).
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from argparse import ArgumentParser, RawTextHelpFormatter
from packages.trajectories import read_trajectory, get_rdf

__author__ = "Johan G. McQuillan"
__email__ = "johan.mcquillan.13@ucl.ac.uk"

parser = ArgumentParser(formatter_class=RawTextHelpFormatter)
parser.add_argument('traj', type=str, nargs='+', help='Trajectory files')
parser.add_argument('--start', type=int, default=0, help='Starting frame')
parser.add_argument('--end', type=int, default=-1, help='End frame (negative for all frames)')
parser.add_argument('--stride', type=int, default=1, help='Frame stride')
parser.add_argument('-S', type=str, default='O', help='Atomic species')
parser.add_argument('-d', metavar='dr', type=float, default=0.1, help='Width of histogram bin in Angstroms')
parser.add_argument('-m', metavar='r_max', type=float, default=10., help='Maximum distance in Angstroms')
parser.add_argument('-v', action='store_true', help='Verbose')
args = parser.parse_args()

# Initialise plot and horizontal line at y=1
fig, ax = plt.subplots(1)
ax.axhline(1, color='k', alpha=0.5)

for traj in args.traj:
    positions, frames, cell = read_trajectory(traj, start=args.start, end=args.end, stride=args.stride, cell=True, species=args.S, throb=args.v)
    N = positions.shape[1]
    positions = np.swapaxes(positions, 1, 2)
    
    if args.v:
        sys.stdout.write('{} Particles over {} Frames\n'.format(positions.shape[2], positions.shape[0]))
        sys.stdout.write('Calculating Pair Distances...\n')
        sys.stdout.flush()
    
    r, g, i_max1, i_max2, coord = get_rdf(positions, cell, dr=args.d, r_max=args.m, prog_bar=args.v, verbose=args.v)
    label = os.path.splitext(os.path.basename(traj))[0]
    ax.plot(r, g, label=label)

ax.set_ylabel('$g_{OO}^{xy}(r)$')
ax.set_xlabel(r'$r$ [$\AA$]')

if len(args.traj) > 1:
    plt.legend()
plt.show()

if args.v:
    sys.stdout.write('Done.\n')
    sys.stdout.flush()
