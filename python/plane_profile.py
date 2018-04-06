#!/usr/bin/env python

"""Calculate the distribution of in-plane and out-of-plane water molecules."""

import matplotlib.pyplot as plt
import os
from argparse import ArgumentParser, RawTextHelpFormatter
from packages.trajectories import read_trajectory, plane_distribution

__author__ = "Johan G. McQuillan"
__email__ = "johan.mcquillan.13@ucl.ac.uk"

parser = ArgumentParser(formatter_class=RawTextHelpFormatter)
parser.add_argument('traj', nargs='+', type=str, help='Input trajectory (.pdb or .xyz)')
parser.add_argument('--start', type=int, default=0, help='Starting frame')
parser.add_argument('--end', type=int, default=-1, help='End frame (negative for all frames)')
parser.add_argument('--stride', type=int, default=1, help='Frame stride')
args = parser.parse_args()

# Initialise plot
fig, ax = plt.subplots(1)

for traj in args.traj:
    hydrogens, frames = read_trajectory(traj, start=args.start, end=args.end, stride=args.stride,
                                        axes=(0, 1, 2), species=('H',))
    hist, bins = plane_distribution(hydrogens)
    label = os.path.splitext(os.path.basename(traj))[0]
    ax.plot(bins[1:], hist, label=label)

ax.set_xlabel(r'$\theta$ $[\deg]$')
if len(args.traj) > 1:
    ax.legend()
plt.show()
