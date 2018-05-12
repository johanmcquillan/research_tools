#!/usr/bin/env python2.7

"""Plot the z-density profile averaged over a trajectory."""

import numpy as np
import matplotlib.pyplot as plt
import os
from argparse import ArgumentParser, RawTextHelpFormatter
from potentials import help_text, get_potential
from packages.trajectories import read_trajectory

__author__ = "Johan G. McQuillan"
__email__ = "johan.mcquillan.13@ucl.ac.uk"

L = 6.022E23    # Avogadro's number
M = 18.         # Molar mass of H2O
m = M / L       # Molecular mass of H2O

parser = ArgumentParser(formatter_class=RawTextHelpFormatter)
parser.add_argument('traj', type=str, nargs='+',
                    help='Trajectory files')
parser.add_argument('--start', type=int, default=0,
                    help='Starting frame')
parser.add_argument('--end', type=int, default=-1,
                    help='End frame (negative for no stop)')
parser.add_argument('--stride', type=int, default=1,
                    help='Frame stride length')
parser.add_argument('-S', type=str, default='O', nargs='+',
                    help='Atomic species')
parser.add_argument('-V', type=str, nargs=2, default=None,
                    help=help_text())
parser.add_argument('-b', type=int, default=100,
                    help='Number of bins for histogram')
args = parser.parse_args()

fig, axes = plt.subplots(1)
ax1 = axes

# Get PotentialEnergySurface
if args.V is not None:
    PES = get_potential(args.V[0])
    W = float(args.V[1])
    c = 100.
    V = PES(W, c, au=False)
    
    # Plot effective confinement width and potential shape
    w_eff = V.w_eff_su
    x = np.linspace((c - w_eff)/2., (c + w_eff)/2., 100)
    v = V.potential_su(x)
    ax1.plot(x-c/2, v)
    ax1.axvline(x=-w_eff/2, color='k', alpha=0.2)
    ax1.axvline(x=+w_eff/2, color='k', alpha=0.2)
    ax1.set_xlim([-W/2., +W/2.])

for traj in args.traj:
    z_coords, frames = read_trajectory(traj, start=args.start, end=args.end, stride=args.stride,
                                       axes=2, species=args.S)
    N = z_coords.shape[1]   # Number of atoms
    
    hist, bins = np.histogram(z_coords, bins=args.b, normed=True)
    avg_z = np.mean(z_coords)   # Used to centre plot around 0
    
    # Calculate mass density
    rho_z = m * hist / ((bins[-1] - bins[0])/args.b)
    
    label = os.path.splitext(os.path.basename(traj))[0]
    ax1.plot(bins[1:]-avg_z, hist, label=label)

if len(args.traj):
    ax1.legend()
ax1.set_xlabel(r'$z$ $[\AA]$')

plt.show()
