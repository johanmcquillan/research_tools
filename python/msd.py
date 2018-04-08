#!/usr/bin/env python

"""Calculate the two-dimensional mean-square displacement of trajectories."""

import matplotlib.pyplot as plt
import os
import sys
from argparse import ArgumentParser
from packages.trajectories import read_trajectory, get_msd

parser = ArgumentParser()
parser.add_argument('traj', type=str, nargs='+', help='Trajectory files')
parser.add_argument('--start', type=int, default=1, help='Starting frame')
parser.add_argument('--end', type=int, default=-1, help='End frame (negative for no end)')
parser.add_argument('--stride', type=int, default=1, help='Frame stride length')
parser.add_argument('-S', type=str, default='O', nargs='+', help='Atomic species')
parser.add_argument('--delta-t', type=float, default=0.2, help='Time step in femtoseconds')
args = parser.parse_args()

__author__ = "Johan G. McQuillan"
__email__ = "johan.mcquillan.13@ucl.ac.uk"

dt = 5*1E-4     # ps
d_regime = 100  # ps

D = {}

fig, ax = plt.subplots(1)

for traj in args.traj:
    positions, steps = read_trajectory(traj, start=args.start, end=args.end, stride=args.stride,
                                       axes=(0, 1), species=args.S)
    time = steps*args.delta_t*1E-3

    msd, a, a_error, b, b_error = get_msd(positions, time, d_regime)

    label = os.path.splitext(os.path.basename(traj))[0]
    ax.plot(time, msd, label=label)
    # ax.plot(steps[intervals]*dt, msd_fft(z_coords[:, 0])*2)
    # ax.errorbar(steps[intervals]*dt, msd[], label=r'${}$K'.format(pdb[-8:-5]))
    ax.plot(time, a+b*time, alpha=0.3, color='k')
    
    D[label] = [b/2, b_error/2]

temp = sorted(D.keys())
diff = [D[T][0] for T in temp]
error = [D[T][1] for T in temp]

for T in temp:
    sys.stdout.write('{}:  {:.2e} +/- {:.2e}\n'.format(T, D[T][0], D[T][1]))
    sys.stdout.flush()

# Plot rdf
sys.stdout.write('Plotting 2D MSD\n')
sys.stdout.flush()

ax.set_xlabel(r'$t$ / $ps$')
ax.set_ylabel(r'$\langle r^2 (t) \rangle$ / $\AA^2$')
plt.legend()

plt.savefig('msd.png')
plt.show()

fig, ax = plt.subplots(1)
ax.scatter(temp, diff, marker='x', color='k')
ax.set_ylabel(r'$D$ / $\AA^2 ps^{-1}$')
ax.set_xlabel(r'$T$ / $K$')

plt.savefig('diff.png')
plt.show()

sys.stdout.write('Done.\n')
sys.stdout.flush()

