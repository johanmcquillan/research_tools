#!/usr/bin/env python

"""Calculate the two-dimensional mean-square displacement of trajectories.

Assign one file per process, assuming there are more processes than files.
There must be at least one process dedicated to collating data and plotting.
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from mpi4py import MPI
from argparse import ArgumentParser
from packages.trajectories import read_trajectory, get_msd

__author__ = "Johan G. McQuillan"
__email__ = "johan.mcquillan.13@ucl.ac.uk"

# MPI initialisation
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
name = MPI.Get_processor_name()
assert size > 1     # Only works for more than one core

parser = ArgumentParser()
parser.add_argument('traj', type=str, nargs='+', help='Trajectory files')
parser.add_argument('--start', type=int, default=1, help='Starting frame')
parser.add_argument('--end', type=int, default=-1, help='End frame (negative for no end)')
parser.add_argument('--stride', type=int, default=1, help='Frame stride length')
parser.add_argument('-S', type=str, default='O', nargs='+', help='Atomic species')
parser.add_argument('--delta-t', type=float, default=0.2, help='Time step in femtoseconds')
args = parser.parse_args()

# Initialise
throbber = {0:  '-',
            1:  '\\',
            2:  '|',
            3:  '/',
            4:  '-',
            5:  '\\',
            6:  '|',
            7:  '/'}
throb_damp = 5E2

labels = []
dt = 5*1E-4     # ps
d_regime = 100  # ps

# Include log-log plot for determining transport regimes
log_plot = True
if log_plot:
    msd_cols = 1
else:
    msd_cols = 2

if rank == 0:
    # Initialise plot
    ax1 = plt.subplot2grid((2, 2), (0, 0), colspan=msd_cols)
    if log_plot:
        ax2 = plt.subplot2grid((2, 2), (0, 1))
    ax3 = plt.subplot2grid((2, 2), (1, 0), colspan=2)

    ax1.set_xlabel(r'$t$ [ps]')
    ax1.set_ylabel(r'$\langle r^2 (t) \rangle$  [$\AA^2$]')
    ax1.set_xlim([0, 800])
    ax1.set_ylim([0, 20])

    if log_plot:
        ax2.set_xlabel(r'$t$ [ps]')
        ax2.set_xscale('log')
        ax2.set_yscale('log')
        ax2.set_xlim([1E-2, 1E3])
        ax2.set_ylim([1E-3, 1E2])

    ax3.set_ylabel(r'$D$ [$\AA^2 ps^{-1}$]')
    ax3.set_xlabel(r'$T$ [K]')
    ax3.set_ylim([0, 0.02])

for i in range(len(args.traj)):
    traj = args.traj[i]
    
    # Use filename without extension as label
    # Files should ideally be named by temperature
    T = os.path.splitext(os.path.basename(traj))[0]
    
    # For all non-master processes, calculate msd
    if rank != 0:
        # Assign multiple files to one process if too few processes
        if (i - rank) % (size - 1) == 0 and rank <= len(args.traj):
            positions, steps = read_trajectory(traj, start=args.start, end=args.end,
                                               stride=args.stride, axes=(0, 1), species=args.S)
            time = steps * args.delta_t
            msd, a, a_error, b, b_error = get_msd(positions, time, d_regime)

            # Diffusion is half of gradient for diffusion in one dimension
            #  (or diffusion in two dimensions, where one particle in two dimensions
            #   is treated as two particles in one dimension)
            D = b/2
            De = b_error/2
            sys.stdout.write('{}K:  {:.2e} +/- {:.2e}\n'.format(T, D, De))
            sys.stdout.flush()
            
            # Send data to master process
            comm.isend(np.array([time, msd]), dest=0, tag=2*i)
            comm.isend(D, dest=0, tag=2*i+1)
    else:
        # Receive and plot data
        time, msd = comm.irecv(2**20, tag=2*i).wait()
        D = comm.irecv(tag=2*i+1).wait()
        
        ax1.plot(time, msd, label=r'{}K'.format(T))
        if log_plot:
            ax2.plot(time, msd, label=r'{}K'.format(T))
        ax3.scatter(T, D, marker='x', color='k')
        # ax.plot(steps[intervals]*dt, msd_fft(z_coords[:, 0])*2)
        # ax.errorbar(steps[intervals]*dt, msd[], label=r'${}$K'.format(pdb[-8:-5]))
        # ax.plot(time, a+b*time, alpha=0.3, color='k')
        # labels.append(pdb[:4])

# Show plot
if rank == 0:
    if log_plot:
        ax2.legend(loc=1)
    else:
        ax1.legend(loc=1)
    plt.savefig('msd.png')
    plt.show()

