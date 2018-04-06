#!/usr/bin/env python

"""Calculate the velocity density of states and the power spectrum from a velocity trajectory"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from argparse import ArgumentParser, RawTextHelpFormatter
from packages.cli_animations import update_bar
from packages.trajectories import read_trajectory, autocorrelation

parser = ArgumentParser(formatter_class=RawTextHelpFormatter)
parser.add_argument('vtraj', type=str, nargs='+', help='Veloctiy trajectory files')
parser.add_argument('--start', type=int, default=0, help='Starting frame')
parser.add_argument('--end', type=int, default=-1, help='End frame (negative for all frames)')
parser.add_argument('--stride', type=int, default=1, help='Frame stride')
parser.add_argument('-d', metavar='dr', type=float, default=0.1, help='Width of histogram bin in Angstroms')
parser.add_argument('-m', metavar='r_max', type=float, default=10., help='Maximum distance in Angstroms')
parser.add_argument('-v', action='store_true', help='Verbose')
args = parser.parse_args()

fig, axes = plt.subplots(2)
ax1, ax2 = axes

dt = 0.2*10*1E-3    # Time step in femtoseconds

# Mass of elements in kg
mass = {'O':    2.6566962E-26,
        'H':    1.67E-30}
mass['H2O'] = mass['O'] + 2*mass['H']

for traj in args.vtraj:
    velocities, frames = read_trajectory(traj, start=args.start, end=args.end, stride=args.stride,
                                         axes=(0, 1, 2), species=('O', 'H'), throb=args.v)
    norm_atoms = velocities / np.sum(velocities**2, axis=2)[..., None]
    
    if args.v:
        sys.stdout.write('{} Particles over {} Frames\n'.format(velocities.shape[2],
                                                                velocities.shape[0]))
        sys.stdout.write('Calculating autocorrelation functions...\n')
        sys.stdout.flush()
    
    # Number of each atomic species
    N = {'O': velocities.shape[1]/3,
         'H': velocities.shape[1]*2/3,
         'total': velocities.shape[1]}
    
    # Split velocity array by species
    vels = {}
    vels['O'] = velocities[:, ::3]
    vels['H'] = np.zeros((velocities.shape[0], N['H'], velocities.shape[2]))
    vels['H'][:, 0::2] += velocities[:, 1::3]
    vels['H'][:, 1::2] += velocities[:, 2::3]
    
    # Initialise ACF
    vacf = np.zeros(velocities.shape[0], dtype=np.complex)
    
    # Calculate ACF for each species and dimension and combine into total ACF
    bars_done = 0
    for s in vels:
        vel = vels[s]
        # Weight by the relative mass
        m = mass[s] / mass['H2O']
        for i in range(N[s]):
            xacf = autocorrelation(vel[:, i, 0], frames)
            yacf = autocorrelation(vel[:, i, 1], frames)
            zacf = autocorrelation(vel[:, i, 2], frames)
            vacf += m*(xacf + yacf + zacf)
            if args.v:
                update_bar(i, N[s], bars_done)
        if args.v:
            sys.stdout.write('\n')
            sys.stdout.flush()
    
    # Normalise by number of atoms
    vacf /= N['total']
    
    # Get the total velocity density of states and power spectrum
    vdos = np.fft.fft(vacf)
    power_spectrum = vdos * np.conj(vdos)
    
    # Get x axis arrays
    time = frames*dt
    freq = np.fft.fftfreq(vdos.shape[0])/dt

    # Plot results
    label = os.path.splitext(os.path.basename(traj))[0]
    ax1.plot(time, np.real(vacf), label=label)
    ax2.plot(freq, power_spectrum)

ax1.set_xlim([0, 10])
ax1.set_ylabel(r'VACF$(t)$')
ax2.set_xlabel(r'$t$ [ps]')

if len(args.vtraj) > 1:
    ax1.legend()
plt.show()

if args.v:
    sys.stdout.write('Done.\n')
    sys.stdout.flush()
