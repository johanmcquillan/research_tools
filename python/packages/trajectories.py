"""Trajectory analysis tools for i-PI simulations

This module provides methods for analysing .pdb or .xyz trajectories
from i-PI simulations (https://github.com/cosmo-epfl/i-pi-dev), as well
as output files (that is, files containing the thermodynamic
properties).

This is developed as part of my PhD and so the methods of analysis
implemented here are geared towards applications regarding quasi-two-
dimensional water.
"""

import numpy as np
import os
import sys
import warnings
from collections import Iterable


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
bars_total = 20
bar_char = '>'


def read_pdb(pdb_path, start=0, end=-1, stride=1, cell=False, spec_list=False,
             species=('O',), axes=(0, 1), throb=False):
    """Parse .pdb trajectory.

    Args:
        pdb_path (str): Path to the .pdb trajectory.
        start (int, opt): Starting frame (time step). Defaults to 0.
        end (int, opt): End frame (time step).
            A negative value means do not stop until EOF.
            Defaults to -1.
        stride (int): Number of time steps to skip between saving.
            Defaults to 1.
        cell (bool): If True, return cell dimensions. Defaults to False.
        spec_list (bool): If True, return corresponding list of
            atomic species. Defaults to False.
        species (tuple, str): Species to read. Defaults to ('O',).
        axes (tuple, int): Axes to read.
            0, 1, 2 for x, y, z respectively. Defaults to all.
        throb (bool): If True, display a throbber. Defaults to True.

    Returns:
        coords (array, float): Coordinates of trajectory.
            Indexed by [step, atom_index, axis].
        frames (array, int): Time steps.
        cell_array (array, float): Dimensions of supercell.
            Returned only if cell argument is True.
        species_list (list, str): List of atomic species.
            Returned only if spec_list argument is True.
    """
    # Initialise
    coords = []
    frames = []
    cells = []
    species_list = []
    warned = False
    first = True
    
    # Check axes argument is either an int or an Iterable of ints between 0
    # and 2
    if isinstance(axes, int):
        if 0 <= axes <= 2:
            axes = [axes]
        else:
            raise ValueError('Axis must be 0, 1, or 2 (for x, y, or z)')
    elif isinstance(axes, Iterable):
        if not all([isinstance(x, int) and 0 <= x <= 2 for x in axes]):
            raise ValueError('Axes must be 0, 1, or 2 (for x, y, or z)')
    else:
        raise ValueError('Axes must be an integer or list of integers')
    # Convert to ndarray array for advanced indexing
    axes = np.array(axes, dtype=int)
    
    if 0 <= end < start:
        raise ValueError('Starting frame must be the same as or before the final frame')
    
    if isinstance(species, list) and 'H' in species:
        for H in ['H1', 'H2']:
            if H not in species:
                species.append(H)
    elif isinstance(species, str):
        species = species.split(' ')
    elif not isinstance(species, Iterable):
        species = [str(species)]
    
    # Initialise throbber
    if throb:
        sys.stdout.write('Reading {}   '.format(pdb_path))
        sys.stdout.flush()
    
    # Begin parsing file
    with open(pdb_path, 'r') as pdb:
        # Placeholder variables
        step = 0
        line = 'Null'
        line_split = [line]
        try:
            # Skip lines until we get to start frame
            while step < start:
                line = pdb.next()
                line_split = line.split()
                # Skip lines until we get to TITLE line, which lists the step
                #  number
                while line_split[0] != 'TITLE':
                    line = pdb.next()
                    line_split = line.split()
                step = int(line_split[2])
            
            last_added_step = step - stride
            # Start saving data
            # Loop ends when step > end or EOF
            while step <= end or end < 0:
                # Update throbber
                if throb:
                    sys.stdout.write('\rReading {}  {} '.format(pdb_path, throbber[int(step / throb_damp % 8)]))
                    sys.stdout.flush()
                
                # Skip lines until we get to cell dimension data
                while line_split[0] != 'CRYST1':
                    line = pdb.next()
                    line_split = line.split()
                
                # Parse cell dimensions and angles (latter not saved)
                cell_dims = np.array([float(x) for x in [line[6 + y * 9:6 + (y + 1) * 9] for y in range(3)]])
                cell_angs = np.array([float(x) for x in [line[33 + y * 7:33 + (y + 1) * 7] for y in range(3)]])
                
                line = pdb.next()
                line_split = line.split()
                
                # Save frame if stride reached
                diff = step - last_added_step
                if diff >= stride:
                    # If difference between frames larger than stride, warn user and continue
                    if diff != stride and not warned:
                        warned = True
                        warnings.warn(
                                'From step {}, the stride length may not be equal.'.format(step) +
                                'BE CAREFUL')
                    last_added_step = step
                    add_step = True
                    frames.append(step - start)
                    coords.append([])
                    if cell:
                        cells.append(cell_dims[axes])
                else:
                    add_step = False
                
                # Parse coordinates for each atom
                while line_split[0] != 'END':
                    spec = line_split[2]
                    if add_step and spec in species:
                        coords[-1].append([float(x) for x in [line[30 + y * 8:30 + (y + 1) * 8] for y in axes]])
                        if spec_list and first:
                            species_list.append(spec)
                    line = pdb.next()
                    line_split = line.split()
                if add_step:
                    first = False
                
                line = pdb.next()
                line_split = line.split()
                step = int(line_split[2])
        except (EOFError, StopIteration):
            pass
        finally:
            # Clean throbber
            if throb:
                sys.stdout.write('\rReading {}   \n'.format(pdb_path))
                sys.stdout.flush()
    
    # Cast to ndarray and return appropriate arrays
    coords = np.array(coords)
    frames = np.array(frames)
    if not cell:
        if not spec_list:
            return coords, frames
        else:
            return coords, frames, cells
    else:
        cells = np.array(cells)
        if not spec_list:
            return coords, frames, cells
        else:
            return coords, frames, cells, species_list


def read_xyz(xyz_path, start=0, end=-1, stride=1, cell=False, spec_list=False,
             species=('O',), axes=(0, 1), throb=False):
    """Parse .xyz or .pdb trajectories.

    Args:
        xyz_path (str): Path to the .xyz trajectory.
        start (int, opt): Starting frame (time step). Defaults to 0.
        end (int, opt): End frame (time step).
            A negative value means do not stop until EOF.
            Defaults to -1.
        stride (int): Number of time steps to skip between saving.
            Defaults to 1.
        cell (bool): If True, return cell dimensions. Defaults to False.
        spec_list (bool): If True, return corresponding list of
            atomic species. Defaults to False.
        species (tuple, str): Species to read. Defaults to ('O',).
        axes (tuple, int): Axes to read.
            0, 1, 2 for x, y, z respectively. Defaults to all.
        throb (bool): If True, display a throbber. Defaults to True.

    Returns:
        coords (array, float): Coordinates of trajectory.
            Indexed by [step, atom_index, axis].
        frames (array, int): Time steps.
        cell_array (array, float): Dimensions of supercell.
            Returned only if cell argument is True.
        species_list (list, str): List of atomic species.
            Returned only if spec_list argument is True.
    """
    # Initialise
    coords = []
    frames = []
    cells = []
    species_list = []
    warned = False
    first = True

    # Check axes argument is either an int or an Iterable of ints between 0 and 2
    if isinstance(axes, int):
        if 0 <= axes <= 2:
            axes = [axes]
        else:
            raise ValueError('Axis must be 0, 1, or 2 (for x, y, or z)')
    elif isinstance(axes, Iterable):
        if not all([isinstance(x, int) and 0 <= x <= 2 for x in axes]):
            raise ValueError('Axes must be 0, 1, or 2 (for x, y, or z)')
    else:
        raise ValueError('Axes must be an integer or list of integers')
    # Convert to ndarray array for advanced indexing
    axes = np.array(axes, dtype=int)
    
    if 0 <= end < start:
        raise ValueError('Starting frame must be the same as or after the final frame')
    
    if isinstance(species, list) and 'H' in species:
        for H in ['H1', 'H2']:
            if H not in species:
                species.append(H)
    elif isinstance(species, str):
        species = species.split(' ')
    elif not isinstance(species, Iterable):
        species = [str(species)]

    # Initialise throbber
    if throb:
        sys.stdout.write('Reading {}   '.format(xyz_path))
        sys.stdout.flush()

    # Begin parsing file
    with open(xyz_path, 'r') as xyz:
        # Placeholder variables
        step = -1
        line = 'Null'
        line_split = [line]
        try:
            # Skip lines until we to start frame
            while step < start:
                # Skip lines until we get to comment line, which detail step number
                while line_split[0] != '#':
                    line = xyz.next()
                    line_split = line.split()
                step = int(line_split[9])
                if step < start:
                    line = xyz.next()
                    line_split = line.split()

            last_added_step = step - stride
            # Start saving data
            # Loop ends when step > end or EOF
            while step <= end or end < 0:
                if throb:
                    sys.stdout.write('\rReading {} {} '.format(xyz_path, throbber[int(step/throb_damp % 8)]))

                # Parse cell dimensions and angles (latter not saved)
                cell_dims = np.array([float(x) for x in np.array(line_split)[2+axes]])
                cell_angs = np.array([float(x) for x in np.array(line_split)[5+axes]])

                line = xyz.next()
                line_split = line.split()

                # Save frame if stride reached
                diff = step - last_added_step
                if diff >= stride:
                    # If difference between frames larger than stride, warn user and continue
                    if diff != stride and not warned:
                        warned = True
                        warnings.warn('From step {}, the stride length may not be equal. BE CAREFUL'.format(step))
                    last_added_step = step
                    add_step = True
                    frames.append(step - start)
                    coords.append([])
                    if cell:
                        cells.append(cell_dims[axes])
                else:
                    add_step = False

                # Parse coordinates for each atom
                while len(line_split) < 5:
                    spec = line_split[0]
                    if add_step and spec in species:
                        coords[-1].append([float(x) for x in np.array(line_split)[1+axes]])
                        if spec_list and first:
                            species_list.append(spec)
                    line = xyz.next()
                    line_split = line.split()
                if add_step:
                    first = False
                step = int(line_split[9])
        except (EOFError, StopIteration):
            pass
        except IndexError:
            if line.split() == ['']:
                pass
            else:
                raise IOError('Problem reading step {}'.format(step))
        finally:
            # Clean throbber
            sys.stdout.write('\rReading {}   \n'.format(xyz_path))
            sys.stdout.flush()

    # Cast to ndarray and return appropriate arrays
    coords = np.array(coords)
    frames = np.array(frames)
    if not cell:
        if not spec_list:
            return coords, frames
        else:
            return coords, frames, cells
    else:
        cells = np.array(cells)
        if not spec_list:
            return coords, frames, cells
        else:
            return coords, frames, cells, species_list


def write_pdb(path, positions, cell, species, frames=None):
    """Write trajectory to .pdb file.

    Args:
        path (str): Path to file.
        positions (ndarray, float): Atomic coordinates.
            Indexed by [step, atom_index, dimension].
        cell (ndarray, float): Cell dimensions.
        species (list, str): Atomic species.
        frames (ndarray, int): If given, labels frames with step number.
            Defaults to None, so frames are labelled in increments of 1.
    """
    if os.path.splitext(path)[1] != '.pdb':
        path += '.pdb'

    with open(path, 'w') as pdb:
        # Iterate over frames
        for i in range(positions.shape[0]):
            if frames is not None:
                step = frames[i]
            else:
                step = i

            # Write header
            pdb.write('TITLE     Step: {:>11d}  Bead:       0 positions{{angstrom}}  cell{{angstrom}}\n'.format(step))
            pdb.write('CRYST1 {:<8.5f} {:<8.5f} {:<8.5f} 90.00  90.00  90.00 P 1           1\n'.format(*cell[i]))

            # Iterate over z_coords
            for j in range(positions.shape[1]):
                spec = species[j]
                x, y, z = positions[i, j]

                # Format coords to pdb's very specific length
                x_str = '{:.3f}'.format(x)[:7]
                y_str = '{:.3f}'.format(y)[:7]
                z_str = '{:.3f}'.format(z)[:7]
                pdb.write('ATOM  {:>5} {:>3}    1   {:>3}    {:<7} {:<7} {:<7} 1.00 0.00\n'.format(j+1, spec, 1, x_str, y_str, z_str))
            pdb.write('ENDMDL\n')
