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

__author__ = "Johan G. McQuillan"
__email__ = "johan.mcquillan.13@ucl.ac.uk"

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
        positions (array, float): Coordinates of trajectory.
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
        positions (array, float): Coordinates of trajectory.
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


def read_trajectory(path, **kwargs):
    """Parse .xyz or .pdb trajectories.

    Args:
        path (str): Path to the trajectory.
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
        positions (array, float): Coordinates of trajectory.
            Indexed by [step, atom_index, axis].
        frames (array, int): Time steps.
        cell_array (array, float): Dimensions of supercell.
            Returned only if cell argument is True.
        species_list (list, str): List of atomic species.
            Returned only if spec_list argument is True.
    """
    # Check extension to decide which method to use
    ext = os.path.splitext(path)[1]
    if ext == '.pdb':
        return read_pdb(path, **kwargs)
    elif ext == '.xyz':
        return read_xyz(path, **kwargs)
    else:
        raise IOError(os.path.basename(path)+' - File must be a pdb or xyz')


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

                # Format positions to pdb's very specific length
                x_str = '{:.3f}'.format(x)[:7]
                y_str = '{:.3f}'.format(y)[:7]
                z_str = '{:.3f}'.format(z)[:7]
                pdb.write('ATOM  {:>5} {:>3}    1   {:>3}    {:<7} {:<7} {:<7} 1.00 0.00\n'.format(j+1, spec, 1, x_str, y_str, z_str))
            pdb.write('ENDMDL\n')


def write_xyz(path, positions, cell, species, frames=None):
    """Write trajectory to .xyz file.

    Args:
        path (str): Path to file.
        positions (ndarray, float): Atomic coordinates.
            Indexed by [step, atom_index, dimension].
        cell (ndarray, float): Cell dimensions.
        species (list, str): Atomic species.
        frames (ndarray, int): If given, labels frames with step number.
            Defaults to None, so frames are labelled in increments of 1.
    """
    if os.path.splitext(path)[1] != '.xyz':
        path += '.xyz'

    with open(path, 'w') as xyz:
        # Iterate over frames
        for i in range(positions.shape[0]):
            if frames is not None:
                step = frames[i]
            else:
                step = i

            # Print number of z_coords
            xyz.write('{}\n'.format(positions.shape[1]))

            # Print comment line in i-PI format
            a, b, c = cell[i]
            xyz.write('# CELL(abcABC):   {:<9}   {:<9}   {:<9}    90.00000    90.00000    90.00000  Step: {:>11}  Bead:      0 positions{{angstrom}}  cell{{angstrom}}\n'.format(a, b, c, step))
            for j in range(positions.shape[1]):
                s = species[j]
                x, y, z = positions[i, j]

                # Though xyz positions only need to be separated by whitespace some
                #  programs can only pass them if they fit into specific columns,
                coords = ''
                for u in [x, y, z]:
                    if u > 0:
                        coords += ' '
                    coords += ' {:10.5e}'.format(u)
                xyz.write('{:3}         {}  \n'.format(s, coords))


def write_trajectory(path, *args, **kwargs):
    """Write trajectory to .pdb or .xyz file.

    Args:
        path (str): Path to file.
        positions (ndarray, float): Atomic coordinates.
            Indexed by [step, atom_index, dimension].
        cell (ndarray, float): Cell dimensions.
        species (list, str): Atomic species.
        frames (ndarray, int): If given, labels frames with step number.
            Defaults to None, so frames are labelled in increments of 1.
    """
    # Check extension to decide which method to use
    ext = os.path.splitext(path)[1]
    if ext == '.pdb':
        write_pdb(path, *args, **kwargs)
    elif ext == '.xyz':
        write_xyz(path, *args, **kwargs)
    else:
        raise IOError(os.path.basename(path)+' - File must be a pdb or xyz')


def get_rdf(positions, cell, dr=0.1, r_max=10., prog_bar=False, verbose=False):
    """Return a the 2D radial distribution function of a trajectory.

    Args:
        positions (ndarray, float): Atomic coordinates.
            Indexed by [step, dimension, atom].
        cell (ndarray, float): Cell dimensions.
        dr (float): Bin width for histogram. Defaults to 0.1.
        r_max (float): Maximum range to consider. Defaults to 10.
        prog_bar (bool): If True, show a progress bar.
            Defaults to False.
        verbose (bool): If True, display verbose output.

    Returns:
        bins (ndarray, float): Distances.
        histogram (ndarray, float): RDF histogram.
        i_max1 (int): Index of first maximum.
        i_max2 (int): Index of second maximum.
        coordination (float): Coordination number.
            Integral of RDF up to second maximum.
    """
    # Get average surface number density for each snapshot
    density = positions.shape[2] / (cell[:, 0] * cell[:, 1])
    
    # Get displacements between every pair
    xy = np.zeros((positions.shape[0], positions.shape[1], positions.shape[2], positions.shape[2]))
    for i in range(positions.shape[2]):
        xy[:, :2, i] = positions[:, :2] - positions[:, :2, i, None]
    
    # Shift displacements to nearest pbc image
    while np.any(xy >= cell[:, :, None, None] / 2):
        xy -= cell[:, :, None, None] * (xy >= cell[:, :, None, None] / 2).astype(int)
    while np.any(xy < -cell[:, :, None, None] / 2):
        xy += cell[:, :, None, None] * (xy < -cell[:, :, None, None] / 2).astype(int)
    
    # Calculate absolute distances
    r = np.sqrt(xy[:, 0, :, :] ** 2 + xy[:, 1, :, :] ** 2)
    
    # Initialise histogram
    if verbose or prog_bar:
        sys.stdout.write('Filling Histogram...\n')
        sys.stdout.flush()
    bins = np.linspace(dr / 2, r_max + dr / 2, int(r_max / dr), endpoint=False)
    hist = np.zeros((density.shape[0], bins.shape[0]))
    
    # Fill bins, ignoring zeroth
    bars_done = 0
    for i in range(1, bins.shape[0]):
        hist[:, i] += (1. / (2 * np.pi * i * dr ** 2) *
                       np.sum((r >= i * dr).astype(int) * (r < (i + 1) * dr).astype(int),
                              axis=(1, 2)))
        
        # Update progress bar
        progress = float(i) / bins.shape[0] * bars_total
        if prog_bar and progress >= bars_done:
            bars_done = int(np.round(progress))
            sys.stdout.write('\r  [{:<{}}] '.format(bar_char * bars_done, bars_total))
            sys.stdout.write(' {:>3.0f}%'.format(np.ceil(float(i) / bins.shape[0] * 100)))
            sys.stdout.flush()
    if prog_bar:
        sys.stdout.write('\n')
        sys.stdout.flush()
    
    # Normalise rdf by density
    hist /= density[:, None] * (positions.shape[2] - 1)
    # Average over all frames
    histogram = np.sum(hist, axis=0) / hist.shape[0]
    
    # Get first maximum (assume it is the highest peak)
    i_max1 = np.argmax(histogram)
    
    minima_2_found = False
    i = i_max1 + 1
    # Find second minimum
    # Loop until gradient is positive
    while not minima_2_found:
        if histogram[i] > histogram[i - 1]:
            i_min2 = i - 1
            minima_2_found = True
        else:
            i += 1
    
    # Calculate coordination number by integrating from 0 to i_min2
    coordination = 0
    for i in range(0, i_min2):
        coordination += histogram[i] * bins[i] * dr
    # Unnormalise
    coordination *= 2 * np.pi * np.mean(density)
    
    # Find second maxima
    maxima_2_found = False
    i = i_min2 + 1
    # Loop until gradient is negative
    while not maxima_2_found:
        if histogram[i] < histogram[i - 1]:
            i_max2 = i - 1
            maxima_2_found = True
        else:
            i += 1
    
    return bins, histogram, i_max1, i_max2, coordination


def reduce_to_supercell(positions, cell, molecules=False):
    """Use periodic boundary conditions to ensure all coordinates
    are within the first supercell.

    Args:
        positions (ndarray, float): Atomic coordinates.
            Indexed by [step, atom_index, dimension].
        cell (ndarray, float): Cell dimensions.
        molecules (bool): If True, then keep molecules together.
            Reduce oxygen positions to first supercell, then translate
            hydrogen positions to stay connected, even if out of supercell.
            Assumes species order in array is [O, H, H, O, H, H, ...].

    Returns:
        new_positions (ndarray, float): Reduced atomic coordinates.
    """
    new_positions = positions
    # Pad cell array for every atom
    cell_array = cell[:, None, :]
    if not molecules:
        # Reduce all atomic new_positions
        while np.any(new_positions >= cell_array):
            new_positions -= cell_array * (new_positions >= cell_array).astype(int)
        while np.any(new_positions < 0.):
            new_positions += cell_array * (new_positions < 0.).astype(int)
    else:
        # Only check for oxygens out of cell
        while np.any(new_positions[:, ::3] >= cell_array):
            # Find which oxygens are out of cell
            truth_grid = new_positions[:, ::3] >= cell_array

            # Duplicate above grid to apply to each group of OHH
            mask = np.zeros(new_positions.shape, dtype=bool)
            mask[:, 0::3] = truth_grid
            mask[:, 1::3] = truth_grid
            mask[:, 2::3] = truth_grid

            new_positions -= cell_array * mask.astype(int)
        while np.any(new_positions[:, ::3] < 0.):
            # Find which oxygens are out of cell
            truth_grid = new_positions[:, ::3] < 0.

            # Duplicate above grid to apply to each group of OHH
            mask = np.zeros(new_positions.shape, dtype=bool)
            mask[:, 0::3] = truth_grid
            mask[:, 1::3] = truth_grid
            mask[:, 2::3] = truth_grid
            new_positions += cell_array * mask.astype(int)
    return new_positions


def plane_distribution(hydrogens, bins=50):
    """Return histogram of the plane distributions of water molecules.

    Finds the distribution of the angle that the vector from one H to other H
    in a molecule makes with the z axis. 0 degrees means both H are in
    the x-y plane; 90 degrees means both H are out of plane.

    Args:
        hydrogens (ndarray, float): Coordinates of hydrogen z_coords.
            Must be even number of z_coords and consecutive z_coords must
            pertain to same molecule.
            Indexed by [step, atom, dimensions]
        bins (int): Number of bins for histogram. Defaults to 50.
    Returns:
        hist (ndarray, float): Normalised frequency.
        degrees (ndarray, float): Angle in degrees.
    """
    if hydrogens.shape[1] % 2 != 0:
        raise ValueError('Odd number of hydrogens given. Have you included oxygens?')
    
    # Split each pair of hydrogen into two arrays
    h1 = hydrogens[:, 0::2]
    h2 = hydrogens[:, 1::2]
    
    # Get the unit vector from one H to other
    d = h2[:, :] - h1[:, :]
    d_mag = np.linalg.norm(d, axis=2)
    d_norm = d / d_mag[:, :, None]
    
    # Calculate the angle with z axis
    theta = 90. - np.arccos(d_norm[:, :, -1])*180./np.pi
    
    # Return histogram using absolute value of theta as hydrogens are identical
    degrees, hist = np.histogram(abs(theta), bins=bins, normed=True)
    return degrees, hist
