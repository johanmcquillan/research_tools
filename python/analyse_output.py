#!/usr/bin/env python2.7

import matplotlib.pyplot as plt
from argparse import ArgumentParser
from packages.trajectories import read_output

__author__ = "Johan G. McQuillan"
__email__ = "johan.mcquillan.13@ucl.ac.uk"

parser = ArgumentParser()
parser.add_argument('path', type=str, help='Path to i-PI output file')
parser.add_argument('x', type=str, help='Quantity for x-axis')
parser.add_argument('y', type=str, help='Quantity for y-axis')
args = parser.parse_args()

data, units = read_output(args.path)

# Check if valid quantities
for quantity1 in [args.x, args.y]:
    if quantity1 not in data:
        error_string = 'The quantity "{}" was not recorded (check spelling or input file).\n'.format(
                quantity1)
        error_string += '    The quantities recorded in {} include:'.format(args.path)
        for quantity2 in sorted(data):
            error_string += '\n        {}'.format(quantity2)
        raise LookupError(error_string)

x = data[args.x]
y = data[args.y]

plt.plot(x, y)
plt.show()
