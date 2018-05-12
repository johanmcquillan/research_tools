#!/usr/bin/env python2.7

from argparse import ArgumentParser
from packages import filer, ipi_builder

__author__ = "Johan G. McQuillan"
__email__ = "johan.mcquillan.13@ucl.ac.uk"

# Parse arguments
parser = ArgumentParser()
parser.add_argument('prefix', type=str, help='Prefix for i-PI output files')
parser.add_argument('xyz_path', metavar='xyz', type=str, help='Path to .xyz input')
parser.add_argument('h', metavar=('a', 'b', 'c'), type=float, nargs=3, help='Cell dimensions in Angstroms')
parser.add_argument('ip', type=str, nargs='?', default=None, help='IP address if using INET port')
parser.add_argument('wport', type=str, nargs='?', default='water', help='INET port or UNIX address for water-water force clients')
parser.add_argument('xport', type=str, nargs='?', default=None, help='INET port or UNIX address for external potential')
parser.add_argument('alpha', type=float, nargs='?', default=0.29, help='Screening parameter for Ewald summation')
parser.add_argument('k', type=list, nargs='?', default=[12, 12, 1], help='Max K-points for Ewald summation')
args = parser.parse_args()

# Parse IP, port, and vector variables
if args.ip == 'None':
    args.ip = None
if args.xport == 'None':
    args.xport = None
args.h = [float(x) for x in args.h]
args.k = [int(k) for k in args.k]

# Get number of molecules
n = ipi_builder.get_n_from_xyz(args.xyz_path)

input_file = 'nvt.xml'
field_file = 'FIELD'
config_file = 'CONFIG.01'
control_file = 'CONTROL'

# Write files
with filer.safe_open(input_file, 'w') as f:
    contents = ipi_builder.mbpol_xml_min(n, args.prefix, args.xyz_path, args.h, args.ip, args.wport, args.xport)
    f.write(contents)
with filer.safe_open(field_file, 'w') as f:
    contents = ipi_builder.mbpol_field(n)
    f.write(contents)
with filer.safe_open(config_file, 'w') as f:
    contents = ipi_builder.mbpol_config(n, args.h)
    f.write(contents)
with filer.safe_open(control_file, 'w') as f:
    contents = ipi_builder.mbpol_control(args.ip, args.wport, args.alpha, args.k)
    f.write(contents)

