#!/usr/bin/env python2.7

from argparse import ArgumentParser
from packages import filer, ipi_builder

# Parse arguments
parser = ArgumentParser()
parser.add_argument('prefix', type=str,
                    help='Prefix for i-PI output files')
parser.add_argument('xyz_path', metavar='xyz', type=str,
                    help='Path to .xyz input')
parser.add_argument('a', type=float,
                    help='Cell X dimensions in Angstroms')
parser.add_argument('b', type=float,
                    help='Cell Y dimensions in Angstroms')
parser.add_argument('c', type=float,
                    help='Cell Z dimensions in Angstroms')
parser.add_argument('s', type=int,
                    help='Total time steps')
parser.add_argument('T', type=float,
                    help='Temperature in Kelvin')
parser.add_argument('dt', type=float, default=0.2,
                    help='Time step in femtoseconds')
parser.add_argument('tau', type=float, nargs='?', default=50.0,
                    help='Langevin time constant in femtoseconds')
parser.add_argument('S', type=int, nargs='?', default=1,
                    help='Stride for printing frame data')
parser.add_argument('ip', type=str, nargs='?', default=None,
                    help='IP address if using INET port')
parser.add_argument('wport', type=str, nargs='?', default='water',
                    help='INET port or UNIX address for water-water force clients')
parser.add_argument('xport', type=str, nargs='?', default=None,
                    help='INET port or UNIX address for external potential')
parser.add_argument('--alpha', type=float, nargs='?', default=0.29,
                    help='Screening parameter for Ewald summation')
parser.add_argument('-k', type=list, nargs='?', default=[12, 12, 1],
                    help='Max K-points for Ewald summation')
parser.add_argument('--no-pdb', dest='pdb', action='store_false',
                    help='Do not output PDB trajectory')
parser.add_argument('--vhigh', dest='high', action='store_true',
                    help='Set i-PI verbosity to high')
args = parser.parse_args()

# Parse IP, port, and vector variables
if args.ip == 'None':
    args.ip = None
if args.xport == 'None':
    args.xport = None
h = [args.a, args.b, args.c]
args.k = [int(k) for k in args.k]

# Get number of molecules
n = ipi_builder.get_n_from_xyz(args.xyz_path)

input_file = 'nvt.xml'
field_file = 'FIELD'
config_file = 'CONFIG.01'
control_file = 'CONTROL'

# Write files
with filer.safe_open(input_file, 'w') as f:
    contents = ipi_builder.mbpol_xml_nvt(n, args.prefix, args.xyz_path, h, args.s, args.T, args.dt,
                                         args.tau, args.ip, args.wport, args.xport, stride=args.S,
                                         out_pdb=args.pdb, vhigh=args.high)
    f.write(contents)
with filer.safe_open(field_file, 'w') as f:
    contents = ipi_builder.mbpol_field(n)
    f.write(contents)
with filer.safe_open(config_file, 'w') as f:
    contents = ipi_builder.mbpol_config(n, h)
    f.write(contents)
with filer.safe_open(control_file, 'w') as f:
    contents = ipi_builder.mbpol_control(args.ip, args.wport, args.alpha, args.k)
    f.write(contents)
