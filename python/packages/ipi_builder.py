#!/usr/bin/env python

__author__ = "Johan G. McQuillan"
__email__ = "johan.mcquillan.13@ucl.ac.uk"


def mbpol_xml_nvt(n, prefix, structure_file, h, s, temp, dt, tau, ip, wport, xport, stride=50,
                  out_pdb=True, out_vel=True, vhigh=False):
    """Generate i-PI XML input file contents for a NVT MB-pol simulation.
    
    Args:
        prefix (str): Prefix for output file names.
        structure_file (str): Path to input co-ordinate file.
        h (iterable(float)): Super cell dimensions in Angstroms (assume orthorhombic cell).
        s (int): Total number of steps.
        temp (float): Temperature in Kelvin.
        dt (float): Duration of time step in femtoseconds.
        tau (float): Time constant for Langevin integrator in femtoseconds.
        ip (str): IP of head node. None if running with UNIX ports.
        wport (int): Port number for water potential calculations.
        xport (int): Port number for external potential calculations.
        stride (int): Time steps to skip before printing data.
        out_pdb (bool): If True, output a PDB as well as XYZ trajectory. Else, just XYZ.
        out_vel (bool): If True, output velocity trajectories.
        vhigh (bool): If True, set i-PI to high verbosity.
    
    Returns:
        str: The XML contents.
    """
    
    ext = structure_file.split('.')[-1]
    fstride = 10
    if ip is not None:
        sock = 'inet'
        address1 = '''<address>{}</address>
    <port>{}</port>'''.format(ip, wport)
        if xport is not None:
            address2 = '''<address>{}</address>
    <port>{}</port>'''.format(ip, xport)
        else:
            address2 = ''
    else:
        sock = 'unix'
        address1 = '<address>{}</address>'.format(wport)
        if xport is not None:
            address2 = '<address>{}</address>'.format(xport)
        else:
            address2 = ''
    sockets = '''<ffsocket mode='{}' name='mbpol'>
    {}
  </ffsocket>'''.format(sock, address1)
    if xport is not None:
        sockets += '''
  <ffsocket mode='{}' name='external'>
    {}
  </ffsocket>'''.format(sock, address2)
        extforce = '\n      <force forcefield=\'external\'/>'
    else:
        extforce = ""
    if out_pdb:
        pdb_line = "<trajectory filename='pdb' stride='{}' format='pdb' cell_units='angstrom'> positions{{angstrom}} </trajectory>\n    ".format(stride)
    else:
        pdb_line = ""
    if out_vel:
        vel_line = "<trajectory filename='velocity' stride='{}' format='xyz' cell_units='angstrom'> velocities{{m/s}} </trajectory>\n    ".format(fstride)
    else:
        vel_line = ""
    if vhigh:
        string = "<simulation verbosity='high'>\n"
    else:
        string = "<simulation verbosity='low'>\n"
    string += '''  <output prefix='{}'>
    <properties filename='out' stride='{}'> [step, time{{picosecond}}, temperature{{kelvin}}, kinetic_md{{millielectronvolt}}, pot_component(index=0){{millielectronvolt}}, pot_component(index=1){{millielectronvolt}}, pressure_md{{pascal}}, stress_md{{pascal}}] </properties>
    {}<trajectory filename='xyz' stride='{}' format='xyz' cell_units='angstrom'> positions{{angstrom}} </trajectory>
    {}<trajectory filename='force' stride='{}' cell_units='angstrom'> forces{{piconewton}} </trajectory>
    <checkpoint stride='1000'/>
  </output>
  <step> 0 </step>
  <total_steps> {} </total_steps>
  {}
  <system>
    <initialize nbeads='1'>
      <file mode='{}' units='angstrom'> {} </file>
      <cell mode='abc' units='angstrom'> [{}, {}, {}] </cell>
      <velocities mode='thermal' units='kelvin'> {} </velocities>
    </initialize>
    <motion mode='dynamics'>
      <dynamics mode='nvt'>
        <timestep units='femtosecond'> {} </timestep>
        <thermostat mode='langevin'>
          <tau units='femtosecond'> {} </tau>
        </thermostat>
      </dynamics>
    </motion>
    <ensemble>
      <temperature units='kelvin'> {} </temperature>
    </ensemble>
    <forces>
      <force forcefield='mbpol'/>{}
    </forces>
  </system>
</simulation>'''.format(prefix, stride, pdb_line, stride, vel_line, fstride, s, sockets, ext, structure_file, h[0], h[1], h[2], temp, dt, tau, temp, extforce)
    return string


def mbpol_xml_min(n, prefix, structure_file, h, ip, wport, xport, stride=1):
    """Generate i-PI XML input file contents for an MB-pol geometry optimisation.

    Args:
        prefix (str): Prefix for output file names.
        structure_file (str): Path to input co-ordinate file.
        h (iterable(float)): Super cell dimensions in Angstroms (assume orthorhombic cell).
        ip (str): IP of head node. None if running with UNIX ports.
        wport (int): Port number for water potential calculations.
        xport (int): Port number for external potential calculations.
        stride (int): Time steps to skip before printing data.

    Returns:
        str: The XML contents.
    """
    
    ext = structure_file.split('.')[-1]
    if ip is not None:
        sock = 'inet'
        address1 = '''<address>{}</address>
    <port>{}</port>'''.format(ip, wport)
        if xport is not None:
            address2 = '''<address>{}</address>
    <port>{}</port>'''.format(ip, xport)
        else:
            address2 = ''
    else:
        sock = 'unix'
        address1 = '<address>{}</address>'.format(wport)
        if xport is not None:
            address2 = '<address>{}</address>'.format(xport)
        else:
            address2 = ''
    sockets = '''<ffsocket mode='{}' name='mbpol'>
    {}
  </ffsocket>'''.format(sock, address1)
    if xport is not None:
        sockets += '''
  <ffsocket mode='{}' name='external'>
    {}
  </ffsocket>'''.format(sock, address2)
        extforce = '\n      <force forcefield=\'external\'/>'
    else:
        extforce = ''
    return '''<simulation verbosity='low'>
  <output prefix='{}'>
    <properties filename='out' stride='{}'> [step, time{{picosecond}}, temperature{{kelvin}}, kinetic_md{{millielectronvolt}}, pot_component(index=0){{millielectronvolt}}, pot_component(index=1){{millielectronvolt}}, pressure_md{{pascal}}, stress_md{{pascal}}] </properties>
    <trajectory filename='pdb' stride='{}' format='pdb' cell_units='angstrom'> positions{{angstrom}} </trajectory>
    <trajectory filename='xyz' stride='{}' format='xyz' cell_units='angstrom'> positions{{angstrom}} </trajectory>
    <trajectory filename='force' stride='{}' cell_units='angstrom'> forces{{piconewton}} </trajectory>
    <checkpoint stride='10000'/>
  </output>
  {}
  <system>
    <initialize nbeads='1'>
      <file mode='{}' units='angstrom'> {} </file>
      <cell mode='abc' units='angstrom'> [{}, {}, {}] </cell>
    </initialize>
    <motion mode='minimize'>
    </motion>
    <forces>
      <force forcefield='mbpol'/>{}
    </forces>
  </system>
</simulation>'''.format(prefix, stride, stride, stride, stride, sockets, ext, structure_file, h[0], h[1], h[2], extforce)


def mbpol_field(n):
    """Return DLPOLY2 FIELD file contents for an MB-pol simulation of n water molecules"""
    
    return '''MBPOL H2O [CONTROL file must include 'vttm4']
units kcal
molecular types 2
H2O
nummols    {}
z_coords      3
O       15.9949    0.000000    1.310       1.310       1
H        1.0079    0.500000    0.294       0.294       1
H        1.0079    0.500000    0.294       0.294       1
bonds            3
mors    1    2  120.876713 0.951960716 2.587949758
mors    1    3  120.876713 0.951960716 2.587949758
expo    2    3  15319.9337 12.66426998
angles     1
schw    2    1    3
finish
M site
nummols    {}
z_coords      1
M        0.0000     -1.000000  1.310       0.0         1
finish
vdw        3
H      H      tt68 20.09358600 9.406475170 0.0         1.0
O      H      tt68 83.49556670 9.775202425 0.0         1.0
O      O      tt68 237.3212214 9.295485815 0.0         1.0
close
'''.format(n, n)


def mbpol_config(n, h):
    """Return DLPOLY2 CONFIG.01 file contents for an MB-pol simulation.
    
    Args:
        n (int): Number of water molecules
        h (interable(float)): Dimensions of super cell. Orthorhmobic cell is assumed.
    
    Returns:
        str: CONFIG.01 contents.
    """
    
    a, b, c = h
    output = '\n         0         2         0      0.00000000    \n'
    output += '       {:0<13}       {:0<13}       {:0<13}\n'.format(a, 0., 0.)
    output += '       {:0<13}       {:0<13}       {:0<13}\n'.format(0., b, 0.)
    output += '       {:0<13}       {:0<13}       {:0<13}\n'.format(0., 0., c)
    for i in range(0, 3*n):
        if i % 3 == 0:
            s = 'O'
        else:
            s = 'H'
        output += '{}       {:>10}\n'.format(s, i+1)
        output += '          {:0<16}          {:0<16}          {:0<16}\n'.format(0., 0., 0.)
    for i in range(0, n):
        output += 'M       {:>10}\n'.format(3*n+i+1)
        output += '          {:0<16}          {:0<16}          {:0<16}\n'.format(0., 0., 0.)
    return output


def mbpol_control(ip, wport, alpha, k):
    """Return DLPOLY2 CONTROL file contents for an MB-pol simulation.
    
    Args:
        ip (str): IP of head node. None if running with UNIX ports.
        wport (int): Port number for water potential calculations.
        alpha (float): Ewald summation screening parameter.
        k (iterable(int)): Maximum k-points for Fourier part of Ewald sum, eg [kx, ky, kz].
    
    Returns:
        str: CONTROL file contents.
    """
    if ip is not None:
        socket_line = 'ipi      inet:{}:{}'.format(ip, wport)
    else:
        socket_line = 'ipi      unix:{}:32848'.format(wport)
    output = '''MBPol potential -- water 

{}

# PREDEFINED MODELS.
vmbpol            ! MB-pol (requires polarizability)

# SYSTEM CUTOFFS
cutoff           9.00        ! electrostatics cutoff
rvdw             9.00        ! nonbonded cutoff
delr width       0.00        ! for Verlet list


# EWALD
ewald sum {} {} {} {}      ! alpha & k-vectors for bulk simulations 
#ewald sum 1.0e-16 0 0 0      ! alpha & k-vectors for clusters simulations 
#ewald precision  1.0E-8     ! dlpoly format two floats to same precision define alpha & k-vectors based on the precision

# POLARIZABLE SIMULATIONS.
polarizability           ! define a polarizable model
aspc                     ! always stable predictor-corrector algorithm
#iteration 1.0e-8         ! iterative (self-consistent) algorithm 
#car-parrinello 5.0e-8    ! Car-Parrinello type propagation for the dipoles 

#***********
#Ignore the following keywords. Do not use the following keywords to set 
# the temperature, etc of the system. You must
# specify them in i-pi input file. 



# SYSTEM TARGET TEMPERATURE & PRESSURE
#pressure                0.001    ! katm (must be defined for NPT simulations)
temperature             300.00   ! Kelvin
timestep                0.0002   ! time step in ps
# CLOSING FILES
job time           1.0D+9
close time         1.0D+2

finish

'''.format(socket_line, alpha, *k)
    return output


def get_n_from_xyz(path):
    """Return the number of molecules from an XYZ file."""
    
    with open(path, 'r') as xyz:
        N = int(xyz.next().split()[0])  # Number of atoms
        assert N % 3 == 0
        return N / 3
