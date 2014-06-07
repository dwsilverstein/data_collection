#! /usr/bin/env python

from __future__ import division, print_function
import sys, os
import numpy as np
from numpy import array, append, average, std
from collect import collect

def main():
    """\
    Script for reading a LAMMPS files.  These can be output (*.out), log 
    (log.lammps, for instance), or data files.
    """

    from argparse import ArgumentParser, RawDescriptionHelpFormatter
    from textwrap import dedent
    parser = ArgumentParser(description=dedent(main.__doc__),
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('-o', '--output', help='Output the data to a table.',
                        action='store_true', default=False)    
    parser.add_argument('-p', '--plot', help='Plot the data.',
                        action='store_true', default=False)    
    parser.add_argument('-a', '--atomtypes', help='The atom type numbers you want converted.'
                        ' These should be written in "# element" format.', nargs='+')
    parser.add_argument('-k', '--cp2k', help='Create coordinate blocks for CP2K.',
                        action='store_true', default=False)
    parser.add_argument('files', help='The files that you want converted.'
                        , nargs='+')
    args = parser.parse_args()

    # Store command line arguments as convenient variables
    fh = args.files[0]
    output = args.output
    plot = args.plot
    atable = args.atomtypes
    cp2k = args.cp2k
    
    if output is False and plot is False:
        msg = "Please request output or plotting.  Use -h to see what is available"
        error = '\n' + len(msg)*'%' + '\n' + msg + '\n' + len(msg)*'%' + '\n'
        sys.exit(error)

    # Convert the table to a dictionary
    if cp2k:
        convert = {}
        for item in atable:
            item = item.split()
            convert[item[0]] = item[1]

    # Collect all data into memory for data retention
    thermodata = collect(fh)

    # Strings to search for
    table = {'Define variables' : 'VARIABLES',
             'Define units' : 'UNITS',
             'Define LJ potential input' : 'LJPOTENTIAL',
             'Define Potential Types' : 'POTENTIALTYPES',
             'Define Electrostatics' : 'ELECTROSTATICS',
             'Read Data for System' : 'READDATA',
             'Coefficients for Bond Potentials' : 'BONDS',
             'Coefficients for Angle Potentials' : 'ANGLES',
             'Coefficients for LJ Potential' : 'LJCOEFF',
             'Parameters for Building Pairwise Neighbor Lists' : 'NEIGHBORS',
             'Define Timestep Size for MD' : 'TIMESTEP',
             'Define the Thermodynamics Output' : 'THERMOUT',
             'Define Operations Applied to the System' : 'FIXES',
             'Define Extra Output' : 'DUMPS',
             'Restart File Information' : 'RESTART',
             'Run Dynamics' : 'RUNDYN',
             'End of Job' : 'END'}

    # Determine which units were used.  This is important for correct output.
    units_table = {'REAL'       : { 'Mass'              : r'g/mol',
                                    'Distance'          : r'\AA',
                                    'Time'              : r'fs',
                                    'Energy'            : r'kcal/mol',
                                    'Velocity'          : r'\AA/fs',
                                    'Force'             : r'kcal/(mol*\AA)',
                                    'Torque'            : r'kcal/mol',
                                    'Temperature'       : r'K',
                                    'Pressure'          : r'atm',
                                    'Dynamic Viscosity' : r'Poise',
                                    'Charge'            : r'e',
                                    'Dipole'            : r'e*\AA',
                                    'Electric Field'    : r'V/\AA',
                                    'Density'           : r'g/cm$^3$' },
                   'METAL'      : { 'Mass'              : r'g/mol',
                                    'Distance'          : r'\AA',
                                    'Time'              : r'ps',
                                    'Energy'            : r'eV',
                                    'Velocity'          : r'\AA/ps',
                                    'Force'             : r'eV/\AA',
                                    'Torque'            : r'eV',
                                    'Temperature'       : r'K',
                                    'Pressure'          : r'bars',
                                    'Dynamic Viscosity' : r'Poise', 
                                    'Charge'            : r'e',
                                    'Dipole'            : r'e*\AA', 
                                    'Electric Field'    : r'V/\AA', 
                                    'Density'           : r'g/cm$^3$' },
                   'SI'         : { 'Mass'              : r'kg',
                                    'Distance'          : r'm',
                                    'Time'              : r's',
                                    'Energy'            : r'J',
                                    'Velocity'          : r'm/s',
                                    'Force'             : r'N',
                                    'Torque'            : r'N*m',
                                    'Temperature'       : r'K',
                                    'Pressure'          : r'Pa',
                                    'Dynamic Viscosity' : r'Pa*s', 
                                    'Charge'            : r'C',
                                    'Dipole'            : r'C*m', 
                                    'Electric Field'    : r'V/m', 
                                    'Density'           : r'kg/cm$^3$' },
                   'CGS'        : { 'Mass'              : r'g',
                                    'Distance'          : r'cm',
                                    'Time'              : r's',
                                    'Energy'            : r'ergs',
                                    'Velocity'          : r'cm/s',
                                    'Force'             : r'dynes',
                                    'Torque'            : r'dyne*cm',
                                    'Temperature'       : r'K',
                                    'Pressure'          : r'dyne/cm^2',
                                    'Dynamic Viscosity' : r'Poise', 
                                    'Charge'            : r'esu',
                                    'Dipole'            : r'esu*cm', 
                                    'Electric Field'    : r'dyne/esu', 
                                    'Density'           : r'g/cm$^3$' },
                   'ELECTRON'   : { 'Mass'              : r'amu',
                                    'Distance'          : r'Bohr',
                                    'Time'              : r'fs',
                                    'Energy'            : r'Hartrees',
                                    'Velocity'          : r'Bohr/(atomic time unit)',
                                    'Force'             : r'Hartree/Bohr',
                                    'Temperature'       : r'K',
                                    'Pressure'          : r'Pa',
                                    'Charge'            : r'e',
                                    'Dipole'            : r'Debye', 
                                    'Electric Field'    : r'V/cm', },
                  }

    # Table of quantities for performing averages
    comp_table = ( 'Temp', 'TotEng', 'PotEng', )

    # Output some information to screen
    if output:
        if thermodata['FILETYPE'] == 'LAMMPS Logfile':
            # Define the units
            units = thermodata['UNITS']

            # thermotable converts between the keys in thermodata to column headings.
            thermotable = { 'Step'   : 'Time (' + units['Time'] + ')', 
                            'Temp'   : 'Temperature (' + units['Temperature'] + ')', 
                            'PotEng' : 'Potential Energy (' + units['Energy'] + ')', 
                            'TotEng' : 'Total Energy (' + units['Energy'] + ')', }

            # Generate the column heading format and data format
            fmt = '|{0:20}|'
            datafmt = '|{0:20}|'
            count = 1
            for i in range(len(thermodata.keys())):
                if thermodata.keys()[i] == 'Step':
                    count = count
                elif thermodata.keys()[i] in comp_table:
                    fmt += '{' + str(count) + ':^30}|'
                    datafmt += '{' + str(count) + ':^30.3f}|'
                    count += 1

            # Generate the column headings
            head    = ['']
            avg     = ['Average']
            stdev   = ['Standard Deviation']
            for key in thermodata.keys():
                if key in comp_table:
                    head.append(thermotable[key])
                    avg.append(thermodata[key][2])
                    stdev.append(thermodata[key][3])

            print()
            print('='*len(fmt.format(*head)))        
            print(fmt.format(*head))
            print('='*len(fmt.format(*head)))        
            print(datafmt.format(*avg))
            print(datafmt.format(*stdev))
            print('='*len(fmt.format(*head)))        
            print()
        elif thermodata['FILETYPE'] == 'LAMMPS Datafile':
            coords = []
            if 'Velocities' in thermodata.keys():
                velocities = []
            natoms = thermodata['Num_atoms'] 
            atoms = []
            for i in range(natoms):
                ix = str(i+1)
                atoms.append(convert[str(thermodata['Coords'][ix][0])])
                coords.append([ thermodata['Coords'][ix][1],
                                thermodata['Coords'][ix][2], 
                                thermodata['Coords'][ix][3] ])
                if 'Velocities' in thermodata.keys():
                    velocities.append([ thermodata['Velocities'][ix][0],
                                        thermodata['Velocities'][ix][1], 
                                        thermodata['Velocities'][ix][2] ])

            # Print coordinates
            print('&COORD')
            lnstyle = ' {0:2} {1:13.8f} {2:13.8f} {3:13.8f}'
            for i in range(natoms):
                print(lnstyle.format(atoms[i],
                                     coords[i][0], coords[i][1], coords[i][2]))
            print('&END COORD')

            # Print velocities
            lenconv = 0.5291772          # Angstrom / Bohr
            timeconv = 2.418884326505e-2 # fs / a.t.u
            conv = timeconv / lenconv
            if 'Velocities' in thermodata.keys():
                # Velocity units are by default atomic units.  We need
                # to convert from values in LAMMPS as a result for CP2K.
                print('&VELOCITY')
                lnstyle = ' {0:13.8f} {1:13.8f} {2:13.8f}'
                for i in range(natoms):
                    print(lnstyle.format(conv*velocities[i][0], 
                                         conv*velocities[i][1], 
                                         conv*velocities[i][2]))
                print('&END VELOCITY')

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
