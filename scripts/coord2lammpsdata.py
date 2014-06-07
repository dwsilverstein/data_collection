#! /usr/bin/env python

from __future__ import division, print_function
import sys
from collect import collect

def main():
    """\
    Script for putting coordinates from an XYZ file into a template LAMMPS data file.
    This works by taking the coordinates from the XYZ file and transplanting them
    into the already made LAMMPS data file.

    """

    from argparse import ArgumentParser, RawDescriptionHelpFormatter
    from textwrap import dedent
    parser = ArgumentParser(description=dedent(main.__doc__),
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('-t', '--template', help='Template LAMMPS data file.',
                        required=True)
    parser.add_argument('files', help='The XYZ file(s).', nargs='+')
    args = parser.parse_args()

    # Store command line arguments as convenient variables
    fh = args.files
    template = args.template

    # Read template LAMMPS data file into memory 
    with open(template) as fl:
        tf = tuple([line.rstrip() for line in fl])

    natoms = int(tf[2].split()[0]) # Number of atoms
    nlines = len(tf)               # Number of lines

    # Determine where to insert coordinates in the template file
    for i,ln in enumerate(tf):
        if 'Atoms' in ln:
            s = i + 2
            e = s + natoms
            break

    fmtln = '{0:>4} {1:>4} {2:>2} {3:>12} {4:11.4e} {5:11.4e} {6:11.4e}'

    # Add the coordinates to the template file for each XYZ file
    for f in fh:
        xyz = collect(f)
        coords = xyz['COORDS']
        fname = f.split('.')[0] + '.data'
        fl = open(fname,'w')

        # Print header
        for i in range(s):
            print(tf[i], file=fl)
        # Print atomblock
        j = 0
        for i in range(s,e):
            tmp = tf[i].split()
            print( fmtln.format( tmp[0], tmp[1], tmp[2], tmp[3],
                                 coords[j][0],
                                 coords[j][1],
                                 coords[j][2] ), file=fl )
            j += 1
            print(j)
        # Print end
        for i in range(e,nlines):
            print(tf[i], file=fl)
        
        fl.close()

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
