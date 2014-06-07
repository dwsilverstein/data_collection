#! /usr/bin/env python

from __future__ import division, print_function
import sys
import numpy as np
from numpy import append

def main():
    """\
    Script for generating a PSF file using a LAMMPS datafile. 
    """

    from argparse import ArgumentParser, RawDescriptionHelpFormatter
    from textwrap import dedent
    parser = ArgumentParser(description=dedent(main.__doc__),
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('files', help='The input files that you want converted.',
                        nargs='+')
    parser.add_argument('-ns', '--numsolute', help='Number of atoms in solute molecule',
                        required=True)
    parser.add_argument('-a', '--atomtypes', help='Atom type conversion',
                        nargs='+')
    args = parser.parse_args()

    # Store command line arguments as convenient variables
    fh = args.files[0]
    nsolute = int(args.numsolute)
    convert = {}
    for item in args.atomtypes:
        item = item.split()
        convert[item[0]] = item[1]

    # Collect the file into memory for retention
    with open(fh) as fl:
        f = [line.rstrip() for line in fl]

    # Find specific lines in the file
    lines = ( 'atoms', 'bonds', 'angles', 'dihedrals', 'impropers',
              'Atoms', 'Bonds', 'Angles', 'Dihedrals', 'Impropers',
            )

    table = { 'NATOMS'        : 0,
              'NBONDS'        : 0,
              'NANGLES'       : 0,
              'NDIHEDRALS'    : 0,
              'NIMPROPERS'    : 0,
              'ATOMSTART'     : 0,
              'BONDSTART'     : 0,
              'ANGLESTART'    : 0,
              'DIHEDRALSTART' : 0,
              'IMPROPERSTART' : 0,
              'MASSLIST'      : 0,
            }

    for i, ln in enumerate(f):
        tmp = ln.split()
        if len(tmp) == 2 and tmp[1] in lines:
            key = 'N' + tmp[1].upper()
            table[key] = int(tmp[0]) 
        elif len(tmp) == 1 and tmp[0] in lines:
            key = tmp[0].upper() + 'TART'
            table[key] = i + 2
        elif 'Masses' in tmp:
            table['MASSLIST'] = i + 2
    
    # Get the list of masses and atom indices
    masses = {}
    s = table['MASSLIST']
    e = [ s + i for i, x in enumerate(f[s:]) if x == ''][0]
    for item in f[s:e]:
        tmp = item.split()
        masses[tmp[0]] = float(tmp[1])

    # Blocks in a PSF file
    # 1. PSF
    # 2. !NTITLE
    # 3. !NATOM
    # 4. !NBONDS
    # 5. !NTHETA
    # 6. !NPHI
    # 7. !NIMPHI
    psfhead = [ 'PSF', 
                '',
                '         0 !NTITLE',
                '',
              ] 
    
    # Atoms part
    natoms = table['NATOMS']
    fmt = '{0:>8} {1}'
    atomshead = fmt.format( str(natoms), '!NATOM')
    s = table['ATOMSTART']
    e = s + natoms
    count = 1
    atomblock = []
    fmt = '{0:>8} {1:<4} {2:<4} {3:<4} {4:<4} {5:<4} {6:>10.6f}       {7:>7.4f} {8:>11}' 
    for i in range(s,e):
        tmp = f[i].split()
        mol = int(tmp[1])
        atom = convert[tmp[2]] 
        charge = float(tmp[3])
        # Molecule name
        if count <= nsolute:
            mname = 'SOL' # Solute molecule
            aname = atom
            atype = atom
        else:
            mname = 'SPC' # Water molecule
            if atom == 'O' or atom == 'OW':
                aname = 'O'
                atype = 'OW'
            else:
                aname = 'H'
                atype = 'HW'
        tmp2 = fmt.format( str(count), mname, str(mol),
                           mname, aname, atype,
                           charge, masses[tmp[2]], '0' )
        atomblock.append(tmp2) 
        count += 1

    # Bonds part
    nbonds = table['NBONDS']
    fmt = '{0:>8} {1}'
    bondshead = fmt.format( str(nbonds), '!NBONDS')
    s = table['BONDSTART']
    e = s + nbonds
    bondblock = []
    lenline = 64 # Length of a line in the bond block
    ln = ''
    fmt = '{0:>4}'
    #remainder = nbonds / 4.0
    for i in range(s,e):
        tmp = f[i].split()
        one = fmt.format(str(tmp[2])) 
        two = fmt.format(str(tmp[3])) 
        ln += '    ' + one + '    ' + two 
        if len(ln) == lenline:
            bondblock.append(ln)
            ln = ''
        elif i == e-1:
            # This part handles the remainder of bonds
            # if we don't have a complete list for the final line.
            bondblock.append(ln)
    
    # Angles part
    nangles = table['NANGLES']
    fmt = '{0:>8} {1}'
    angleshead = fmt.format( str(nangles), '!NTHETA')
    s = table['ANGLESTART']
    e = s + nangles
    angleblock = []
    lenline = 72 # Length of a line in the angle block
    ln = ''
    fmt = '{0:>4}'
    for i in range(s,e):
        tmp = f[i].split()
        one = fmt.format(str(tmp[2]))
        two = fmt.format(str(tmp[3]))
        three = fmt.format(str(tmp[4]))
        ln += '    ' + one + '    ' + two + '    ' + three
        if len(ln) == lenline:
            angleblock.append(ln)
            ln = ''
        elif i == e-1:
            # This part handles the remainder of angles
            # if we don't have a complete list for the final line.
            angleblock.append(ln)

    # Dihedrals part
    ndihedrals = table['NDIHEDRALS']
    fmt = '{0:>8} {1}'
    dihedralshead = fmt.format( str(ndihedrals), '!NPHI')
    s = table['DIHEDRALSTART']
    e = s + ndihedrals
    dihedralblock = []
    lenline = 64 # Length of a line in the angle block
    ln = ''
    fmt = '{0:>4}'
    for i in range(s,e):
        tmp = f[i].split()
        one = fmt.format(str(tmp[2]))
        two = fmt.format(str(tmp[3]))
        three = fmt.format(str(tmp[4]))
        four = fmt.format(str(tmp[5]))
        ln += '    ' + one + '    ' + two + '    ' + three + '    ' + four
        if len(ln) == lenline:
            dihedralblock.append(ln)
            ln = ''
        elif i == e-1:
            # This part handles the remainder of angles
            # if we don't have a complete list for the final line.
            dihedralblock.append(ln)

    # Impropers part
    nimpropers = table['NIMPROPERS']
    fmt = '{0:>8} {1}'
    impropershead = fmt.format( str(nimpropers), '!NIMPHI')
    s = table['IMPROPERSTART']
    e = s + nimpropers
    improperblock = []
    lenline = 64 # Length of a line in the angle block
    ln = ''
    fmt = '{0:>4}'
    for i in range(s,e):
        tmp = f[i].split()
        one = fmt.format(str(tmp[2]))
        two = fmt.format(str(tmp[3]))
        three = fmt.format(str(tmp[4]))
        four = fmt.format(str(tmp[5]))
        ln += '    ' + one + '    ' + two + '    ' + three + '    ' + four
        if len(ln) == lenline:
            improperblock.append(ln)
            ln = ''
        elif i == e-1:
            # This part handles the remainder of angles
            # if we don't have a complete list for the final line.
            improperblock.append(ln)

    # Write the PSF file
    psfname = fh.split('.')[0] + '.psf'
    psf = open(psfname, 'w')

    for i in psfhead:
        print(i, file=psf)
    print(atomshead, file=psf)
    for i in atomblock:
        print(i, file=psf)
    print('', file=psf)
    print(bondshead, file=psf)
    for i in bondblock:
        print(i, file=psf)
    print('', file=psf)
    print(angleshead, file=psf)
    for i in angleblock:
        print(i, file=psf)
    print('', file=psf)
    print(dihedralshead, file=psf)
    for i in dihedralblock:
        print(i, file=psf)
    print('', file=psf)
    print(impropershead, file=psf)
    for i in improperblock:
        print(i, file=psf)
    print('', file=psf)

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
