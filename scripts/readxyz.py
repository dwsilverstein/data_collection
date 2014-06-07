#! /usr/bin/env python

from __future__ import division, print_function
import sys
from collect import collect

def main():
    """\
    Script for reading XYZ files (*.xyz).  This is designed to work with the output
    of VMD solvation.

    CHANGELOG

    12-28-2013 DWS v1.3
    Modified for generating the environment file

    12-11-2013 DWS v1.2
    Modified for rearranging molecules based on proximity

    10-2-2013 DWS v1.1
    Modified to work for CO3, HCO3, H2CO3
    
    9-25-2013 DWS v1.0
    Initial build.  
    """

    from argparse import ArgumentParser, RawDescriptionHelpFormatter
    from textwrap import dedent
    parser = ArgumentParser(description=dedent(main.__doc__),
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('-m', '--molecules', help='Molecules in the simulation, including .  Read in '
                        'as "M1" "M2" ... list', nargs='+')
    parser.add_argument('-t', '--types', help='Atom types in the simulation.  Read in as '
                        '"# Atom" in a space separated list', nargs='+')
    parser.add_argument('-r', '--rearrange', help='Rearrange atoms in XYZ file based on proximity to solute.',
                        action='store_true', default=False)
    parser.add_argument('-w', '--within', help='Split the output based on the number of molecules within a '
                        'user defined distance to the solute.', required=False)
    parser.add_argument('-e', '--env', help='Generate an env file for the QM region in a CP2K QMMM '
                        'calculation.', action='store_true', required=False)
    parser.add_argument('-ns', '--nsolute', help='Number of atoms in the solute molecule.',
                        required=False)
    parser.add_argument('files', help='The files that you want converted.',
                        nargs='+')
    args = parser.parse_args()

    # Store command line arguments as convenient variables
    fh = args.files[0]
    molecules = args.molecules
    types = args.types
    rearrange = args.rearrange
    env = args.env
    nsolute = int(args.nsolute)
    if args.within:
        within = float(args.within)
        # Also store the squared distance because we use that later.
        within2 = within*within
    else:
        within = None

    # Molecule data, including number of atoms, central atom, charges
    moldata = { 'H2O'   : { 'Natoms' : 3, 'Catom' : 'O', 
                            'NBonds' : 2, 'NAngles' : 1, 
                            'Mult_cent' : False,
                            'Atoms' : ['O', 'H', 'H'],
                            'Connect' : ['H', 'H'], 'OW' : -0.82, 'HW' : 0.41 },
                'H3O'   : { 'Natoms' : 4, 'Catom' : 'O', 
                            'NBonds' : 3, 'NAngles' : 3,
                            'Mult_cent' : False,
                            'Atoms' : ['O', 'H', 'H', 'H'],
                            'Connect' : ['H', 'H', 'H'], 'OH' : -0.50,   'HH' : 0.50 },
                'CO2'   : { 'Natoms' : 3, 'Catom' : 'C', 
                            'NBonds' : 2, 'NAngles' : 1,
                            'Mult_cent' : False,
                            'Atoms' : ['C', 'O', 'O'],
                            'Connect' : ['O', 'O'], 'O'  : -0.3256, 'C'  : 0.6512 },
                'CO3'   : { 'Natoms' : 4, 'Catom' : 'C', 
                            'NBonds' : 3, 'NAngles' : 3,
                            'Mult_cent' : False,
                            'Atoms' : ['C', 'O', 'O', 'O'],
                            'Connect' : ['O', 'O', 'O'], 'O'  : -1.0000, 'C'  : 1.0000 },
                'HCO3'  : { 'Natoms' : 5, 'Catom' : 'C', 
                            'NBonds' : 4, 'NAngles' : 4,
                            'Mult_cent' : True,
                            'Atoms' : ['C', 'O', 'O', 'O', 'H'],
                            'Connect' : ['O', 'O', 'O'], 'O'  : -1.0000, 'C'  : 1.0000, 'H' : 1.0000 },
                'H2CO3' : { 'Natoms' : 6, 'Catom' : 'C', 
                            'NBonds' : 5, 'NAngles' : 5,
                            'Mult_cent' : True,
                            'Atoms' : ['C', 'O', 'O', 'O', 'H', 'H'],
                            'Connect' : ['O', 'O', 'O'], 'O'  : -1.0000, 'C'  : 1.0000, 'H' : 1.0000 },
              }

    # H2O parameters are from: jcp_124_024503
    # H3O parameters are from: jpcb_112_467
    # CO2 parameters are from: cjce_17_268
    # CO3 parameters are from: http://cinjweb.umdnj.edu/~kerrigje/data/par_all36_cgenff.prmi, I made up charges here
    # HCO3 parameters are made up
    # H2CO3 parameters are made up

    # Reformat the types
    if not rearrange and not env:
        atypes = {}
        for i in types:
            tmp = i.split()
            atypes[tmp[1]] = tmp[0] 

    # /\/\/\/\/\/\/\/\
    # Collect XYZ file 
    # /\/\/\/\/\/\/\/\

    xyz = collect(fh)

    if rearrange:
        # Number of atoms
        natoms = xyz['NATOMS']
        # Determine location of central atoms
        central = [[0,0]] # Index of central atom, number of atoms in molecule
        for i,atom in enumerate(xyz['ATOMS']):
            if i > 3: # Assumes the species is a carbonate
                if atom is 'O':
                     central.append([i,3]) # Assumes that subsequent molecules are H2O
        nmol = len(central)
        central[0][1] = central[1][0] - central[0][0]
        # Determine distance between the central atom in solute and solvent central atoms.
        # Note that this is actually the squared distance, to avoid a square root.
        distance = []
        for i,l in enumerate(central):
            if i == 0:
                continue
            else:
                j = l[0]
                d = ( ( xyz['COORDS'][j][0] - xyz['COORDS'][0][0] )**2
                    + ( xyz['COORDS'][j][1] - xyz['COORDS'][0][1] )**2
                    + ( xyz['COORDS'][j][2] - xyz['COORDS'][0][2] )**2 )
                distance.append(d)
        # Number of distances
        ndist = len(distance)
        # Connection of indices between the original coordinates and sorted coordinates
        indices = sort_coords(distance,nmol,ndist)
        # Sort coordinates
        satoms = [0.0 for i in range(natoms)]
        scoords = [0.0 for i in range(natoms)]
        for i in range(len(central)):
            m = indices[i][0]
            n = indices[i][1]
            j = central[n][0] 
            k = j + central[n][1] 
            p = central[m][0]
            for o in range(j,k):
                satoms[o] = xyz['ATOMS'][p]
                scoords[o] = xyz['COORDS'][p]
                p += 1
        # Sort the distances
        sdistance = sorted(distance)
        # Print values to file 
        if within:
            fmt = '{0} {1:14.8f} {2:14.8f} {3:14.8f}'
            tmp = str(within).split('.')
            # Determine the number of atoms in each region
            nwithin = 1
            noutside = 0
            for dist in sdistance:
                if dist <= within2:
                    nwithin += 1
                else:
                    noutside += 1
            natin = 6 + 3*nwithin
            natout = natoms - natin 
            # Check that the numbers of molecules in both regions is correct
            if nwithin + noutside != nmol:
                sys.exit('Number of atoms in each region is incorrect')
            # Name and open files
            fname1 = fh.split('.')[0] + '_leq' + tmp[0] + 'pt' + tmp[1] + '.xyz'
            fname2 = fh.split('.')[0] + '_gt' + tmp[0] + 'pt' + tmp[1] + '.xyz'
            f1 = open(fname1, 'w')
            f2 = open(fname2, 'w')
            # Write data to files
            print(natin, file=f1)
            print(natout, file=f2)
            print('', file=f1)
            print('', file=f2)
            for i in range(natin):
                print(fmt.format(satoms[i], scoords[i][0],
                                 scoords[i][1], scoords[i][2]), file=f1)
            for i in range(natin,natoms):
                print(fmt.format(satoms[i], scoords[i][0],
                                 scoords[i][1], scoords[i][2]), file=f2)
        else:
            fmt = '{0} {1:14.8f} {2:14.8f} {3:14.8f}'
            fname = 'newcoords.xyz'
            f = open(fname, 'w')
            print(natoms, file=f)
            print('', file=f)
            for i in range(natoms):
                print(fmt.format(satoms[i], scoords[i][0],
                                 scoords[i][1], scoords[i][2]), file=f)
    elif env:
        # /\/\/\/\/\/\/\/\/\/\/\/\/\
        # Write a CP2K QMMM env file 
        # /\/\/\/\/\/\/\/\/\/\/\/\/\
        
        # Conversion table
        table = { 'H' : 'Hqm',
                  'C' : 'Cqm',
                  'O' : 'Oqm',
                }
        
        # Convert QM atoms to CP2K types
        cp2katoms = []
        for atom in xyz['ATOMS']:
            cp2katoms.append(table[atom])

        # Write the env file with the QM atom identification
        envname = fh.split('.')[0] + '.env'
        envfile = open(envname, 'w')
        for i in range(len(cp2katoms)):
            print('&QM_KIND ' + cp2katoms[i], file=envfile)
            print('MM_INDEX ' + str(i+1), file=envfile)
            print('&END QM_KIND', file=envfile)
    else:
        # /\/\/\/\/\/\/\/\/\/\/
        # Output formatted data
        # /\/\/\/\/\/\/\/\/\/\/

        # Atoms block
        print('Atoms')
        print()

        # Locate the central atom.  This defines the start of a molecule 
        dcentral ={}

        # Break down the element list
        newstart = 0
        for m in molecules:
            # Store the central atom locations
            dcentral[m] = []    
            natoms = moldata[m]['Natoms']
            catom = moldata[m]['Catom']
            # Check if there are multiple "central" atoms
            mcent = moldata[m]['Mult_cent']
            for i in range(newstart,len(xyz['ATOMS']),natoms):
                c = xyz['ATOMS'][i]
                # For molecules with multiple centers, do this
                if mcent:
                    dcentral[m].append(i)
                    newstart = i + natoms
                    break
                else:
                    if c == catom:
                        keep = []
                        # Determine terminal atoms.
                        k = 0
                        for j in range(i+1,i+natoms):
                            if xyz['ATOMS'][j] == moldata[m]['Connect'][k]:
                                keep.append(True)
                            else:
                                keep.append(False)
                            k += 1
                        # Check if we have a molecule
                        if False in keep:
                            continue
                        else:
                            dcentral[m].append(i)
                            newstart = i + natoms
                            # Change atom types for easier sorting later
                            if m == 'H2O':
                                xyz['ATOMS'][i] = 'OW'
                                xyz['ATOMS'][i+1] = 'HW'
                                xyz['ATOMS'][i+2] = 'HW'
                            elif m == 'H3O':
                                xyz['ATOMS'][i] = 'OH'
                                xyz['ATOMS'][i+1] = 'HH'
                                xyz['ATOMS'][i+2] = 'HH'
                                xyz['ATOMS'][i+3] = 'HH'

        # Start printing data
        fmt = '{0:>4} {1:>4} {2:>2} {3:11.4e} {4:11.4e} {5:11.4e} {6:11.4e}'
        counter = 1 # Counter for atoms
        mnum = 1    # Counter for molecules

        morder = [] # Order of molecules
        ncentral = {} # New index of central atom.
        for m in molecules:
            morder.append(m)
            ncentral[m] = []
            for i in dcentral[m]:
                ln = []
                # Add counter for atoms
                ln.append(str(counter))
                # Append to the new counter
                ncentral[m].append(counter)
                # Add counter for molecules
                ln.append(str(mnum))
                # Add atomtype
                ln.append(atypes[xyz['ATOMS'][i]])
                # Add in atomic charge
                ln.append(moldata[m][xyz['ATOMS'][i]])
                # Add in coordinates
                ln.append(xyz['COORDS'][i][0])
                ln.append(xyz['COORDS'][i][1])
                ln.append(xyz['COORDS'][i][2])
                # Print the central atom
                print(fmt.format(ln[0], ln[1], ln[2], ln[3], ln[4], ln[5], ln[6]))
                # Repeat the process for the terminal atoms
                for j in range(i+1,i+moldata[m]['Natoms']):
                    counter += 1
                    ln = []
                    ln.append(str(counter))
                    ln.append(str(mnum))
                    ln.append(atypes[xyz['ATOMS'][j]])
                    ln.append(moldata[m][xyz['ATOMS'][j]])
                    ln.append(xyz['COORDS'][j][0])
                    ln.append(xyz['COORDS'][j][1])
                    ln.append(xyz['COORDS'][j][2])
                    print(fmt.format(ln[0], ln[1], ln[2], ln[3], ln[4], ln[5], ln[6]))
                # Change molecule number
                mnum += 1
                # Change atom number
                counter += 1

        # Bonds block 
        # Note that I assume bond type 1 belongs to water.
        print()
        print('Bonds')
        print()

        fmt = '{0:>4} {1:>4} {2:>4} {3:>4}'

        bnum = 1
        for m in morder:
            if m == 'H2O':
                btype = 1
            else:
                btype = 2
            for i in ncentral[m]:
                ln = []
                ln.append(bnum)
                ln.append(btype)
                ln.append(i)
                for j in range(i+1,i+moldata[m]['Natoms']):
                    ln.append(j)
                    print(fmt.format(ln[0], ln[1], ln[2], ln[3]))
                    ln.pop()
                    bnum += 1
                    ln[0] = bnum

        # Angles block
        # Again, I assume that angle type 1 belongs to water
        print()
        print('Angles')
        print()

        # Keep in mind that the order of print out is terminal, central, terminal
        fmt = '{0:>4} {1:>4} {2:>4} {3:>4} {4:>4}'

        anum = 1
        # This part will probably need to be revised for molecules with more than 3 atoms
        for m in morder:
            if m == 'H2O':
                a = 1
            else:
                a = 2
            for i in ncentral[m]:
                ln = [0,0,0,0,0]
                ln[0] = anum
                ln[1] = a
                ln[3] = i
                ln[2] = i+1
                ln[4] = i+2
                #t = 2
                #for j in range(i+1,i+moldata[m]['Natoms']):
                #    ln[t] = j
                #    t += 2
                # Not the smartest way to make this work, but it will do.
                nangles = moldata[m]['NAngles']
                for j in range(nangles):
                    print(fmt.format(ln[0], ln[1], ln[2], ln[3], ln[4]))
                    anum += 1

def sort_coords(dist,nmol,ndist):
    '''Sort coordinates based on distance.'''

    sortval = [[i,0] for i in range(nmol)]

    if ndist == nmol - 1:
        # Life is easy
        sdist = sorted(dist)
        for i in range(1,nmol):
            j = i - 1
            ind = sdist.index(dist[j]) + 1
            sortval[i][1] = ind 
    #else:
    #    # Life is hard

    return sortval

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
