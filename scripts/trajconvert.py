#! /usr/bin/env python

from __future__ import division, print_function
import sys

def main():
    """\
    Script for processing a LAMMPS trajectory file (.lammpstrj).  This can be
    used in two ways:

    1. To convert atom type numbers to element names for generating movies in VMD.

    2. To generate an MS-EVB topology file for RAPTOR using the atom information
    in the trajectory file.
    """

    from argparse import ArgumentParser, RawDescriptionHelpFormatter
    from textwrap import dedent
    parser = ArgumentParser(description=dedent(main.__doc__),
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s 1.1')
    parser.add_argument('-a', '--atomtypes', help='The atom type numbers you want converted.'
                        ' These should be written in "# element" format.', nargs='+')
    parser.add_argument('-n', '--newfile', help='New file name.', nargs='+')
    parser.add_argument('-c', '--convert', help='Convert atom type numbers to element names.',
                        action='store_true', default=False)
    parser.add_argument('-t', '--topology', help='Create an MS-EVB topology file from the final '
                        'coordinates of the trajectory.', action='store_true', default=False)
    parser.add_argument('-e', '--evbin', help='The MS-EVB input file containing atom type information.') 
    parser.add_argument('-d', '--datafile', help='The data file containing atom ID number information.') 
    parser.add_argument('files', help='The files that you want converted.'
                        , nargs='+')
    parser.add_argument('-k', '--cp2k', help='Create coordinate blocks for CP2K.',
                        action='store_true', default=False)
    parser.add_argument('-s', '--snapshot', help='Chooses the snapshot of the system for CP2K output.',
                        default=-1)
    args = parser.parse_args()

    # Store command line arguments as convenient variables
    fh = args.files[0]
    table = args.atomtypes
    conv = args.convert
    topology = args.topology
    evbin = args.evbin
    datafile = args.datafile
    cp2k = args.cp2k
    snapshot = int(args.snapshot)

    # Check that the user specified either convert or topology.  Die here if not.
    if conv is False and topology is False and cp2k is False:
        error = 'Must specify the -c, -t, or -k flag.  Use -h flag for information.'
        sys.exit('\n' + '%'*len(error) + '\n' + error + '\n' + '%'*len(error) + '\n')

    # Check if the user is trying to convert the file.  Exit here if -c is specified
    # without giving the atom type numbers.
    if conv is True and table is None:
        error = 'Must use -a flag when converting a trajectory file.  Use -h flag for information.'
        sys.exit('\n' + '%'*len(error) + '\n' + error + '\n' + '%'*len(error) + '\n')

    # Check if the user is trying to generate a topology file.  Exit here if
    # -t is specified without giving the evbin, par4amber and data files.
    if topology is True and (evbin is None or datafile is None):
        error = 'Must use -e and -d flags when generating a topology file.  Use -h flag for information.'
        sys.exit('\n' + '%'*len(error) + '\n' + error + '\n' + '%'*len(error) + '\n')
    
    # Convert the table to a dictionary
    if conv or cp2k:
        convert = {}
        for item in table:
            item = item.split()
            convert[item[0]] = item[1]

    # Data blocks
    data = ['ITEM: TIMESTEP',
            'ITEM: NUMBER OF ATOMS',
            'ITEM: BOX BOUNDS pp pp pp',
            'ITEM: ATOMS id type x y z fx fy fz',
            'ITEM: ATOMS id type x y z ix iy iz',
            'ITEM: ATOMS id mol type x y z ix iy iz',
            'ITEM: ATOMS id type q x y z', 
            'ITEM: ATOMS id mol type q x y z']

    # Collect all data into memory for data retention
    with open(fh) as fl:
        f = tuple([line.rstrip() for line in fl])

    # Locate data blocks
    # Need to find:
    # 1. Timestep value
    # 2. Number of atoms
    # 3. Box bounds
    # 4. Coordinates/forces for each timestep
    timestep = []
    natoms = 0
    box_bounds = []
    coord_and_force = []
    coord_and_index = []
    coord_and_charge = []
    coord_molid_and_index = []
    coord_and_molid = []
    for i,line in enumerate(f):
        # Timestep value
        if 'TIMESTEP' in line:
            timestep.append(int(f[i+1]))
        # Number of atoms
        if natoms == 0:
            if 'NUMBER OF ATOMS' in line:
                natoms = int(f[i+1])
        # Box bounds
        #if box_bounds == []:
        if ('BOX BOUNDS pp pp pp' in line or
            'BOX BOUNDS ff ff ff' in line or
            'BOX BOUNDS fm fm fm' in line or
            'BOX BOUNDS mm mm mm' in line):
            tmp1 = f[i+1].split()
            tmp2 = f[i+2].split()
            tmp3 = f[i+3].split()
            tmp4 = []
            box_bounds.append([ [float(tmp1[0]),float(tmp1[1])],
                                [float(tmp2[0]),float(tmp2[1])],
                                [float(tmp3[0]),float(tmp3[1])] ])
        # Coordinates and forces
        if ('ITEM: ATOMS id type x y z fx fy fz' in line or
            'ITEM: ATOMS id type xu yu zu fx fy fz' in line):
            # Temporary list
            tmp = []
            for ix in range(i+1,i+1+natoms):
                ln = f[ix] # Grab the current line
                ln = ln.split() # Split on whitespace
                if conv or cp2k:
                    ln[1] = convert[ln[1]] # Use the conversion table to convert the label to element
                else:
                    ln[1] = int(ln[1])
                ln[2] = float(ln[2])
                ln[3] = float(ln[3])
                ln[4] = float(ln[4])
                ln[5] = float(ln[5])
                ln[6] = float(ln[6])
                ln[7] = float(ln[7])
                tmp.append(ln) # Add the line to tmp
            coord_and_force.append(tmp) # Make coord_and_force a list of lists.
        # Coordinates and box indices
        if 'ITEM: ATOMS id type x y z ix iy iz' in line:
            tmp = []
            for ix in range(i+1,i+1+natoms):
                ln = f[ix] # Grab the current line
                ln = ln.split() # Split on whitespace
                if conv or cp2k:
                    ln[1] = convert[ln[1]] # Use the conversion table to convert the label to element
                else:
                    ln[1] = int(ln[1])
                ln[2] = float(ln[2])
                ln[3] = float(ln[3])
                ln[4] = float(ln[4])
                ln[5] = int(ln[5])
                ln[6] = int(ln[6])
                ln[7] = int(ln[7])
                tmp.append(ln) # Add the line to tmp
            coord_and_index.append(tmp) # Make coord_and_force a list of lists.
        # Coordinates, charges
        if 'ITEM: ATOMS id type q x y z' in line:
            tmp = []
            for ix in range(i+1,i+1+natoms):
                ln = f[ix] # Grab the current line
                ln = ln.split() # Split on whitespace
                if conv or cp2k:
                    ln[1] = convert[ln[1]] # Use the conversion table to convert the label to element
                else:
                    ln[1] = int(ln[1])
                ln[2] = float(ln[2])
                ln[3] = float(ln[3])
                ln[4] = float(ln[4])
                ln[5] = float(ln[5])
                tmp.append(ln) # Add the line to tmp
            coord_and_charge.append(tmp) # Make coord_and_force a list of lists.
        # Coordinates, mol_id and box indices
        if 'ITEM: ATOMS id mol type x y z ix iy iz' in line:
            tmp = []
            for ix in range(i+1,i+1+natoms):
                ln = f[ix] # Grab the current line
                ln = ln.split() # Split on whitespace
                if conv or cp2k:
                    ln[2] = convert[ln[2]] # Use the conversion table to convert the label to element
                else:
                    ln[2] = int(ln[2])
                ln[3] = float(ln[3])
                ln[4] = float(ln[4])
                ln[5] = float(ln[5])
                ln[6] = int(ln[6])
                ln[7] = int(ln[7])
                ln[8] = int(ln[8])
                tmp.append(ln) # Add the line to tmp
            coord_molid_and_index.append(tmp) # Make coord_molid_and_force a list of lists.
        # Coordinates and mol_id
        if 'ITEM: ATOMS id mol type q x y z' in line:
            tmp = []
            for ix in range(i+1,i+1+natoms):
                ln = f[ix] # Grab the current line
                ln = ln.split() # Split on whitespace
                if conv or cp2k:
                    ln[2] = convert[ln[2]] # Use the conversion table to convert the label to element
                else:
                    ln[2] = int(ln[2])
                ln[3] = float(ln[3])
                ln[4] = float(ln[4])
                ln[5] = float(ln[5])
                ln[6] = float(ln[6])
                tmp.append(ln) # Add the line to tmp
            coord_and_molid.append(tmp) # Make coord_molid_and_force a list of lists.

    # Close the trajectory file and empty the tuple
    f = () 

    # Convert atom type numbers to element names.
    if conv:
        # Number of timesteps
        ntstep = len(timestep)

        # Print data to a new file.
        nfile = []
        if args.newfile is not None:
            for item in args.newfile:
                nfile.append(args.newfile)
        else:
            nfile.append(fh.split('.')[0] + '_new.lammpstrj')

        lstyle1 = '{0:>9.5f} {1:>9.5f}'
        lstyle2 = '{0:>4} {1:>4} {2:>9.5f} {3:>9.5f} {4:>9.5f} {5:>9.5f} {6:>9.5f} {7:>9.5f}'
        lstyle3 = '{0:>4} {1:>4} {2:>9.5f} {3:>9.5f} {4:>9.5f} {5:>4} {6:>4} {7:>4}'
        lstyle4 = '{0:>4} {1:>4} {2:>2} {3:>9.5f} {4:>9.5f} {5:>9.5f} {6:>4} {7:>4} {8:>4}'
        lstyle5 = '{0:>4} {1:>4} {2:>5.2f} {3:>9.5f} {4:>9.5f} {5:>9.5f}'
        lstyle6 = '{0:>4} {1:>4} {2:>2} {3:>5.2f} {4:>9.5f} {5:>9.5f} {6:>9.5f}'
        for item in nfile:
            fl = open(item, 'w')
            for ix in range(ntstep):
                print(data[0], file=fl)
                print(timestep[ix], file=fl)
                print(data[1], file=fl)
                print(natoms, file=fl)
                print(data[2], file=fl)
                print(lstyle1.format(box_bounds[ix][0][0], box_bounds[ix][0][1]), file=fl)
                print(lstyle1.format(box_bounds[ix][1][0], box_bounds[ix][1][1]), file=fl)
                print(lstyle1.format(box_bounds[ix][2][0], box_bounds[ix][2][1]), file=fl)
                if coord_and_force != []:
                    print(data[3], file=fl)
                    for ln in coord_and_force[ix]:
                        print(lstyle2.format(ln[0], ln[1], ln[2], ln[3], ln[4], ln[5], ln[6], ln[7]), file=fl)
                elif coord_and_index != []:
                    print(data[4], file=fl)
                    for ln in coord_and_index[ix]:
                        print(lstyle3.format(ln[0], ln[1], ln[2], ln[3], ln[4], ln[5], ln[6], ln[7]), file=fl)
                elif coord_and_charge != []:
                    print(data[6], file=fl)
                    for ln in coord_and_charge[ix]:
                        print(lstyle5.format(ln[0], ln[1], ln[2], ln[3], ln[4], ln[5]), file=fl)
                elif coord_molid_and_index != []:
                    print(data[5], file=fl)
                    for ln in coord_molid_and_index[ix]:
                        print(lstyle4.format(ln[0], ln[1], ln[2], ln[3], ln[4], ln[5], ln[6], ln[7], ln[8]), file=fl)
                elif coord_and_molid != []:
                    print(data[7], file=fl)
                    for ln in coord_and_molid[ix]:
                        print(lstyle6.format(ln[0], ln[1], ln[2], ln[3], ln[4], ln[5], ln[6]), file=fl)
    # Generate a topology file for RAPTOR.
    elif topology:
        # Trajectory data
        traj_data = coord_and_force[-1]

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
        # Read the input file for the MS-EVB code
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
        # This tells us the atom types in RAPTOR related to the LAMMPS atom types
        with open(evbin) as data:
            evb = tuple([line.rstrip() for line in data])

        # Read the LAMMPS types and molecule numbers from the MS-EVB input
        lammps_type = {}
        reaction_mol_num = {}
        for i,line in enumerate(evb):
            if 'MODEL' in line:
                s1 = i + 3
            if 'Reaction' in line:
                m1 = i + 1
            if ': atom type' in line:
                s2 = i + 2
            if ': bond type' in line:
                e2 = i - 1
        
        # Store the reaction type(s)
        for data in evb[m1:]:
            if '#define' in data:
                tmp = data.split()
                reaction_mol_num[tmp[1]] = {}
            # Stop the loop if a blank line is encountered
            if data == '':
                break
        
        # Store the molecule names and numbers belonging to a reaction type
        # I assume that the naming convention of reactions has correspondence
        # with the naming of the molecular model
        for i,data in enumerate(evb[s1:m1]):
            # This part is pretty convoluted.  At least it works.
            if '#define' in data:
                count = s1 + i + 1
                tmp = data.split()
                # Here we loop through the reaction types we stored above
                for reaction in reaction_mol_num.keys():
                    if tmp[1] in reaction.split('_'):
                        reaction_mol_num[reaction][tmp[1]] = {}
                        nloop = 0
                        # Here we look through lines greater and stash the molecule
                        # number information. 
                        for line in evb[count:]:
                            mol = line.split()
                            if '#define' in mol:
                                break
                            elif mol == []:
                                break
                            else:
                                reaction_mol_num[reaction][tmp[1]][mol[1]] = int(mol[2])

        # Store the LAMMPS types
        for data in evb[s2:e2]:
            if '#define' in data:
                tmp = data.split()
                lammps_type[tmp[1]] = int(tmp[2])

        # Clear out the tuple (free memory)
        evb = ()

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
        # Read the data file from the LAMMPS output
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
        # This tells us the atom topology information
        with open(datafile) as data:
            df = tuple([line.rstrip() for line in data])

        # Collect the 'Atoms' information from the data file.
        for i,line in enumerate(df):
            if 'Atoms' in line:
                s = i + 2
            elif 'Velocities' in line:
                e = i - 1
            elif 'Bonds' in line:
                s1 = i + 2
            elif 'Angles' in line:
                e1 = i - 1

        # Store the atom number and atom type.
        natoms = len(df[s:e])
        atomdata = {}
        for atom in df[s:e]:
            tmp = atom.split()
            atomdata[tmp[0]] = [int(tmp[2]),[]]

        # Also store the bond topology.  Note that only only the "central"
        # atom in the molecule has bond information.
        for atom in df[s1:e1]:
            tmp = atom.split()
            atomdata[tmp[2]][1].append(int(tmp[3]))        

        # Clear out the tuple (free memory)
        df = ()

        # What we have now
        # 1. lammps_type is a dictionary containing the correspondence between 
        #    RAPTOR atom types and LAMMPS atom ID numbers
        # 2. reaction_mol_num defines the correspondence between reaction names,
        #    molecule names, and molecule numbers ("kernels").
        # 3. atomdata defines the atom number, its LAMMPS atom ID, and which atoms
        #    it is bonded to (if it is a central atom in the molecule).
        # 
        # This information can be used to print the bond topology file.  The columns
        # in the topology file are:
        # atom number | molecule number | RAPTOR atom type

        # Molecule to atom name correspondence
        mol2atom = { 'OW' : 'H2O',
                     'HW' : 'H2O',
                     'OH' : 'H3O', 
                     'HO' : 'H3O' }

        # Based on the atomtype of the central atom, it is now possible to determine
        # the types of the terminal atoms.
        ltop = [[0,0,0] for row in range(natoms)]
        for key in atomdata.keys():
            ltop[int(key)-1][0] = int(key)

        #fi = open('top.tmp', 'w')

        # Based on the bonding topology, we can determine the molecule number ("kernel").
        for key in atomdata.keys():
            # Identify central atoms.  These contain bonding information.
            if atomdata[key][1] != []:
                #print(key, atomdata[key], file=fi)
                # List of bonded partners
                tmp = atomdata[key][1]
                # Type of the central atom
                ctype = atomdata[key][0]
                # Central atom
                # Here, we match ctype with what is contained in lammps_type.
                for key2 in lammps_type.keys():
                    if lammps_type[key2] == ctype:
                        mol = mol2atom[key2]
                        # Subtraction is accounting for Python numbering.
                        ltop[int(key)-1][1] = mol
                        # Central atom is assumed to always have index 1 for convenience.
                        # If this wasn't true, this part would be a lot more complicated.
                        ltop[int(key)-1][2] = 1
                # Terminal atoms.  Start the counter at 2 since the central atom is 1.
                count = 2
                for item in tmp:
                    # Account for Python numbering
                    i = item - 1
                    # Terminal atom type to be compared with lammps_type
                    ttype = atomdata[str(item)][0]
                    for key3 in lammps_type.keys():
                        if lammps_type[key3] == ttype:
                            mol = mol2atom[key3]
                            # Temporarily place the molecule name in the list.
                            ltop[i][1] = mol
                            # Place the counter in the list.  This is trivial as
                            # long as the terminal atoms are identical.  We must
                            # be careful in situations where this isn't true.
                            ltop[i][2] = count
                            count += 1

        # Finally we convert molecule strings to molecule numbers
        for i in range(len(ltop)):
            for reaction in reaction_mol_num.keys():
                for key in reaction_mol_num[reaction].keys():
                    ltop[i][1] = reaction_mol_num[reaction][key][ltop[i][1]] 

        # Output data to a file
        stop = 'lammps_msevb.top'
        top = open(stop, 'w')
        
        fmt = '{0:>4} {1:>2} {2:>2}'
        for i in range(len(ltop)):
            print(fmt.format(ltop[i][0], ltop[i][1], ltop[i][2]), file=top)
    elif cp2k:
        # Trajectory data
        traj_data = coord_and_force[snapshot]

        # Print out box bounds
        bb = box_bounds[snapshot]
        x = bb[0][1] - bb[0][0]
        y = bb[1][1] - bb[1][0]
        z = bb[2][1] - bb[2][0]
        fmt = ' ABC {0:7.5f} {1:7.5f} {2:7.5f}'
        print('&CELL')
        print(fmt.format(x,y,z))
        print('&END CELL')
        print('#')
        print('&COORD')
        lnstyle = ' {0:2} {1:10.5f} {2:10.5f} {3:10.5f}'
        for ln in traj_data:
            print(lnstyle.format(ln[1],ln[2],ln[3],ln[4]))
        print('&END COORD')

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
