from __future__ import division, print_function
import sys,os
import units_and_strings as uas
from numpy import array, append, average, std

def collect(fh):
    '''Collect data based on the type of file present.'''

    ###########################
    # Identify the type of file
    ###########################

    # Initialize logical variables
    xyz_file = False
    g09_input = False
    lammps_output = False
    lammps_data = False
    cp2k_output = False
    evb_output = False
    g09_output = False

    # Collect all data into memory for data retention
    with open(fh) as fl:
        f = tuple([line.rstrip() for line in fl])

    # Assuming we continually update RAPTOR for LAMMPS, we'll use newer versions
    # of the latter.  I store all versions I've used here for backwards 
    # compatability.
    lammps_version = ('LAMMPS (24 Apr 2013)',
                      'LAMMPS (5 Sep 2013)')

    # Determine if the file is an XYZ file
    if fh.split('.')[1] == 'xyz':
        xyz_file = True
    # Determine if the file is a Gaussian input
    elif fh.split('.')[1] == 'g09':
        g09_input = True
    # Check for Gaussian output
    elif 'Entering Gaussian System' in f[0]:
        g09_output = True
    # Check for LAMMPS data file
    elif 'LAMMPS data file' in f[0]:
        lammps_data = True
    # Check for EVB output
    elif 'REACTION_CENTER_LOCATION [ center_id | molecule_id | location_rank ]' in f:
        evb_output = True
    # Check for CP2K output
    elif 'CP2K' in f[18]:
        cp2k_output = True
    elif 'CP2K' in f[19]:
        cp2k_output = True
    elif 'CP2K' in f[20]:
        cp2k_output = True
    # Check for LAMMPS output
    else:
        for item in lammps_version:
            if item in f:
                lammps_output = True
                break

    ##################
    # Collect the data
    ##################

    # Initialize the data dictionary
    data = {}

    # Initialize the indices dictionary
    indices = {}

    # LAMMPS Output
    if lammps_output:
        # Strings to search for
        table = uas.lammps_output_strings() 
        # Determine which units were used.  This is important for correct output
        units_table = uas.lammps_units()
        
        # Read through the file
        for i,line in enumerate(f):
            if '-----' in line:
                tmp = line.split('-----')
                if tmp[1] in table.keys():
                    indices[table[tmp[1]]] = i
            # Thermodynamics data table
            if 'Step' in line:
                indices['THERMOTABLE'] = i

        s = indices['UNITS'] + 1
        for i,line in enumerate(f[s:]):
            if 'units' in line:
                tmp = line.split()
                data['UNITS'] = units_table[tmp[1].upper()]
                break

        # Determine the type of data stored in the thermodynamics table
        thermo = f[indices['THERMOTABLE']].split()
        # Here, data stores the column number, the data from the table,
        # the average, and standard deviation.
        for i,item in enumerate(thermo):
            data[item] = [i,array([],dtype=float),0.0,0.0]

        # Table of quantities for performing averages
        comp_table = ( 'Temp', 'TotEng', 'PotEng', )

        # Collect the data from the thermodynamics table
        s = indices['THERMOTABLE'] + 1
        e = [i+s for i,line in enumerate(f[s:]) if 'Loop time' in line][0]
        for i in range(s,e):
            tmp = f[i].split()
            for key in comp_table:
                data[key][1] = append(data[key][1], float(tmp[data[key][0]])) 

        # Determine averages and standard deviations of quantities
        for key in data.keys():
            if key in comp_table:
                data[key][2] = average(data[key][1])
                data[key][3] = std(data[key][1])

        # Add that this is a LAMMPS logfile
        data['FILETYPE'] = 'LAMMPS Logfile'
    elif lammps_data:
        props = ( 'atoms', 'bonds', 'angles', 'dihedrals', 'impropers', 
                  'atom types', 'bond types', 'angle types' 'dihedral types',
                  'improper types', )

        # Number of atoms, bonds, angles
        # Number of atomtypes, bond types, angle types
        for i in range(2,13):
            for item in props:
                if item in f[i]:
                    if 'types' in item:
                        tmp = item.split()
                        k = 'Num_' + tmp[0] + tmp[1]
                        data[k] = int(f[i].split()[0])
                    else:        
                        data['Num_' + item] = int(f[i].split()[0])

        # Locate atoms, bonds, and angles
        natoms = data['Num_atoms']
        for i,ln in enumerate(f):
            if 'Atoms' in ln:
                indices['ATOMS'] = [i+2,i+2+natoms]
            elif 'Bonds' in ln:
                indices['BONDS'] = [i+2,i+2+natoms]
            elif 'Angles' in ln:
                indices['ANGLES'] = [i+2,i+2+natoms]
            elif 'Velocities' in ln:
                indices['VELOCITIES'] = [i+2,i+2+natoms]

        # Collect atoms and charges
        data['Coords'] = {}
        data['Charges'] = {}
        s = indices['ATOMS'][0]
        e = indices['ATOMS'][1]
        for i in range(s,e):
            tmp = f[i].split()
            data['Coords'][tmp[0]] = [int(tmp[2]), float(tmp[4]), float(tmp[5]), float(tmp[6])]
            data['Charges'][tmp[2]] = float(tmp[3])

        # Collect velocities, if they are present
        if 'VELOCITIES' in indices.keys():
            data['Velocities'] = {}
            s = indices['VELOCITIES'][0]
            e = indices['VELOCITIES'][1]
            for i in range(s,e):
                tmp = f[i].split()
                data['Velocities'][tmp[0]] = [float(tmp[1]), float(tmp[2]), float(tmp[3])] 

        # Add that this is a LAMMPS datafile
        data['FILETYPE'] = 'LAMMPS Datafile'
    # RAPTOR Output
    elif evb_output:
        # Add that this is a EVB output
        data['FILETYPE'] = 'RAPTOR Output'

        # Strings to search for.  Store as line type keys and block type keys.
        # This will need to be reworked for the situation where multiple complexes
        # exist.
        table = uas.raptor_output_strings() 

        # Read through the file
        # Locate lines storing timesteps.  These are stored for the block type keys
        # in indices.  The dictionary data stores information from the RAPTOR
        # simulation.

        # Initialize lists for data storage
        for key in table['LINE'].keys():
            data[table['LINE'][key]] = []
        for key in table['BLOCK'].keys():
            data[table['BLOCK'][key]] = []
            # We need the output file line indices for collecting block type keys
            indices[table['BLOCK'][key]] = []

        # By reading through the entire file 1 time here, we don't spend
        # much time when collecting blocks later.
        for i,line in enumerate(f):
            tmp = line.split()
            if tmp == []:
                continue
            else:
                # Collect data for the line type keys.  These can be easily
                # collected while reading the entire file.
                if 'COMPLEX' in tmp and 'state(s)' in tmp:
                    data[table['LINE'][tmp[3]]].append(float(tmp[2]))
                elif tmp[0] in table['LINE'].keys():
                    data[table['LINE'][tmp[0]]].append(float(tmp[1]))
                # Get the line numbers for the start of a block.
                if tmp[0] in table['BLOCK'].keys():
                    indices[table['BLOCK'][tmp[0]]].append([i+1])
                # Get the line numbers for the end of a block.  This is added
                # to the list generated for the start of a block.
                if tmp[0] in table['B_END'].keys():
                    indices[table['B_END'][tmp[0]]][-1].append(i)

        # Collect block information piece by piece.
        # Reaction center location
        if 'RXNCENTER' in indices.keys():
            for block in indices['RXNCENTER']:
                s = block[0]
                e = block[1]
                for line in f[s:e]:
                    tmp = line.split() 
                    ltmp = [int(tmp[0]), int(tmp[1]), int(tmp[2])]
                    data['RXNCENTER'].append(ltmp)

        # Environment energy decomposition
        if 'ENV_ENERGY_DECOMP' in indices.keys():
            for block in indices['ENV_ENERGY_DECOMP']:
                s = block[0]
                e = block[1]
                for line in f[s:e]:
                    tmp = line.split() 
                    ltmp = [float(tmp[0]), float(tmp[1]), float(tmp[2]),
                            float(tmp[3]), float(tmp[4]), float(tmp[5]), float(tmp[6])]
                    data['ENV_ENERGY_DECOMP'].append(ltmp)

        # EVB states
        if 'STATES' in indices.keys():
            for block in indices['STATES']:
                s = block[0]
                e = block[1]
                for line in f[s:e]:
                    tmp = line.split() 
                    ltmp = [int(tmp[0]), int(tmp[1]), int(tmp[2]),
                            int(tmp[3]), int(tmp[4]), int(tmp[5]),
                            int(tmp[6]), int(tmp[7])]
                    data['STATES'].append(ltmp)

        # System energy (diagonal portion)
        if 'ENERGY_DIAGONAL' in indices.keys():
            for block in indices['ENERGY_DIAGONAL']:
                s = block[0]
                e = block[1]
                data['ENERGY_DIAGONAL'].append([])
                for line in f[s:e]:
                    tmp = line.split() 
                    ltmp = [int(tmp[0]), float(tmp[1]), float(tmp[2]),
                            float(tmp[3]), float(tmp[4]), float(tmp[5]),
                            float(tmp[6]), float(tmp[7]), float(tmp[8]),
                            float(tmp[9])]
                    data['ENERGY_DIAGONAL'][-1].append(ltmp)

        # System energy (off-diagonal portion)
        if 'ENERGY_OFF_DIAGONAL' in indices.keys():
            for block in indices['ENERGY_OFF_DIAGONAL']:
                s = block[0]
                e = block[1]
                data['ENERGY_OFF_DIAGONAL'].append([])
                for line in f[s:e]:
                    tmp = line.split() 
                    ltmp = [int(tmp[0]), float(tmp[1]), float(tmp[2]),
                            float(tmp[3]), float(tmp[4])]
                    data['ENERGY_OFF_DIAGONAL'][-1].append(ltmp)

        # Ignore extra coupling for the moment.

        # CI vector
        if 'CI_VECTOR' in indices.keys():
            for block in indices['CI_VECTOR']:
                s = block[0]
                e = block[1]
                for line in f[s:e]:
                    tmp = line.split() 
                    # We use a list comprehension because the number of CI coefficients
                    # depends on the number of MS-EVB states, and therefore has variable
                    # length.  We store the CI coefficients and number of MS-EVB states
                    # with coefficients greater than 0.001.
                    ltmp = [float(tmp[i]) for i in range(len(tmp)) if float(tmp[i]) >= 0.001]
                    data['CI_VECTOR'].append([len(ltmp),ltmp])
 
        # Center of excess charge (CEC)
        if 'CEC' in indices.keys:
            for block in indices['CEC']:
                s = block[0]
                e = block[1]
                for line in f[s:e]:
                    tmp = line.split() 
                    ltmp = [float(tmp[0]), float(tmp[1]), float(tmp[2])]
                    data['CEC'].append(ltmp)

        # CEC V2 coordinate
        if 'CEC_V2' in indices.keys():
            for block in indices['CEC_V2']:
                s = block[0]
                e = block[1]
                for line in f[s:e]:
                    tmp = line.split() 
                    ltmp = [float(tmp[0]), float(tmp[1]), float(tmp[2])]
                    data['CEC_V2'].append(ltmp)
    # XYZ file
    elif xyz_file:
        # Add that this is a XYZ file
        data['FILETYPE'] = 'XYZ file'

        # Number of atoms
        data['NATOMS'] = int(f[0])

        # Coordinates
        data['ATOMS'] = []
        data['COORDS'] = []
        for i in range(2,len(f)):
            ln = f[i].split()
            data['ATOMS'].append(ln[0])
            data['COORDS'].append([ float(ln[1]), float(ln[2]), float(ln[3]) ])
    # CP2K output
    elif cp2k_output:
        # Add that this is a CP2K output
        data['FILETYPE'] = 'CP2K Output'

        # Locate the data.
        for i,ln in enumerate(f):
            if 'CELL| Volume [angstrom^3]:' in ln:
                indices['CELL'] = i
            elif 'DFT| Spin restricted Kohn-Sham (RKS) calculation' in ln:
                indices['DFT_INFO'] = i
            elif 'FUNCTIONAL|' in ln:
                if 'DFT_FXN' not in indices.keys():
                    indices['DFT_FXN'] =[]
                    indices['DFT_FXN'].append(i)
                else:
                    indices['DFT_FXN'].append(i)
            elif 'MODULE QUICKSTEP:  ATOMIC COORDINATES IN angstrom' in ln:
                indices['COORD_INITIAL_START'] = i + 4
            elif 'SCF PARAMETERS' in ln:
                indices['COORD_INITIAL_END'] = i - 4
            elif 'MD_ENERGIES| Initialization proceeding' in ln:
                indices['MD_1STSTEP'] = i+4
            elif 'ENSEMBLE TYPE' in ln:
                if 'TIMESTEP' not in indices.keys():
                    indices['TIMESTEP'] = []
                    indices['TIMESTEP'].append(i+2)
                else:
                    indices['TIMESTEP'].append(i+2)
            #elif 'Total energy:' in ln:
            #    # Don't get the first time, which corresponds to the
            #    # atomic density guess.
            #    if 'ENERGIES' not in indices.keys():
            #        indices['ENERGIES'] = []
            #    else:
            #        indices['ENERGIES'].append(i-6)
            elif 'SCF run converged in' in ln:
                # Get energies from FORCE_EVAL
                indices['ENERGIES'] = i+8
                data['ENERGIES'] = {}
            elif 'Total FORCE_EVAL ( QS )' in ln:
                # Total QM energy in FORCE_EVAL calculation
                indices['TOTAL_QM_ENERGY'] = i
            elif 'Total FORCE_EVAL ( QMMM )' in ln:
                # Total QMMM energy in FORCE_EVAL calculation
                indices['TOTAL_QMMM_ENERGY'] = i
            elif 'MULLIKEN POPULATION ANALYSIS' in ln:
                # Don't get the first time, which corresponds to the
                # atomic density guess.
                if 'MULLIKEN' not in indices.keys():
                    indices['MULLIKEN'] = []
                else:
                    indices['MULLIKEN'].append(i+3)
            elif 'ATOMIC FORCES in [a.u.]' in ln:
                # Don't get the first time, which corresponds to the
                # atomic density guess.
                if 'FORCES' not in indices.keys():
                    indices['FORCES'] = []
                else:
                    indices['FORCES'].append(i+3)

        # Collect the cell parameters
        if 'CELL' in indices.keys():
            data['CELL'] = {}
            ind = indices['CELL']
            data['CELL']['VOLUME'] =  float(f[ind].split()[3])
            data['CELL']['A_VEC'] =  [ float(f[ind+1].split()[4]),
                                       float(f[ind+1].split()[5]),
                                       float(f[ind+1].split()[6]),
                                       float(f[ind+1].split()[9]) ]
            data['CELL']['B_VEC'] =  [ float(f[ind+2].split()[4]),
                                       float(f[ind+2].split()[5]),
                                       float(f[ind+2].split()[6]),
                                       float(f[ind+2].split()[9]) ]
            data['CELL']['C_VEC'] =  [ float(f[ind+3].split()[4]),
                                       float(f[ind+3].split()[5]),
                                       float(f[ind+3].split()[6]),
                                       float(f[ind+3].split()[9]) ]
            data['CELL']['BC_ANGLE'] = float(f[ind+4].split()[5])
            data['CELL']['AC_ANGLE'] = float(f[ind+5].split()[5])
            data['CELL']['AB_ANGLE'] = float(f[ind+6].split()[5])

        # Collect the DFT information
        if 'DFT' in indices.keys():
            data['DFT'] = {}
            ind = indices['DFT_INFO']
            data['DFT']['TYPE'] = f[ind].split()[6]
            data['DFT']['MULT'] = int(f[ind+1].split()[2])
            data['DFT']['CHARGE'] = int(f[ind+3].split()[2])
            data['DFT']['SIC'] = f[ind+4].split()[4]
            data['DFT']['DENSITY_SMOOTHING'] = f[ind+9].split()[4]
            data['DFT']['XC_DERIV'] = f[ind+10].split()[3]

            # Collect the XC-functional
            fxn_list = { 'BECKE88' : 'B', 
                         'LYP'     : 'LYP', 
                         'P86C'    : 'P', 
                         'PADE'    : 'Pade LDA', 
                         'PBE'     : 'PBE', }
            fnl = ''
            for i in indices['DFT_FXN']:
                tmp = f[i].split()[1][:-1]
                if tmp in fxn_list.keys():
                    fnl += fxn_list[tmp]
    
            data['DFT']['FUNCTIONAL'] = fnl

        # Collect the initial coordinates
        if 'COORD_INITIAL_START' in indices.keys():
            s = indices['COORD_INITIAL_START']
            e = indices['COORD_INITIAL_END']
            natoms = e - s
            data['NATOMS'] = natoms

            data['INITIAL_COORDS'] = []
            for i in range(s,e):
                ln = f[i].split()
                data['INITIAL_COORDS'].append([ ln[2], float(ln[4]), float(ln[5]), float(ln[6]) ])

        # Initialize the timestep information
        if 'TIMESTEP' in indices.keys():
            data['TS'] = {}
            # Initial step is unique
            data['TS']['0'] = { 'TIME' : 0.000, }
            for i in indices['TIMESTEP']:
                step = f[i-1].split()[3]
                time = float(f[i].split()[3])
                data['TS'][step] = { 'TIME' : time, }

        # Collect energy information
        if 'TIMESTEP' in indices.keys() and 'ENERGIES' in indices.keys():
            ix = 0
            for i in indices['ENERGIES']:
                data['TS'][str(ix)]['ENERGY'] = []
                ovlp = float(f[i].split()[7])
                selfe = float(f[i+1].split()[7])
                core = float(f[i+2].split()[3])
                hartree = float(f[i+3].split()[2])
                xc = float(f[i+4].split()[2])
                total = float(f[i+6].split()[2])
            
                data['TS'][str(ix)]['ENERGY'].append([ ovlp, selfe, core, hartree, xc, total ])

                ix += 1
        else:
            s = indices['ENERGIES']
            e = [i+s+1 for i,ln in enumerate(f[s:]) if 'Total energy:' in ln][0]
            data['ENERGIES']['OVERLAP_ENERGY'] = float(f[s].split()[7])
            data['ENERGIES']['CORE_SELF_ENERGY'] = float(f[s+1].split()[7])
            data['ENERGIES']['CORE_HAMILTONIAN'] = float(f[s+2].split()[3])
            data['ENERGIES']['COULOMB_ENERGY'] = float(f[s+3].split()[2])
            data['ENERGIES']['DFT_XC_ENERGY'] = float(f[s+4].split()[2])
            # For functionals with HF exchange
            if 'Hartree-Fock Exchange' in f[s+5]:
                data['ENERGIES']['HFX_ENERGY'] = float(f[s+5].split()[3])
                # If dispersion energy is included
                if 'Dispersion' in f[s+6]:
                    data['ENERGIES']['DISPERSION_ENERGY'] = float(f[s+6].split()[2])
                    if 'QM/MM Electrostatic' in f[s+7]:
                        data['ENERGIES']['QMMM_ELECTROSTATIC_ENERGY'] = float(f[s+7].split()[3])
                elif 'QM/MM Electrostatic' in f[s+6]:
                    data['ENERGIES']['QMMM_ELECTROSTATIC_ENERGY'] = float(f[s+6].split()[3])
            # For functionals with empirical dispersion
            elif 'Dispersion' in f[s+5]: 
                data['ENERGIES']['DISPERSION_ENERGY'] = float(f[s+5].split()[2])
                if 'QM/MM Electrostatic' in f[s+6]:
                    data['ENERGIES']['QMMM_ELECTROSTATIC_ENERGY'] = float(f[s+6].split()[3])
            # For QMMM with no HF exchange or empirical dispersion
            elif 'QM/MM Electrostatic' in f[s+5]:
                data['ENERGIES']['QMMM_ELECTROSTATIC_ENERGY'] = float(f[s+5].split()[3])
            # QM system total energy
            i = indices['TOTAL_QM_ENERGY']
            data['ENERGIES']['QM_TOTAL_ENERGY'] = float(f[i].split()[8])
            # QMMM total energy
            if 'TOTAL_QMMM_ENERGY' in indices.keys():
                i = indices['TOTAL_QMMM_ENERGY']
                data['ENERGIES']['QMMM_TOTAL_ENERGY'] = float(f[i].split()[8])

        # Collect Mulliken charges
        if 'MULLIKEN' in indices.keys():
            ix = 0
            for i in indices['MULLIKEN']:
                data['TS'][str(ix)]['MULLIKEN'] = []

                s = i
                e = i + natoms
                for j in range(s,e):
                    ln = f[j].split()
                    data['TS'][str(ix)]['MULLIKEN'].append([ ln[1], float(ln[3]), float(ln[4]) ])

                ix += 1

        # Collect Forces
        if 'FORCES' in indices.keys():
            ix = 0
            for i in indices['FORCES']:
                data['TS'][str(ix)]['FORCE'] = []

                # Get the forces on all atoms
                s = i
                e = i + natoms
                for j in range(s,e):
                    ln = f[j].split()
                    data['TS'][str(ix)]['FORCE'].append([ ln[2], float(ln[3]), float(ln[4]), float(ln[5]) ])
                
                # Get the sum of forces on all atoms and its magnitude
                j = e 
                ln = f[j].split()
                data['TS'][str(ix)]['FORCE_SUM'] = [ float(ln[4]), float(ln[5]), float(ln[6]), float(ln[7]) ]

                ix += 1

        # Store the filename
        data['FILENAME'] = fh
    # Gaussian input
    elif g09_input:
        # Note the filetype
        data['FILETYPE'] = 'G09_INPUT'
        # Determine the locations data in the file
        # Grab the type of job
        data['JOBTYPE'] = [ln for ln in f if '#' in ln]
        # Grab the title
        if f[2] != '':
            data['TITLE'] = f[2]
        else:
            data['TITLE'] = f[3]
        for i,ln in enumerate(f):
            tmp = ln.split()
            try:
                if ln != '':
                    one = int(tmp[0])
                    two = int(tmp[1])
                    s = i + 1
                    e = [s+j for j,l in enumerate(f[s:]) if l == ''][0]
                    indices['COORD_START'] = s
                    indices['COORD_END'] = e
            except ValueError:
                continue

        # Collect coordinates
        s = indices['COORD_START']
        e = indices['COORD_END']

        # Collect charge and multiplicity
        data['CHARGE_MULT'] = f[s-1]

        # Define some relevant parameters
        data['NATOMS'] = e - s 
        data['ATOMS'] = []
        data['COORDS'] = []
        # Decide if the file is written in Cartesian format or Z-matrix format.
        tmp = f[indices['COORD_START']].split()
        data['ZMAT'] = False
        if len(tmp) == 1:
            data['ZMAT'] = True
            data['BONDS'] = []
            data['ANGLES'] = []
            data['DIHEDRALS'] = []

        for ln in f[s:e]:
            tmp = ln.split()
            # Z-matrix
            if data['ZMAT']:
                data['ATOMS'].append(tmp[0])
                if len(tmp) == 1:
                    data['COORDS'].append([])
                else:
                    data['BONDS'].append(float(tmp[2]))
                    # Line contains only bonding
                    if len(tmp) == 3:
                        data['COORDS'].append([ tmp[1], float(tmp[2]) ])  
                    # Line contains bond and angle          
                    elif len(tmp) == 5:
                        data['COORDS'].append([ tmp[1], float(tmp[2]), tmp[3], float(tmp[4]) ])            
                        data['ANGLES'].append(float(tmp[4]))
                    # Line contains bond, angle, and dihedral
                    elif len(tmp) == 7:
                        data['COORDS'].append([ tmp[1], float(tmp[2]), 
                                                tmp[3], float(tmp[4]),
                                                tmp[5], float(tmp[6]) ])            
                        data['ANGLES'].append(float(tmp[4]))
                        data['DIHEDRALS'].append(float(tmp[6]))
            # Cartesian
            else:
                data['ATOMS'].append(tmp[0])
                data['COORDS'].append([ float(tmp[1]), float(tmp[2]), float(tmp[3]) ])
    # Gaussian output
    elif g09_output:
        # Note the filetype
        data['FILETYPE'] = 'G09_OUTPUT'

        # Store the charges of the system
        data['CHARGE_MODEL'] = {}

        # Store atomic dipoles
        data['ATOMIC_DIPOLE'] = {}

        # Store the filename
        data['FILENAME'] = fh

        # Determine locations of data
        for i,ln in enumerate(f):
            # Total energy from Hartree-Fock or DFT.
            if 'SCF Done:' in ln:
                if 'SCF_ENERGY' not in indices.keys():
                    indices['SCF_ENERGY'] = []
                    indices['SCF_ENERGY'].append(i)
                else:
                    indices['SCF_ENERGY'].append(i)
            # Bond lengths and angles in any calculation
            elif 'Symbolic Z-matrix:' in ln:
                indices['STRUCT_PARAMS'] = []
                indices['STRUCT_PARAMS'].append(i+2)
            # Bond lengths and angles in geometry optimization
            elif 'Optimized Parameters' in ln:
                indices['STRUCT_PARAMS_OPT'] = True 
            # MP2 energy
            elif 'EUMP2' in ln:
                if 'MP2_ENERGY' not in indices.keys():
                    indices['MP2_ENERGY'] = []
                    indices['MP2_ENERGY'].append(i)
                else:
                    indices['MP2_ENERGY'].append(i)
            # Mulliken charges
            elif 'Mulliken atomic charges:' in ln:
                if 'MULLIKEN_CHARGE' not in indices.keys():
                    indices['MULLIKEN_CHARGE'] = []
                    indices['MULLIKEN_CHARGE'].append(i+2)
                    data['CHARGE_MODEL']['Mulliken'] = []
                else:
                    indices['MULLIKEN_CHARGE'].append(i+2)
            # NPA charges
            elif 'Summary of Natural Population Analysis:' in ln:
                if 'NPA_CHARGE' not in indices.keys():
                    indices['NPA_CHARGE'] = []
                    indices['NPA_CHARGE'].append(i+6)
                    data['CHARGE_MODEL']['NPA'] = []
                else:
                    indices['NPA_CHARGE'].append(i+6)
            # Hirshfeld charges
            elif 'Hirshfeld spin densities, charges and dipoles' in ln:
                if 'HIRSHFELD_CHARGE' not in indices.keys():
                    indices['HIRSHFELD_CHARGE'] = []
                    indices['HIRSHFELD_CHARGE'].append(i+2)
                    data['CHARGE_MODEL']['HIRSHFELD'] = []
                    data['ATOMIC_DIPOLE']['HIRSHFELD'] = []
                else:
                    indices['HIRSHFELD_CHARGE'].append(i+2)
            # ESP charges
            elif 'Fitting point charges to electrostatic potential' in ln:
                if 'ESP_CHARGE' not in indices.keys():
                    indices['ESP_CHARGE'] = []
                    indices['ESP_CHARGE'].append(i+4)
                else:
                    indices['ESP_CHARGE'].append(i+4)
            # Hu-Lu-Yang ESP model
            elif 'Generate Potential Derived Charges using the Hu-Lu-Yang model:' in ln:
                data['CHARGE_MODEL']['Hu-Lu-Yang'] = []
            # Merz-Singh-Kollman ESP model
            elif 'Merz-Kollman atomic radii used.' in ln:
                data['CHARGE_MODEL']['Merz-Singh-Kollman'] = []
            # Reoriented coordinates
            #elif 'Standard orientation:' in ln:
            elif 'Input orientation:' in ln:
                if 'CART_COORDS' not in indices.keys():
                    indices['CART_COORDS'] = []
                    indices['CART_COORDS'].append(i+5)
                else:
                    indices['CART_COORDS'].append(i+5)
            # Finite field
            elif 'The following finite field(s) will be applied:' in ln:
                indices['FIELD'] = i + 1
            # Polarizable Continuum Model
            elif 'After PCM corrections,' in ln:
                indices['PCM'] = i
            # Atomic forces
            elif 'Forces (Hartrees/Bohr)' in ln:
                if 'FORCES' not in indices.keys():
                    indices['FORCES'] = []
                    indices['FORCES'].append(i + 3)
                else:
                    indices['FORCES'].append(i + 3)

        # Get the structural parameters (bonds, angles, dihedrals)
        if 'STRUCT_PARAMS_OPT' in indices.keys():
            # Gaussian prints by atomic numbers, rather than element names.
            atomconv = { '1'  : 'H',
                         '2'  : 'He',
                         '3'  : 'Li',
                         '4'  : 'Be',
                         '5'  : 'B',
                         '6'  : 'C',
                         '7'  : 'N',
                         '8'  : 'O',
                         '9'  : 'F',
                         '10' : 'Ne',
                       }

            s = indices['CART_COORDS'][-1]
            strsearch = '---------------------------------------------------------------------'
            e = [ s+i for i,ln in enumerate(f[s:]) if strsearch in ln ][0]
            data['ATOMS'] = []
            data['COORDS'] = []
            for i in range(s,e):
                tmp = f[i].split()
                data['ATOMS'].append(atomconv[tmp[1]])
                data['COORDS'].append([ float(tmp[3]), float(tmp[4]), float(tmp[5]) ])
            data['NATOMS'] = len(data['ATOMS'])
        else:
            s = indices['STRUCT_PARAMS'][-1]
            # Determine where the block ends
            for i,ln in enumerate(f[s:]):
                e = s + i
                if ln == '':
                    break
            # Collect atoms and structural parameters
            data['ATOMS'] = []
            data['COORDS'] = []
            data['BONDS'] = []
            data['ANGLES'] = []
            data['DIHEDRALS'] = []
            tmp = f[s].split()
            if len(tmp) == 1:
                data['ZMATRIX'] = True
            else:
                data['ZMATRIX'] = False

            for i in range(s,e):
                tmp = f[i].split()
                if data['ZMATRIX'] == True:
                    data['ATOMS'].append(tmp[0])
                    if len(tmp) == 3:
                        data['COORDS'].append([ tmp[1], float(tmp[2]) ])
                        data['BONDS'].append(float(tmp[2]))
                    elif len(tmp) == 5:
                        data['COORDS'].append([ tmp[1], float(tmp[2]),
                                                tmp[3], float(tmp[4]) ])
                        data['BONDS'].append(float(tmp[2]))
                        data['ANGLES'].append(float(tmp[4]))
                    elif len(tmp) == 8:
                        data['COORDS'].append([ tmp[1], float(tmp[2]),
                                                tmp[3], float(tmp[4]),
                                                tmp[5], float(tmp[6]) ])
                        data['BONDS'].append(float(tmp[2]))
                        data['ANGLES'].append(float(tmp[4]))
                        data['DIHEDRALS'].append(float(tmp[6]))
                else:
                    data['ATOMS'].append(tmp[0])
                    data['COORDS'].append([ float(tmp[1]), float(tmp[2]), float(tmp[3]) ])
            data['NATOMS'] = len(data['ATOMS'])

        # Get the Cartesian coordinates
        data['CART_COORDS'] = []
        s = indices['CART_COORDS'][-1]
        # Determine the end point of the search
        check = '-'*69
        for i,ln in enumerate(f[s:]):
            if check in ln: break
            tmp = ln.split()
            data['CART_COORDS'].append([float(tmp[3]), float(tmp[4]), float(tmp[5])])

        # Get the total SCF energy of the system.
        data['ENERGY'] = {}
        ln = f[indices['SCF_ENERGY'][-1]].split()
        data['ENERGY']['SCF'] = float(ln[4])

        # Get the energy from the post-HF calculation.
        if 'MP2_ENERGY' in indices.keys():
            ln = f[indices['MP2_ENERGY'][-1]].split()
            ln[2] = ln[2].replace('D','E')
            data['ENERGY']['CORR'] = float(ln[2])
            ln[5] = ln[5].replace('D','E')
            data['ENERGY']['MP2'] = float(ln[5])

        # Get the Mulliken charges of the system
        if 'MULLIKEN_CHARGE' in indices.keys():
            s = indices['MULLIKEN_CHARGE'][-1]
            e = s + data['NATOMS']
            for ln in f[s:e]:
                tmp = ln.split()
                data['CHARGE_MODEL']['Mulliken'].append(float(tmp[2]))

        # Get the NPA charges of the system
        if 'NPA_CHARGE' in indices.keys():
            s = indices['NPA_CHARGE'][-1]
            e = s + data['NATOMS']
            for ln in f[s:e]:
                tmp = ln.split()
                data['CHARGE_MODEL']['NPA'].append(float(tmp[2]))

        # Get the Hirshfeld charges of the system
        if 'HIRSHFELD_CHARGE' in indices.keys():
            s = indices['HIRSHFELD_CHARGE'][-1]
            e = s + data['NATOMS']
            for ln in f[s:e]:
                tmp = ln.split()
                data['CHARGE_MODEL']['HIRSHFELD'].append(float(tmp[3]))
                data['ATOMIC_DIPOLE']['HIRSHFELD'].append([ float(tmp[4]), float(tmp[5]), float(tmp[6]) ]) 

        # Get the ESP charges of the system
        if 'ESP_CHARGE' in indices.keys():
            s = indices['ESP_CHARGE'][-1]
            e = s + data['NATOMS']
            for ln in f[s:e]:
                tmp = ln.split()
                if 'Hu-Lu-Yang' in data['CHARGE_MODEL'].keys():
                    data['CHARGE_MODEL']['Hu-Lu-Yang'].append(float(tmp[2]))
                elif 'Merz-Singh-Kollman' in data['CHARGE_MODEL'].keys():
                    data['CHARGE_MODEL']['Merz-Singh-Kollman'].append(float(tmp[2]))

        # Collect the external electric field
        if 'FIELD' in indices.keys():
            i = indices['FIELD']
            tmp = f[i].split()
            data['FIELD'] = [ float(tmp[4]), float(tmp[5]), float(tmp[6]) ] 
        else:
            data['FIELD'] = [ 0.000, 0.000, 0.000 ] 

        # Collect the PCM energy
        if 'PCM' in indices.keys():
            i = indices['PCM']
            tmp = f[i].split()
            data['ENERGY']['PCM'] = float(tmp[6])

        # Collect the atomic forces (in Hartrees/Bohr)
        if 'FORCES' in indices.keys():
            # Get the first force since that is what we need for force matching
            s = indices['FORCES'][0]
            e = s + data['NATOMS'] 
            data['FORCES'] = [] 
            for ln in f[s:e]:
                tmp = ln.split()
                data['FORCES'].append( [ float(tmp[2]), 
                                         float(tmp[3]),
                                         float(tmp[4]) ] )

    return data
