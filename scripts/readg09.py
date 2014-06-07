#! /usr/bin/env python

from __future__ import division, print_function
import sys
import numpy as np
from collect import collect
from textwrap import dedent

def main():
    """\
    Script for reading Gaussian input and output files.
    """

    from argparse import ArgumentParser, RawDescriptionHelpFormatter
    from textwrap import dedent
    parser = ArgumentParser(description=dedent(main.__doc__),
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('files', help='The files that you want converted.',
                        nargs='+')
    #parser.add_argument('-cb', '--changebond', help='Change bonds in the molecule.')
    #parser.add_argument('-ca', '--changeangle', help='Change angles in the molecule.')
    #parser.add_argument('-cd', '--changeangle', help='Change angles in the molecule.')
    parser.add_argument('-ch', '--charges', help='Output charges.',
                        action='store_true', default=False)
    parser.add_argument('-e', '--energies', help='Output energies.',
                        action='store_true', default=False)
    parser.add_argument('-c', '--coords', help='Output coordinates.',
                        action='store_true', default=False)
    parser.add_argument('-f', '--forces', help='Output forces.',
                        action='store_true', default=False)
    parser.add_argument('-tsc', '--tscoords', help='Output coordinates from transition state search.',
                        action='store_true', default=False)
    parser.add_argument('-u', '--units', help='Units for outputting quantities.',
                        default='au')
    parser.add_argument('-p', '--plot', help='Use this to make a plot of data.',
                        action='store_true', default=False)
    parser.add_argument('-v', '--variable', help='Choose the variable for a plot (x-axis) or a dihedral scan.',
                        default=None)
    parser.add_argument('-ds', '--dihedralscan', help='Choose the min, max, and stride for a dihedral scan.',
                        nargs=3, required=False)
    parser.add_argument('-ad', '--atomicdipole', help='Output the atomic dipoles from Hirshfeld analysis.',
                        action='store_true', default=False)
    parser.add_argument('-ap', '--atomicpol', help='Output the atomic polarizabilities.',
                        action='store_true', default=False)
    parser.add_argument('-fe', '--fitevb', help='Output energy data for FitEVB',
                        action='store_true', default=False)
    parser.add_argument('-ba', '--bondedatoms', help='Read in bonded atoms',
                        nargs=2, required=False)
    args = parser.parse_args()

    # Store command line arguments as convenient variables
    fh = args.files
    charges = args.charges
    energies = args.energies
    coords = args.coords
    forces = args.forces
    tscoords = args.tscoords
    units = args.units
    plot = args.plot
    var = None
    atomdip = args.atomicdipole
    atompol = args.atomicpol
    if args.variable is not None: var = args.variable.upper()
    if args.dihedralscan is not None: 
        dscan = [ float(args.dihedralscan[0]), 
                  float(args.dihedralscan[1]), 
                  float(args.dihedralscan[2]) ]
    else:
        dscan = []
    fitevb = args.fitevb
    bondedatoms = args.bondedatoms

    if (plot == True and var == None) or (dscan != [] and var == None):
        sys.exit('Must use -p and -v flags or -ds and -v flags together!')

    # Options for plotting
    # BL = bond length
    # AN = angle
    # DH = dihedral
    plot_opts = ( 'BL1', 'BL2', 'BL3', 'BL4', 'BL5', 'BL6', 'BL7', 
                  'AN1', 'AN2', 'AN3', 'AN4',  
                  'DH1', 'DH2', 'DH3' )

    # Check if the requested plot option makes sense
    if plot:
        if var not in plot_opts:
            sys.exit('Cannot print ' + var + '!  Must choose a bond length, angle, or dihedral!') 
        # Determine the type of data desired
        if 'BL' in var:
            xdata = 'BONDS'
        elif 'AN' in var:
            xdata = 'ANGLES'
        elif 'DH' in var:
            xdata = 'DIHEDRALS'
        # Make lists for storing data
        x = []
        y = []

    # Variables for storing all quantities related to atomic polarizabilities
    if atompol:
        # Electric fields
        ff = []
        # Atomic dipoles 
        ad = []

    # Energy label conversion dictionary
    table = { 'SCF'  : 'SCF Energy',
              'CORR' : 'Correlation Energy',
              'MP2'  : 'MP2 Energy',
              'PCM'  : 'Total Energy + PCM',
            }

    # Unit conversion dictionary
    conv = { 'au'      : ['Hartree' ,      1.00000000],
             'ev'      : ['eV'      ,     27.21138505],
             'kcalmol' : ['kcal/mol',    627.50947428],
             'kjmol'   : ['kJ/mol'  ,   2625.49964038],
             'wvn'     : ['cm^{-1}' , 219474.63068   ],
           }

    # Conversion between Bohr and Angstroms
    au2angstrom = 0.52917720859

    # Dictionary for converting plot input to code input
    pdata = { 'BL1' : 0,
              'BL2' : 1,
              'BL3' : 2,
              'BL4' : 3,
              'BL5' : 4,
              'AN1' : 0,
              'AN2' : 1,
              'AN3' : 2,
              'AN4' : 3,
              'DH1' : 0,
              'DH2' : 1,
              'DH3' : 2,
            }

    # /\/\/\/\/\/\/\/\/
    # Collect the files 
    # /\/\/\/\/\/\/\/\/

    for f in fh:
        data = collect(f)

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
        # Perform the action requested by the user 
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
        if data['FILETYPE'] == 'G09_INPUT':
            if coords:
                atoms = data['ATOMS']
                coords = data['COORDS']
                bonds = data['BONDS']
                angles = data['ANGLES']
                dihedrals = data['DIHEDRALS']

                # Perform variable scan
                if dscan != []:
                    # Equilibrium dihedral
                    deq = dihedrals[pdata[var]]

                    # Build up lists of dihedrals
                    lscan = dscan[0]
                    hscan = dscan[1]+dscan[2]
                    step = dscan[2]
                    pd = []
                    md = []
                    for d in np.arange(deq, deq+hscan, step):
                        pd.append(d)
                    for d in np.arange(lscan, deq, step):
                        md.append(d)

                    d = md + pd

                    # Print the files in the positive direction
                    for num in d:
                        string = str(num).split('.')
                        if num < 0.0:
                            string[0] = string[0].replace('-', 'm')
                        fname = string[0] + 'pt' + string[1] + '_' + var.lower()
                        fname += '.g09'
                        fd = open(fname, 'w')
                        for item in data['JOBTYPE']:
                            print(item, file=fd)
                        print('', file=fd)
                        print(data['TITLE'], file=fd)
                        print('', file=fd)
                        print(data['CHARGE_MULT'], file=fd)
                        bline = '{0:2} {1:1} {2:8.6f}'
                        aline = '{0:2} {1:1} {2:8.6f} {3:1} {4:8.4f}'
                        dline = '{0:2} {1:1} {2:8.6f} {3:1} {4:8.4f} {5:1} {6:9.4f}'
                        for i in range(data['NATOMS']):
                            if i == 0:
                                print(atoms[i], file=fd)
                            elif i == 1:
                                print(bline.format(atoms[i], coords[i][0], bonds[i-1]), file=fd)
                            elif i == 2:
                                print(aline.format(atoms[i], coords[i][0], bonds[i-1],
                                                   coords[i][2], angles[i-2]), file=fd)
                            else:
                                if i-3 == pdata[var]:
                                    print(dline.format(atoms[i], coords[i][0], bonds[i-1],
                                                       coords[i][2], angles[i-2],
                                                       coords[i][4], num), file=fd)
                                else:
                                    print(dline.format(atoms[i], coords[i][0], bonds[i-1],
                                                       coords[i][2], angles[i-2],
                                                       coords[i][4], dihedrals[i-3]), file=fd)
                        print('', file=fd)
                        print('', file=fd)
                        fd.close()
        elif data['FILETYPE'] == 'G09_OUTPUT':
            # Print energies
            if energies:
                en = data['ENERGY']
                
                if plot:
                    if data['BONDS'] == []:
                        sys.exit('Must use file that has coordinates in Z-matrix format!')
                    else:
                        # For FitEVB, we need data in the form of data IDs, rather than
                        # proton transfer coordinate lengths.
                        if fitevb:
                            tmp = f.split('.')[0].upper()
                            x.append(tmp)
                        else:
                            x.append(data[xdata][pdata[var]])
                        # Determine what the y-axis quantity is
                        if 'PCM' in en.keys():
                            y.append(en['PCM']*conv[units][1])
                        elif 'MP2' in en.keys():
                            y.append(en['MP2']*conv[units][1])
                        else:
                            y.append(en['SCF']*conv[units][1])
                else:
                    # Print the data
                    fmt = '{0:>18} {1:20.8f}' 
                    
                    fname = data['FILENAME']
                    txt = 'Output data for file: ' + fname

                    print()
                    print('='*len(txt))
                    print(txt)
                    print('='*len(txt))

                    print()
                    print('Data is output in units of: ', conv[units][0])
                    print()
                    print(fmt.format(table['SCF'], data['ENERGY']['SCF']*conv[units][1]))
                    # Print post-HF energy
                    if 'MP2' in data['ENERGY'].keys():
                        print(fmt.format(table['MP2'], data['ENERGY']['MP2']*conv[units][1]))
                        print(fmt.format(table['CORR'], data['ENERGY']['CORR']*conv[units][1]))
                    # Print PCM energy
                    if 'PCM' in data['ENERGY'].keys():
                        print(fmt.format(table['PCM'], data['ENERGY']['PCM']*conv[units][1]))
                    print()
            # Print coordinates
            elif coords:
                fname = f.split('.')[0]
                output = fname + '.xyz' 
                foutput = open(output, 'w')

                fmt = '{0:4} {1:10.6f} {2:10.6f} {3:10.6f}'

                print(data['NATOMS'], file=foutput)
                print('', file=foutput)
                for i in range(data['NATOMS']):
                    print(fmt.format(data['ATOMS'][i], data['CART_COORDS'][i][0], 
                                     data['CART_COORDS'][i][1], data['CART_COORDS'][i][2]), file=foutput)

                foutput.close()
            # Print forces
            elif forces:
                fname = f.split('.')[0].upper()
                output = 'REF_' + fname 
                foutput = open(output, 'w')

                fmt = '{0:4} {1:12.6f} {2:12.6f} {3:12.6f}'

                # Conversion between Bohr and Angstroms
                au2angstrom = 0.52917720859

                # Conversion between Hartree and kcal/mol
                hart2kcalmol = 627.50947428

                for i in range(data['NATOMS']):
                    print(fmt.format( i+1, 
                                      data['FORCES'][i][0] * hart2kcalmol / au2angstrom,
                                      data['FORCES'][i][1] * hart2kcalmol / au2angstrom, 
                                      data['FORCES'][i][2] * hart2kcalmol / au2angstrom), file=foutput)

                foutput.close()
            ## Print coordinates from a transition state search
            #elif tscoords:
            #    output = 'tscoordinates.xyz'
            #    foutput = open(output, 'a')

            #    fmt = '{0:4} {1:10.6f} {2:10.6f} {3:10.6f}'


            #    print(data['NATOMS'], file=foutput)
            #    print('', file=foutput)
            #    for i in range(data['NATOMS']):
            #        print(fmt.format(data['ATOMS'][i], data['CART_COORDS'][i][0],
            #                         data['CART_COORDS'][i][1], data['CART_COORDS'][i][2]), file=foutput)

            #    foutput.close()
            # Print atomic dipoles
            elif atomdip:
                field = data['FIELD']
                ad = data['ATOMIC_DIPOLE']['HIRSHFELD']
                natoms = data['NATOMS']
                atoms = data['ATOMS']

                fmthead = '{0} {1:6.3f} {2:6.3f} {3:6.3f}'
                print()
                print(fmthead.format('External field (a.u.) =', field[0], field[1], field[2]))       
                print()

                fmtdip = '{0:2} {1:9.6f} {2:9.6f} {3:9.6f}'
                print('      X         Y         Z')
                for i in range(natoms):
                    print(fmtdip.format(atoms[i], ad[i][0], ad[i][1], ad[i][2]))
                print()
            # Store values for atomic polarizabilities 
            elif atompol:
                ff.append(data['FIELD'])
                ad.append(data['ATOMIC_DIPOLE']['HIRSHFELD'])
                natoms = data['NATOMS']
                atoms = data['ATOMS']
            elif charges:
                # Types of charges
                charge_types = { 'Merz-Singh-Kollman' : 'Merz-Kollman',
                                 'Mulliken'           : 'Mulliken',
                               }
                for item in data['CHARGE_MODEL'].keys():
                    ch_type = charge_types[item]
                    # Stash data in a convenient way
                    chg = data['CHARGE_MODEL'][item]
                    atoms = data['ATOMS']
                    # Heading
                    print()
                    print(ch_type + ' charges')
                    print()
                    fmtchg = '{0:2} {1:9.6f}'
                    for i in range(data['NATOMS']):
                        print(fmtchg.format(atoms[i], chg[i]))
                    print()
        else:
            sys.exit('File must be a G09 input or output file')

    if plot:
        # Plot data depending on the user's request.
        # Set to render text with LaTeX, edit font properties, define figure

        # Write data to file
        if fitevb:
            fh = open('ENERGY', 'w')
            fmt0 = '{0} {1}'
            fmt = '{0} {1:20.8f}'
            print(fmt0.format('ZERO', x[y.index(min(y))]), file=fh)
            for i in range(len(x)):
                print(fmt.format(x[i], y[i]), file=fh)
        else:
            fh = open('potential.data', 'w')
            fmt = '{0:12.6f} {1:20.8f}'
            for i in range(len(x)):
                print(fmt.format(x[i], y[i]), file=fh)

            # x-label and fitting information for bonds and angles
            boa = 'False'
            angle = 'False'
            if xdata == 'BONDS':
                xlabel = r'$r-r_{eq}$ (\AA)'
                boa = 'True'
            elif xdata == 'ANGLES':
                xlabel = r'$\theta-\theta_{eq}$ (Degrees)'
                boa = 'True'
                angle = 'True'
            elif xdata == 'DIHEDRALS':
                xlabel = r'$\phi$ (Degrees)'
    
            # Define units for y-axis
            yunit = conv[units][0] 

            string = dedent('''\
            import matplotlib.pyplot as plt
            from matplotlib.ticker import MultipleLocator
            import numpy as np
            from scipy import stats

            plt.rc('text', usetex=True)
            plt.rc('font', **{'family':'serif', 'serif': ['Times'], 'size': 30})

            #*********************************************************
            # Fill this in if you want values in terms of displacement
            #*********************************************************
            eq = 0.0000
            
            # Load the data
            fh = 'potential.data'
            x,y = np.loadtxt(fh,usecols=(0,1),unpack=True)
            x = x - eq
            y = y - min(y)

            if {boa}:
                x2 = np.linspace(min(x), max(x), 100)
                y2 = np.polyfit(x, y, 2)
                y3 = np.poly1d(y2, variable='x')
                # y3 is the formula of the parabola in the harmonic approach
                print(y3)
                if {angle}:
                    # Convert the force constant to radians and print it
                    deg2rad = np.pi/180.0
                    ktheta = y3.c[0]/(deg2rad*deg2rad)
                    print('{0:17} {1:7.4f}'.format('Force constant = ', ktheta))

                # Coefficients
                # x-axis and y-axis crossing points
                xmin = -y3.c[1]/(2.0*y3.c[0])
                ymin = y3.c[0]*xmin*xmin + y3.c[1]*xmin + y3.c[2] 
                print('{0:13} {1:7.4f}'.format('x-axis min = ', xmin))
                print('{0:13} {1:7.4f}'.format('y-axis min = ', ymin))

                # Determine the r-squared value for the fit
                y4 = []
                for i in range(len(x)):
                    y4.append(y3(x[i]))
                slope, intercept, r_value, p_value, std_err = stats.linregress(y,y4)
                print('{0:12} {1:6.4f}'.format('R-squared = ', r_value*r_value))

            # Make a plot
            fig = plt.figure()
            sub = fig.add_subplot(111)
            sub.plot(x, y, 'b', marker='o', linewidth=0)
            if {boa}:
                sub.plot(x2, y3(x2), 'g', linewidth=2)

            fig.subplots_adjust(left=0.12, right=0.90, bottom=0.1, top=0.9)

            # Title, labels
            sub.set_xlabel(r'{xlabel}')

            # Define units for y-axis
            lab = 'Energy (' + '{yunit}' + ')'
            sub.set_ylabel(lab)

            # Axis limits
            sub.set_xlim([round(min(x),3),round(max(x),3)])

            plt.show()\
            ''')
            string = string.replace('{xlabel}', xlabel)
            string = string.replace('{yunit}', yunit)
            string = string.replace('{boa}', boa)
            string = string.replace('{angle}', angle)

            dataplot = open('plot_potential.mpl.py', 'w')
            print(string, file=dataplot) 

    if atompol:
        # Because the systems are spherical, we only need to account for
        # the diagonal terms of the polarizability

        # Magnitude of the finite field
        mff = []
        for f in ff:
            tmp = f[0]*f[0] + f[1]*f[1] + f[2]*f[2]
            tmp = np.sqrt(tmp)
            mff.append(tmp)

        # Get the zero electric field atomic dipoles
        for i in range(len(mff)):
            lf = mff[i]
            if lf == 0.0:
                zfdip = ad[i]
                s = i + 1

        # Determine the polarizability
        ap = [[0.0, 0.0, 0.0] for i in range(natoms)]
        for i in range(s,len(mff)):
            for a in range(natoms):
                if ff[i][0] != 0.0:
                    ap[a][0] = ( ad[i][a][0] - zfdip[a][0] ) / ff[i][0]
                elif ff[i][1] != 0.0:
                    ap[a][1] = ( ad[i][a][1] - zfdip[a][1] ) / ff[i][1]
                elif ff[i][2] != 0.0:
                    ap[a][2] = ( ad[i][a][2] - zfdip[a][2] ) / ff[i][2]

        # Determine the polarizability average
        pol_avg = []
        for i in range(natoms):
            tmp = ap[i][0]*ap[i][0] + ap[i][1]*ap[i][1] + ap[i][2]*ap[i][2]
            tmp *= 1.0 / 3.0
            pol_avg.append(tmp)

        # Print average polarizabilities
        fmt = '{0:2} {1:6.2f}'
        print()
        print('Atomic polarizabilities in a.u.')
        print()
        for i in range(natoms):
            print(fmt.format(atoms[i], pol_avg[i]))
        print()

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
