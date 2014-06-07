#! /usr/bin/env python

from __future__ import division, print_function
import sys
import numpy as np
from collect import collect

def main():
    """\
    Script for reading a CP2K output (*.out).
    """

    from argparse import ArgumentParser, RawDescriptionHelpFormatter
    from textwrap import dedent
    parser = ArgumentParser(description=dedent(main.__doc__),
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('-e', '--energies', help='Output energies.',
                        action='store_true', default=False)
    parser.add_argument('-p', '--plot', help='Use this to make a plot of data.',
                        action='store_true', default=False)
    parser.add_argument('-u', '--units', help='Units for outputting quantities.',
                        default='au')
    parser.add_argument('-b', '--bondatoms', help='The atoms in the bond for a PT scan.'
                        ' These should be written in a "# #" pair.', required=False)
    parser.add_argument('files', help='The files that you want analyzed.'
                        , nargs='+')
    args = parser.parse_args()

    # Store command line arguments as convenient variables
    fh = args.files
    energy = args.energies
    plot = args.plot
    units = args.units
    if args.bondatoms is not None:
        tmp = args.bondatoms.split()
        bond = [int(tmp[0]), int(tmp[1])]
    else:
        bond = None
    if plot and bond is None:
        sys.exit('Need to specify the plot and bond flags together.')

    # Unit conversion dictionary
    conv = { 'au'      : ['Hartree' ,      1.00000000],
             'ev'      : ['eV'      ,     27.21138505],
             'kcalmol' : ['kcal/mol',    627.50947428],
             'kjmol'   : ['kJ/mol'  ,   2625.49964038],
             'wvn'     : ['cm^{-1}' , 219474.63068   ],
           }
    
    # /\/\/\/\/\/\/\/\/
    # Collect the files 
    # /\/\/\/\/\/\/\/\/

    if energy:
        e = None
        e_qmmm = None
        if plot:
            x = []
            y = []

    for f in fh:
        # Collect all data into memory for data retention
        data = collect(f)

        # Output data
        if energy:
            if plot:
                # This is currently designed to work for proton transfer
                # scans.

                # Get the bond distance.
                fxyz = f.split('.')[0] + '.xyz'
                xyz = collect(fxyz)
                coords1 = xyz['COORDS'][bond[0]-1]
                coords2 = xyz['COORDS'][bond[1]-1]
                d = np.sqrt( ( coords2[0] - coords1[0] ) ** 2
                           + ( coords2[1] - coords1[1] ) ** 2 
                           + ( coords2[2] - coords1[2] ) ** 2 )
                x.append(d)

                # Store the energies
                if 'QMMM_TOTAL_ENERGY' in data['ENERGIES'].keys():
                    if e is None:
                        e = []
                        e.append( data['ENERGIES']['QM_TOTAL_ENERGY']*conv[units][1] )
                        e_qmmm = []
                        e_qmmm.append( data['ENERGIES']['QMMM_TOTAL_ENERGY']*conv[units][1] )
                    else:
                        e.append( data['ENERGIES']['QM_TOTAL_ENERGY']*conv[units][1] )
                        e_qmmm.append( data['ENERGIES']['QMMM_TOTAL_ENERGY']*conv[units][1] )
                else:
                    if e is None:
                        e = []
                        e.append( [ data['ENERGIES']['QM_TOTAL_ENERGY']*conv[units][1], 0.0 ] )
                    else:
                        e.append( [ data['ENERGIES']['QM_TOTAL_ENERGY']*conv[units][1], 0.0 ] )
            else:
                # Table of energy term meanings
                table = { 'OVERLAP_ENERGY'            : 'Overlap energy of the core charge distribution:',
                          'CORE_SELF_ENERGY'          : 'Self energy of the core charge distribution:',
                          'CORE_HAMILTONIAN'          : 'Core Hamiltonian energy:',
                          'COULOMB_ENERGY'            : 'Coulomb energy:',
                          'DFT_XC_ENERGY'             : 'Exchange-correlation energy:',
                          'HFX_ENERGY'                : 'Hartree-Fock Exchange energy:',
                          'DISPERSION_ENERGY'         : 'Dispersion energy:',
                          'QMMM_ELECTROSTATIC_ENERGY' : 'QM/MM Electrostatic energy:',
                          'QM_TOTAL_ENERGY'           : 'Total QM energy:',
                          'QMMM_TOTAL_ENERGY'         : 'Total QM/MM energy:',
                        }                     

                # Print the data
                fmt = '{0:>47} {1:20.8f}'

                fname = data['FILENAME']
                txt = 'Output data for file: ' + fname

                print()
                print('='*len(txt))
                print(txt)
                print('='*len(txt))

                print()
                print('Data is output in units of: ', conv[units][0])
                print()
                print(fmt.format(table['OVERLAP_ENERGY'], data['ENERGIES']['OVERLAP_ENERGY']*conv[units][1]))
                print(fmt.format(table['CORE_SELF_ENERGY'], data['ENERGIES']['CORE_SELF_ENERGY']*conv[units][1]))
                print(fmt.format(table['CORE_HAMILTONIAN'], data['ENERGIES']['CORE_HAMILTONIAN']*conv[units][1]))
                print(fmt.format(table['COULOMB_ENERGY'], data['ENERGIES']['COULOMB_ENERGY']*conv[units][1]))
                print(fmt.format(table['DFT_XC_ENERGY'], data['ENERGIES']['DFT_XC_ENERGY']*conv[units][1]))
                if 'HFX_ENERGY' in data['ENERGIES'].keys():
                    print(fmt.format(table['HFX_ENERGY'], data['ENERGIES']['HFX_ENERGY']*conv[units][1]))
                if 'DISPERSION_ENERGY' in data['ENERGIES'].keys():
                    print(fmt.format(table['DISPERSION_ENERGY'], data['ENERGIES']['DISPERSION_ENERGY']*conv[units][1]))
                print(fmt.format(table['QMMM_ELECTROSTATIC_ENERGY'], data['ENERGIES']['QMMM_ELECTROSTATIC_ENERGY']*conv[units][1]))
                print(fmt.format(table['QM_TOTAL_ENERGY'], data['ENERGIES']['QM_TOTAL_ENERGY']*conv[units][1]))
                print(fmt.format(table['QMMM_TOTAL_ENERGY'], data['ENERGIES']['QMMM_TOTAL_ENERGY']*conv[units][1]))
                print()

    if energy and plot:
        # Make a file containing the data
        fmt = '{0:12.8f} {1:12.8f} {2:12.8f}'
        fname = 'ptcoord.data'
        f = open(fname, 'w')
        for i in range(len(x)):
            if e_qmmm is not None:
                print(fmt.format( x[i], e[i], e_qmmm[i] ), file=f)
            else:
                print(fmt.format( x[i], e[i], 0.0 ), file=f)
        # For plotting the data
        # X-axis label 
        xlabel = r'$R$ (\AA)'
        # Define units for y-axis
        yunit = conv[units][0]

        string = dedent('''\
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MultipleLocator
        import numpy as np
        from scipy import stats

        plt.rc('text', usetex=True)
        plt.rc('font', **{'family':'serif', 'serif': ['Times'], 'size': 30})

        # Load the data
        fh = 'ptcoord.data'
        x,y = np.loadtxt(fh,usecols=(0,1),unpack=True)
        y = y - min(y)
        x2,y2 = np.loadtxt(fh,usecols=(0,2),unpack=True)
        y2 = y2 - min(y2)

        # Make a plot
        fig = plt.figure()
        sub = fig.add_subplot(111)
        sub.plot(x, y, 'b', marker='o', linewidth=1)
        sub.plot(x2, y2, 'g', marker='o', linewidth=1)

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

        dataplot = open('plot_potential.mpl.py', 'w')
        print(string, file=dataplot)
        
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
