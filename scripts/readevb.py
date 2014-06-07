#! /usr/bin/env python

from __future__ import division, print_function
import sys
import numpy as np
from numpy import array, append
import matplotlib.pyplot as plt
from collect import collect

def main():
    """\
    Script for reading a RAPTOR output (evb.out).  This is meant to allow
    for easier analysis of data.
    """

    from argparse import ArgumentParser, RawDescriptionHelpFormatter
    from textwrap import dedent
    parser = ArgumentParser(description=dedent(main.__doc__),
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('-r', '--rxncenter', help='Output the reaction center.',
                        action='store_true', default=False)    
    parser.add_argument('-e', '--energy', help='Output the total environmental energy,'
                        ' complex energy, and total energy (env. + complex).',
                        action='store_true', default=False)    
    parser.add_argument('-eed', '--envenergydecomp', 
                        help='Output the decomposed environmental energy.',
                        action='store_true', default=False)    
    parser.add_argument('-d', '--diagonal', help='Output the complex energy diagonal '
                        'for the reaction center.',
                        action='store_true', default=False)    
    parser.add_argument('-od', '--offdiagonal', help='Output the complex energy off-diagonal ' 
                        'for the reaction center.', 
                        action='store_true', default=False)    
    parser.add_argument('-ev', '--eigenvector', help='Output the complex eigenvector.',
                        action='store_true', default=False)    
    parser.add_argument('-cec', '--ceccoord', help='Output the CEC coordinate.',
                        action='store_true', default=False)    
    parser.add_argument('-p', '--plot', help='Plot the data requested by the user.',
                        action='store_true', default=False)
    parser.add_argument('files', help='The files that you want converted.'
                        , nargs='+')
    args = parser.parse_args()

    # Store command line arguments as convenient variables
    fh = args.files[0]
    rxn_center = args.rxncenter
    energy = args.energy
    decomp_env = args.envenergydecomp
    diagonal = args.diagonal
    offdiagonal = args.offdiagonal
    civec = args.eigenvector
    cec = args.ceccoord
    plot = args.plot

    test = ( rxn_center is False and
             energy is False and
             decomp_env is False and
             diagonal is False and
             offdiagonal is False and
             civec is False and
             cec is False )
    if test:
        msg = "Please request output of specific data.  Use -h to see what is available"
        error = '\n' + len(msg)*'%' + '\n' + msg + '\n' + len(msg)*'%' + '\n'
        sys.exit(error)

    # /\/\/\/\/\/\/\/\/\/\/\/\/\/
    # Collect EVB simulation data 
    # /\/\/\/\/\/\/\/\/\/\/\/\/\/

    evb_data = collect(fh)

    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    # Output data requested by the user 
    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    # Plot the data
    if plot:
        plt.rc('text', usetex=True)
        plt.rc('font', **{'family':'serif', 'serif': ['Times'], 'size': 30})
        fig = plt.figure()
        if rxn_center:
            x = array(evb_data['TIMESTEP'])/1000.0
            y = [i[1] for i in evb_data['RXNCENTER']]
            sub = fig.add_subplot(111)
            sub.plot(x, y, 'b', linewidth=3)

            # Adjust graph scale
            fig.subplots_adjust(left=0.12, right=0.90, bottom=0.1, top=0.9)

            # Title, labels
            sub.set_xlabel(r'time (ps)')
            sub.set_ylabel(r'Reaction Center')

            # Axis limits
            sub.set_xlim([0,x[-1]])
            sub.set_ylim([0,y[-1]+10]) 

            plt.show()
        elif energy:
            x = array(evb_data['TIMESTEP'])/1000.0
            env = evb_data['ENERGY_ENV']
            comp = evb_data['ENERGY_COMPLEX']
            total = evb_data['ENERGY_TOTAL']

            sub = fig.add_subplot(111)
            #sub.plot(x, env, 'b', linewidth=3)
            #sub.plot(x, comp, 'g', linewidth=3)
            sub.plot(x, total, 'r', linewidth=3)

            # Adjust graph scale
            fig.subplots_adjust(left=0.12, right=0.90, bottom=0.1, top=0.9)

            # Title, labels
            sub.set_xlabel(r'time (ps)')
            sub.set_ylabel(r'Potential Energy (kcal/mol)')

            # Axis limits
            sub.set_xlim([0,x[-1]])
            #sub.set_ylim([min(total)+0.20*min(total),0])
            sub.set_ylim([min(total)+0.01*min(total),max(total)-0.01*max(total)])

            plt.show()
        elif decomp_env:
            msg = "Unclear what to plot for environment decomposed energy."
            error = '\n' + len(msg)*'%' + '\n' + msg + '\n' + len(msg)*'%' + '\n'
            sys.exit(msg)
        elif diagonal:
            x = array(evb_data['TIMESTEP'])/1000.0
            y = [i[0][1] for i in evb_data['ENERGY_DIAGONAL']]

            sub = fig.add_subplot(111)
            sub.plot(x, y, 'b', linewidth=3)

            # Adjust graph scale
            fig.subplots_adjust(left=0.12, right=0.90, bottom=0.1, top=0.9)

            # Title, labels
            sub.set_xlabel(r'time (ps)')
            sub.set_ylabel(r'Diagonal Potential Energy Component (kcal/mol)')

            # Axis limits
            sub.set_xlim([0,x[-1]])
            sub.set_ylim([min(y)+0.2*min(y),0])

            plt.show()
        elif offdiagonal:
            x = array(evb_data['TIMESTEP'])/1000.0
            y = [i[0][4] for i in evb_data['ENERGY_OFF_DIAGONAL']]

            sub = fig.add_subplot(111)
            sub.plot(x, y, 'b', linewidth=3)

            # Adjust graph scale
            fig.subplots_adjust(left=0.12, right=0.90, bottom=0.1, top=0.9)

            # Title, labels
            sub.set_xlabel(r'time (ps)')
            sub.set_ylabel(r'Off-Diagonal Potential Energy Component (kcal/mol)')

            # Axis limits
            sub.set_xlim([0,x[-1]])
            sub.set_ylim([min(y)+0.20*min(y),0])

            plt.show()
        elif civec:
            x = array(evb_data['TIMESTEP'])/1000.0
            tmp = [i[1] for i in evb_data['CI_VECTOR']]
            
            y1 = []
            y2 = []
            for l in tmp:
                y1.append(max(l))
                y2.append(second_largest(l))

            sub = fig.add_subplot(111)
            sub.plot(x, y1, 'b', linewidth=3)
            sub.plot(x, y2, 'r', linewidth=3)

            # Adjust graph scale
            fig.subplots_adjust(left=0.12, right=0.90, bottom=0.1, top=0.9)

            # Title, labels
            sub.set_xlabel(r'time (ps)')
            sub.set_ylabel(r'MS-EVB state coefficient, c$_i$')

            # Axis limits
            sub.set_xlim([0,x[-1]])
            sub.set_ylim([0,1.0])

            plt.show()
        elif cec:
            msg = "Unclear what to plot for the center of excess charge."
            error = '\n' + len(msg)*'%' + '\n' + msg + '\n' + len(msg)*'%' + '\n'
            sys.exit(msg)
        else:
            # Code should never get here, but this is here just in case
            msg = "Cannot plot the data you requested.  Use -h to see available options."
            error = '\n' + len(msg)*'%' + '\n' + msg + '\n' + len(msg)*'%' + '\n'
            sys.exit(msg)
    # Output table to screen
    else:
        if rxn_center:
            fmt = '{0:10.1f} {1:2d} {2:5d} {3:3d}'
            for i in range(len(evb_data['TIMESTEP'])):
                print(fmt.format( evb_data['TIMESTEP'][i],
                                  evb_data['RXNCENTER'][i][0],
                                  evb_data['RXNCENTER'][i][1],
                                  evb_data['RXNCENTER'][i][2] ))
        elif energy:
            fmt = '{0:10.1f} {1:10.3f} {2:10.3f} {3:10.3f}'
            for i in range(len(evb_data['TIMESTEP'])):
                print(fmt.format( evb_data['TIMESTEP'][i],
                                  evb_data['ENERGY_ENV'][i],
                                  evb_data['ENERGY_COMPLEX'][i],
                                  evb_data['ENERGY_TOTAL'][i] ))
        elif decomp_env:
            fmt = '{0:10.1f} {1:11.3f} {2:11.3f} {3:11.3f}'
            fmt += ' {4:11.3f} {5:11.3f} {6:11.3f} {7:11.3f}'
            for i in range(len(evb_data['TIMESTEP'])):
                print(fmt.format( evb_data['TIMESTEP'][i],
                                  evb_data['ENV_ENERGY_DECOMP'][i][0],
                                  evb_data['ENV_ENERGY_DECOMP'][i][1],
                                  evb_data['ENV_ENERGY_DECOMP'][i][2],
                                  evb_data['ENV_ENERGY_DECOMP'][i][3],
                                  evb_data['ENV_ENERGY_DECOMP'][i][4],
                                  evb_data['ENV_ENERGY_DECOMP'][i][5],
                                  evb_data['ENV_ENERGY_DECOMP'][i][6] ))
        elif diagonal:
            fmt = '{0:10.1f} {1:5d} {2:11.3f} {3:11.3f} {4:11.3f}'
            fmt += ' {5:11.3f} {6:11.3f} {7:11.3f} {8:11.3f}'
            fmt += ' {9:11.3f}'
            for i in range(len(evb_data['TIMESTEP'])):
                print(fmt.format( evb_data['TIMESTEP'][i],
                                  evb_data['RXNCENTER'][i][1],
                                  evb_data['ENERGY_DIAGONAL'][i][0][1],
                                  evb_data['ENERGY_DIAGONAL'][i][0][2],
                                  evb_data['ENERGY_DIAGONAL'][i][0][3],
                                  evb_data['ENERGY_DIAGONAL'][i][0][4],
                                  evb_data['ENERGY_DIAGONAL'][i][0][5],
                                  evb_data['ENERGY_DIAGONAL'][i][0][6],
                                  evb_data['ENERGY_DIAGONAL'][i][0][7],
                                  evb_data['ENERGY_DIAGONAL'][i][0][8] ))
        elif offdiagonal:
            fmt = '{0:10.1f} {1:5d} {2:11.3f} {3:11.3f} {4:11.3f}'
            fmt += ' {5:11.3f}'
            for i in range(len(evb_data['TIMESTEP'])):
                print(fmt.format( evb_data['TIMESTEP'][i],
                                  evb_data['RXNCENTER'][i][1],
                                  evb_data['ENERGY_OFF_DIAGONAL'][i][0][1],
                                  evb_data['ENERGY_OFF_DIAGONAL'][i][0][2],
                                  evb_data['ENERGY_OFF_DIAGONAL'][i][0][3],
                                  evb_data['ENERGY_OFF_DIAGONAL'][i][0][4] ))
        elif civec:
            # Things are done differently here to adapt to the variable number
            # of MS-EVB states.  Recall that we only store "important" states,
            # with coefficients larger than 0.001.
            for i in range(len(evb_data['TIMESTEP'])):
                fmt = '{0:10.1f} '
                nmsevb = evb_data['CI_VECTOR'][i][0]
                for j in range(nmsevb):
                    fmt += '{1[' + str(j) + ']:8.5f} '
                print(fmt.format( evb_data['TIMESTEP'][i],
                                  evb_data['CI_VECTOR'][i][1] ))
        elif cec:
            fmt = '{0:10.1f} {1:5d} {2:10.6f} {3:10.6f} {4:10.6f}'
            for i in range(len(evb_data['TIMESTEP'])):
                print(fmt.format( evb_data['TIMESTEP'][i],
                                  evb_data['RXNCENTER'][i][1],
                                  evb_data['CEC'][i][0],
                                  evb_data['CEC'][i][1],
                                  evb_data['CEC'][i][2] ))
        else:
            # Code should never get here, but this is here just in case
            msg = "Cannot output the data you requested.  Use -h to see available options."
            error = '\n' + len(msg)*'%' + '\n' + msg + '\n' + len(msg)*'%' + '\n'
            sys.exit(msg)

def second_largest(numbers):
    '''Function for determining the second largest element in an array.'''
    m1, m2 = None, None
    for x in numbers:
        if x >= m1:
            m1, m2 = x, m1
        elif x > m2:
            m2 = x
    return m2
    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
