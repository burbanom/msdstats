#! /usr/bin/env python

# Contributors:
# Benjamin J. Morgan (b.j.morgan@bath.ac.uk)
# Mario Burbano (mario.burbano@upmc.fr)

from __future__ import print_function
import numpy as np
import pandas as pd
import argparse
from subprocess import check_output
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

def parse_commandline_arguments():
    parser = argparse.ArgumentParser( description = 'msd error analysis' )
    parser.add_argument( '--atoms', '-a', metavar = 'N', type=int, required = True, help='set the number of atoms' )
    parser.add_argument( '--frames', '-f', metavar = 'N', type=int, required = True, help='number of displacement frames in the complete displacement file' )
    ###########################################################################################################################################
    # The following two values correspond to tau and delta, respectively in the following article:
    # http://pubs.acs.org/doi/abs/10.1021/acs.jctc.5b00574?journalCode=jctcce
    parser.add_argument( '--slicesize', '-ss', metavar = 'N', type=int, required = True, help='Size of the slice in number of frames' )
    parser.add_argument( '--offset', '-so', metavar = 'N', type=int, required = True, help='Offset (in number of frames) between the slices' )
    ############################################################################################################################################
    parser.add_argument( '--executable', '-x', help='msd code executable name', default='msdconf.x' )
    parser.add_argument( '--displacement-file', '-d', help='complete displacement file filename', default='displong.out' )
    parser.add_argument( '--msd-files', '-m', nargs='+', help='list of msd output files to analyse', required = True )
    return parser.parse_args()

def linear_regression( x, y ):
    # http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.lstsq.html#numpy.linalg.lstsq
    A = np.vstack( [ x, np.ones( len(x) ) ] ).T
    m, c = np.linalg.lstsq( A, y )[0]
    return m, c

def get_slope_from_msd_output( filename ):
    msd_data = np.loadtxt( filename )
    msd_length = msd_data.shape[0]
    long_time_msd_data = msd_data[ msd_length / 5 : -1 ]
    np.savetxt( 'msd_test.tmp', long_time_msd_data )
    slope, intercept = linear_regression( long_time_msd_data[:,0], long_time_msd_data[:,1] )
    return slope

if __name__ == '__main__':
    args = parse_commandline_arguments()

    n_total_atoms = args.atoms
    n_total_disp_frames = args.frames
    n_slicesize_frames = args.slicesize
    delta_slicesize_start = args.offset
    msd_executable = args.executable
    long_disp_filename = args.displacement_file # displacements without the header information
    msd_files = args.msd_files
    number_of_slices = int( ( n_total_disp_frames - n_slicesize_frames ) / delta_slicesize_start ) + 1
    temp_disp_filename = 'dispslice.out'

    complete_disp_data = np.loadtxt( long_disp_filename )
    complete_disp_data = complete_disp_data.reshape( ( n_total_disp_frames, n_total_atoms, 3 ) )

    msd_slope = np.empty( ( len( msd_files ), number_of_slices ) )
    for j, initial_frame in enumerate( range( 0, n_total_disp_frames - n_slicesize_frames + 1, delta_slicesize_start ) ):
        final_frame = initial_frame + n_slicesize_frames
        sliced_disp_data = complete_disp_data[ initial_frame : final_frame ]
        sliced_disp_data = sliced_disp_data.reshape( n_slicesize_frames * n_total_atoms, 3 )
        np.savetxt( temp_disp_filename, sliced_disp_data )
        out = check_output( [ msd_executable ] )
        for i, filename in enumerate( msd_files ):
            msd_slope[i,j] = get_slope_from_msd_output( filename )

    columns = ['Mean','Standard Deviation','sim length','tau','delta','number of slices']
    stats = pd.DataFrame(columns = columns)
    #stats_file = open('stats.dat','w')
    #stats_file.write( '# filename\n' )
    #stats_file.write( '# fileindex mean std nframes\n' )
    for i, filename in enumerate( msd_files ):
	stats.loc[filename] =[ np.mean(msd_slope[i]), np.std(msd_slope[i]), n_total_disp_frames, n_slicesize_frames, delta_slicesize_start, number_of_slices ]
        #stats_file.write( '{}\n'.format( '# '+filename ))
        #stats_file.write( '{} {:f} {:f} {:f}\n'.format(i , np.mean( msd_slope[i] ), np.std( msd_slope[i] ), n_total_disp_frames ))
        #plt.plot(msd_slope[i],'o-',label=filename)
        #plt.legend(loc='best')
    stats.to_csv('stats.csv')
    #stats_file.close()
    #np.savetxt( 'slopes.dat', msd_slope )

    #plt.savefig('slopes.pdf')

