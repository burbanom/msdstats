#! /usr/bin/env python

# Contributors:
# Benjamin J. Morgan (b.j.morgan@bath.ac.uk)
# Mario Burbano (mario.burbano@upmc.fr)

from __future__ import print_function
import numpy as np
import pandas as pd
import argparse
from subprocess import check_output
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
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

def slope_statistics( disp_filename, natoms, nframes, slice_size, slice_offset, files_to_monitor, msd_executable ):

    nslices = int( ( nframes - slice_size ) / slice_offset ) + 1
    temp_disp_filename = 'dispslice.out'

    complete_disp_data = np.loadtxt( disp_filename )
    complete_disp_data = complete_disp_data.reshape( ( nframes, natoms, 3 ) )

    msd_slopes = np.empty( ( len( files_to_monitor ), nslices ) )

    for j, initial_frame in enumerate( range( 0, nframes - slice_size + 1, slice_offset ) ):
        final_frame = initial_frame + slice_size
        sliced_disp_data = complete_disp_data[ initial_frame : final_frame ]
        sliced_disp_data = sliced_disp_data.reshape( slice_size * natoms, 3 )
        np.savetxt( temp_disp_filename, sliced_disp_data )
        out = check_output( [ msd_executable ] )
        for i, filename in enumerate( msd_files ):
            msd_slopes[i,j] = get_slope_from_msd_output( filename )

	columns = ['Mean','Standard Deviation','sim length','tau','delta','number of slices']
	stats = pd.DataFrame( columns = columns, index= msd_files )
	for i, filename in enumerate( msd_files ):
	    stats.loc[filename] =[ np.mean(msd_slopes[i]), np.std(msd_slopes[i]), nframes, slice_size, slice_offset, nslices ]
	stats.to_csv('slope-stats.csv')

    return msd_slopes


if __name__ == '__main__':
    args = parse_commandline_arguments()

    natoms = args.atoms
    nframes_tot = args.frames
    slice_size_frames = args.slicesize
    slice_offset = args.offset
    msd_executable = args.executable
    displacements_file = args.displacement_file # displacements without the header information
    msd_files = args.msd_files

    msd_slopes = slope_statistics( displacements_file, natoms, nframes_tot, slice_size_frames, slice_offset, msd_files, msd_executable ) 

