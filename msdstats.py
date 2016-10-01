#! /usr/bin/env python

# Contributors:
# Benjamin J. Morgan (b.j.morgan@bath.ac.uk)
# Mario Burbano (mario.burbano@upmc.fr)

from __future__ import print_function
import numpy as np
import pandas as pd
import argparse
from subprocess import check_output
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

def read_msd_template(template_file='msd_template'):
    lines = []
    with open(template_file) as f:
        for line in f.readlines():
            lines.append(line)
    return lines

def write_msd_inpt(template,msd_length,n_frames,nprint=1000,timestep=41.3414,dist_warn=3.5):
    f = open('msd.inpt','w')
    for line in template:
        f.write(line)
    f.write(str(msd_length) + '\n')
    f.write('.false.\n')
    f.write('.true.\n')
    f.write(str(nprint) + '\n')
    f.write(str(timestep) + '\n')
    f.write(str(dist_warn) + '\n')
    f.write(str(n_frames*nprint) + '\n')
    f.close()

def plot_slopes(msd_slopes,files_list,tot_sim_length,delta,tau):
    plt.figure()
    plt.title(r'\delta = ' + str(delta) + '\n' + r'\tau = ' + str(tau))
    for i, filename in enumerate( files_list ):
        plt.plot(msd_slopes[i],'o-',label=filename)
        plt.legend(loc='best')

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

def slope_convergence( disp_filename, natoms, nframes, msd_length, slice_offset, files_to_monitor, msd_executable ):
# This function analyses how the slopes change as we increase the simulation time. It creates slices of the 
# displacement file that become increasingly larger. 
    my_template = read_msd_template('msd_template')
    
    temp_disp_filename = 'dispslice.out'

    complete_disp_data = np.loadtxt( disp_filename )
    complete_disp_data = complete_disp_data.reshape( ( nframes, natoms, 3 ) )
    slices = range( slice_offset, nframes + slice_offset, slice_offset )
    nslices = len( slices )
    msd_slopes = np.empty( ( len( files_to_monitor ), nslices ) )
    
    for j, new_length in enumerate( slices ):
        
        sliced_disp_data = complete_disp_data[ 0 : new_length ]
        sliced_disp_data = sliced_disp_data.reshape( new_length * natoms, 3 )
        np.savetxt( temp_disp_filename, sliced_disp_data )
        write_msd_inpt(my_template,msd_length=msd_length,n_frames=new_length)
        out = check_output( [ msd_executable ] )
        for i, filename in enumerate( files_to_monitor ):
            msd_slopes[i,j] = get_slope_from_msd_output( filename )
    columns = slices
    slopes_df = pd.DataFrame( msd_slopes, index=files_to_monitor, columns=columns )
    slopes_df.to_csv('slope-conv.csv')
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
    
    slope_conv = slope_convergence( displacements_file, natoms, nframes_tot, slice_offset, slice_offset, msd_files, msd_executable )
    slope_stats = slope_statistics( displacements_file, natoms, nframes_tot, slice_size_frames, slice_offset, msd_files, msd_executable ) 

