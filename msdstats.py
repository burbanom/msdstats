#! /usr/bin/env python

# Contributors:
# Benjamin J. Morgan (b.j.morgan@bath.ac.uk)
# Mario Burbano (burbano.carmona@gmail.com)

from __future__ import print_function
import numpy as np
import pandas as pd
import argparse
from subprocess import check_output
import matplotlib.pyplot as plt
from os.path import isfile
import sys

def parse_commandline_arguments():
    parser = argparse.ArgumentParser( description = 'msd error analysis' )
    parser.add_argument( '--atoms', '-a', metavar = 'N', type=int, required = True, help='set the number of atoms' )
    parser.add_argument( '--frames', '-f', metavar = 'N', type=int, required = False, help='number of displacement frames in the complete displacement file' )
    ###########################################################################################################################################
    # The following two values correspond to tau and delta, respectively in the following article:
    # http://pubs.acs.org/doi/abs/10.1021/acs.jctc.5b00574?journalCode=jctcce
    parser.add_argument( '--slicesize', '-ss', metavar = 'N', type=int, required = True, help='Size of the slice in number of frames' )
    parser.add_argument( '--offset', '-so', metavar = 'N', type=int, required = True, help='Offset (in number of frames) between the slices' )
    ############################################################################################################################################
    parser.add_argument( '--msdlen', '-ml', metavar = 'N', type=int, required = True, help='Length of msd' )
    parser.add_argument( '--prntfrq', '-p', metavar = 'N', type=int, required = False, default = 1000, help='Print frequency of displacements' )
    parser.add_argument( '--executable', '-x', help='msd code executable name', default='msdconf.x' )
    parser.add_argument( '--displacement-file', '-d', help='complete displacement file filename', default='displong.out' )
    parser.add_argument( '--msd-files', '-m', nargs='+', help='list of msd output files to analyse', required = True )
    parser.add_argument( '--convcalc', '-cc', action='store_true', help='Perform a slope convergence calculation rather than the slope statistics calculation.' )
    return parser.parse_args()

def read_msd_template(template_file='msd_template'):
    from os.path import isfile
    lines = []
    try:
        isfile( template_file ) 
    except:
        sys.exit("File " + template_file + " not found")
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

def slope_statistics( disp_data_np, natoms, nframes, msd_length, nprint, slice_size, slice_offset, files_to_monitor, msd_executable ):

    my_template = read_msd_template('msd_template')
    write_msd_inpt(my_template,msd_length=msd_length,n_frames=slice_size, nprint=nprint)

    nslices = int( ( nframes - slice_size ) / slice_offset ) + 1
    temp_disp_filename = 'dispslice.out'

    complete_disp_data = disp_data_np 
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

def slope_convergence( disp_data_np, natoms, nframes, msd_length, nprint, slice_offset, files_to_monitor, msd_executable ):
# This function analyses how the slopes change as we increase the simulation time. It creates slices of the 
# displacement file that become increasingly larger. 
    my_template = read_msd_template('msd_template')
    
    temp_disp_filename = 'dispslice.out'

    #slices = range( 0, nframes - slice_offset + 1,  slice_offset )
    #nslices = int( ( nframes - slice_offset ) / slice_offset ) + 1
    slices = range( slice_offset, nframes + 1,  slice_offset )
    nslices = len( slices )
    complete_disp_data = disp_data_np 
    complete_disp_data = complete_disp_data.reshape( ( nframes, natoms, 3 ) )
    msd_slopes = np.empty( ( len( files_to_monitor ), nslices ) )
    
    for j, new_length in enumerate( slices ):
        
        sliced_disp_data = complete_disp_data[ 0 : new_length ]
        sliced_disp_data = sliced_disp_data.reshape( new_length * natoms, 3 )
        np.savetxt( temp_disp_filename, sliced_disp_data )
        write_msd_inpt(my_template,msd_length=msd_length,n_frames=new_length, nprint=nprint)
        out = check_output( [ msd_executable ] )
        for i, filename in enumerate( files_to_monitor ):
            msd_slopes[i,j] = get_slope_from_msd_output( filename )
    columns = slices
    slopes_df = pd.DataFrame( msd_slopes, index=files_to_monitor, columns=columns )
    slopes_df.to_csv('slope-conv.csv')
    return msd_slopes

def open_displacements( displacements_file, natoms ):
    import sys
    disp_data = np.loadtxt( displacements_file )
    disp_data_shape = np.shape( disp_data )
    if np.mod( disp_data_shape[0], natoms ) != 0:
        sys.exit('ERROR = The number of frames in the file is not an integer multiple of the number of atoms.')
    nframes = int( disp_data_shape[0] / natoms )
    return disp_data, nframes


if __name__ == '__main__':
    args = parse_commandline_arguments()

    natoms = args.atoms
    #nframes_tot = args.frames
    slice_size_frames = args.slicesize
    slice_offset = args.offset
    msd_executable = args.executable
    displacements_file = args.displacement_file # displacements without the header information
    msd_files = args.msd_files
    msd_length = args.msdlen
    convcalc = args.convcalc
    prntfrq = args.prntfrq
    
    try:
        isfile( displacements_file ) 
    except:
        sys.exit("File " + displacements_file + " not found")

    complete_disp_data, nframes_tot = open_displacements( displacements_file, natoms )

    if not convcalc:
        slope_stats = slope_statistics( complete_disp_data, natoms, nframes_tot, msd_length, prntfrq, slice_size_frames, slice_offset, msd_files, msd_executable ) 
    else:
        slope_conv = slope_convergence( complete_disp_data, natoms, nframes_tot, msd_length, prntfrq, slice_offset, msd_files, msd_executable )
