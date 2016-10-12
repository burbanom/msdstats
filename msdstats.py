#! /usr/bin/env python

# Contributors:
# Benjamin J. Morgan (b.j.morgan@bath.ac.uk)
# Mario Burbano (burbano.carmona@gmail.com)

from __future__ import print_function, division
import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
from os.path import isfile
import sys

"""
Program for calculating statistics of msd/diffusion data
"""
__author__ = "Mario Burbano, Benjamin J. Morgan"
__credits__ = ["Mario Burbano", "Benjamin J. Morgan"]
__maintainer__ = "Mario Burbano"
__license__ = "GPL"
__version__ = "1.0"
__email__ = "burbano.carmona@gmail.com"
__date__ = "11 Oct., 2016"
__status__ = "Production"


def parse_commandline_arguments():
    parser = argparse.ArgumentParser( description = 'msd error analysis. Requires a file system.in file with the number of atoms and their charges' )
    ###########################################################################################################################################
    # The following two values correspond to tau and delta, respectively in the following article:
    # http://pubs.acs.org/doi/abs/10.1021/acs.jctc.5b00574?journalCode=jctcce
    parser.add_argument( '--slicesize', '-ss', metavar = 'N', type=int, required = False, help='Size of the slice in number of frames' )
    parser.add_argument( '--offset', '-so', metavar = 'N', type=int, required = True, help='Offset (in number of frames) between the slices or increment in conv. calc.' )
    ############################################################################################################################################
    parser.add_argument( '--msdlen', '-ml', metavar = 'N', type=int, required = True, help='Length of msd' )
    parser.add_argument( '--prntfrq', '-p', metavar = 'N', type=int, required = False, default = 1000, help='Print frequency of displacements' )
    parser.add_argument( '--timestep', '-t', metavar = 'F', type=float, required = False, default = 41.3414, help='Simulation timestep in a.u.' )
    parser.add_argument( '--displacement-file', '-d', help='complete displacement file filename', default='displong.out.gz' )
    parser.add_argument( '--msd-files', '-m', nargs='+', help='list of msd output files to analyse', required = True )
    parser.add_argument( '--convcalc', '-cc', action='store_true', help='Perform a slope convergence calculation rather than the slope statistics calculation.' )
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
    long_time_msd_data = msd_data[ int( msd_length / 5 ) : -1 ]
    #np.savetxt( 'msd_test.tmp', long_time_msd_data )
    slope, intercept = linear_regression( long_time_msd_data[:,0], long_time_msd_data[:,1] )
    return slope

def slope_statistics( displacements, nspcs, charges, nframes, msd_length, nprint, timestep, slice_size, slice_offset, files_to_monitor ):

    import calcmsds

    natoms = np.sum( nspcs )
    nspecies = len(nspcs)

    nslices = int( ( nframes - slice_size ) / slice_offset ) + 1

    msd_slopes = np.empty( ( len( files_to_monitor ), nslices ) )

    for j, initial_frame in enumerate( range( 0, nframes - slice_size + 1, slice_offset ) ):
        final_frame = initial_frame + slice_size
        sliced_disp_data = displacements[ initial_frame * natoms : final_frame * natoms ]
        calcmsds.calcmsds.z = charges
        calcmsds.calcmsds.numspc = nspcs 
        calcmsds.calcmsds.xdisp_long = sliced_disp_data['x'].values
        calcmsds.calcmsds.ydisp_long = sliced_disp_data['y'].values
        calcmsds.calcmsds.zdisp_long = sliced_disp_data['z'].values
        this_calc = calcmsds.calcmsds.msdconf(num=natoms,nspecies=nspecies,
                nmsdlength=msd_length,nmsdcalltime=nprint,dtime=timestep,nrun= slice_size * nprint)

        for i, filename in enumerate( msd_files ):
            msd_slopes[i,j] = get_slope_from_msd_output( filename )

    columns = ['Mean','Standard Deviation','sim length','tau','delta','number of slices']
    stats = pd.DataFrame( columns = columns, index= msd_files )
    for i, filename in enumerate( msd_files ):
        stats.loc[filename] =[ np.mean(msd_slopes[i]), np.std(msd_slopes[i]), nframes, slice_size, slice_offset, nslices ]
    stats.to_csv('slope-stats.csv')

    return msd_slopes

def slope_convergence( displacements, nspcs, charges, nframes, msd_length, nprint, timestep, slice_offset, files_to_monitor ):
# This function analyses how the slopes change as we increase the simulation time. It creates slices of the 
# displacement file that become increasingly larger. 
    
    import calcmsds

    natoms = np.sum( nspcs )
    nspecies = len( nspcs )

    slices = range( slice_offset, nframes + 1,  slice_offset )
    nslices = len( slices )
    msd_slopes = np.empty( ( len( files_to_monitor ), nslices ) )
    
    for j, new_length in enumerate( slices ):
        
        sliced_disp_data = displacements[ 0 : new_length * natoms ]
        calcmsds.calcmsds.z = charges
        calcmsds.calcmsds.numspc = nspcs 
        calcmsds.calcmsds.xdisp_long = sliced_disp_data['x'].values
        calcmsds.calcmsds.ydisp_long = sliced_disp_data['y'].values
        calcmsds.calcmsds.zdisp_long = sliced_disp_data['z'].values
        this_calc = calcmsds.calcmsds.msdconf(num=natoms,nspecies=nspecies,
                nmsdlength=msd_length,nmsdcalltime=nprint,dtime=timestep,nrun= new_length * nprint)
        for i, filename in enumerate( files_to_monitor ):
            msd_slopes[i,j] = get_slope_from_msd_output( filename )
    columns = slices
    slopes_df = pd.DataFrame( msd_slopes, index=files_to_monitor, columns=columns )
    slopes_df.to_csv('slope-conv.csv')
    return msd_slopes

def disp_slicer( dispfile, tot_rows, row_start, row_end ):
    # This function prevents us from having to open the entire 
    # displacements file.
    header = len( [i for i in range(0, row_start)] )
    footer = len( [i for i in range(row_end, tot_rows) ] )
    return np.genfromtxt(dispfile, skip_header = header, skip_footer = footer, dtype= np.float64)

def my_system( infile = 'system.in' ):
    # This function reads the system information
    # from an input file, whose format should be
    # two comma-separated columns, the first with
    # the number of atoms of each species and the
    # second with their respective charges.
    if (not isfile( infile )):
        sys.exit("File " + infile + " not found. It is required to specify\n\
                the number of atoms and their charges.")
    species = []
    charges = []
    with open( infile ) as f:
        for line in f.readlines():
            if not line[0].startswith('#'):
                species.append(line.split(',')[0])
                charges.append(line.split(',')[1].rstrip())
    species = np.array(species, dtype=np.int64)  
    charges = np.array(charges, dtype=np.float64)
    return species, charges

def wccount(filename):
    import subprocess
    if filename[-3:] == '.gz':
        out = subprocess.Popen(["zcat "+ filename + " | wc -l"],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, shell = True).communicate()[0]
    else:
        out = subprocess.Popen(['wc', '-l', filename],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT
        ).communicate()[0]
    return int(out.partition(b' ')[0])


if __name__ == '__main__':
    args = parse_commandline_arguments()

    slice_size_frames = args.slicesize
    slice_offset = args.offset
    displacements_file = args.displacement_file # displacements without the header information
    msd_files = args.msd_files
    msd_length = args.msdlen
    convcalc = args.convcalc
    prntfrq = args.prntfrq
    timestep = args.timestep
    
    if (not isfile( displacements_file )):
        sys.exit("File " + displacements_file + " not found")

    try:
        nlines = wccount( displacements_file )
    except:
        sys.exit("Problem counting the number of lines in " + displacements_file )

    nspcs, charges = my_system()
    if np.sum( nspcs * charges ) != 0.0:
        sys.exit('System is not charge balanced!')

    natoms = np.sum( nspcs )
    if np.mod( nlines, natoms ) != 0:
        sys.exit('ERROR = The number of frames in the file is not an integer multiple of the number of atoms.')
    nframes_tot = int( nlines / natoms )

    displacements = pd.read_csv( displacements_file, delim_whitespace=True,names= ['x','y','z'], dtype={'x': np.float64, 'y': np.float64, 'z': np.float64} )

    if not convcalc:
        slope_stats = slope_statistics( displacements, nspcs, charges, nframes_tot, msd_length, prntfrq, timestep, slice_size_frames, slice_offset, msd_files ) 
    else:
        slope_conv = slope_convergence( displacements, nspcs, charges, nframes_tot, msd_length, prntfrq, timestep, slice_offset, msd_files )
