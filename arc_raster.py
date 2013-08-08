# arc_raster.py

#	Copyright (C) 2013 Martin D. Hurst
#	This program is free software: you can redistribute it and/or modify it 
#	under the terms of the GNU General Public License as published by the 
#	Free Software Foundation, either version 3 of the License, or (at your 
#	option) any later version. This program is distributed in the hope that 
#	it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
#	warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License along 
#	with this program. If not, see http://www.gnu.org/licenses/.
#
#	Martin D. Hurst mhurst@bgs.ac.uk

"""

Functions to link ARC rasters stored as flt or txt with numpy arrays

Martin Hurst, August 2012

"""

# Import modules
import sys, os
import numpy as np
import matplotlib.pyplot as plt

def read_ascii(ascii_file,datatype=False):

    """
    Imports raster data from ARCMap raster_to_ascii 

    INPUTS:     ascii_file is the file created by ARCMap's
                Raster to ASCii converter

    OUTPUTS:    header_data: a [1,6] vector containing
                [nrows, ncols, xllcorner, yllcorner, cellsize, NODATA_value]
                raster_data: A numpy array of shape [nrows,ncols]
                containing the raster data values

    Martin Hurst, November 2009

    Updated August 2012: Single function for all datatypes, dtype is an argument

    """
    #Check inputs

    print "Reading raster from ascii..."
    
    try: ascii_file
    except:
        sys.exit("Process Terminated: No ascii_file specified.")
    if os.path.exists(ascii_file) is False:
        sys.exit("Process Terminated: Cannot find specified file.")

    if datatype == False:
        datatype = float
        print "\tData type not set, defaulting to floating point,",
    elif datatype == float:
        print "\tData type set to float,",
    elif datatype == int:
        print "\tData type set to interger,",
    else:
        sys.exit("Process Terminated: Invalid data type. Please choose \'int\' or \'float\'")
        
    #Open the text file and read in the information
    f = open(ascii_file, 'r')
    
    # Get header data - first 6 lines of ascii file
    ncols = int((" ".join(f.readline().split())).split(" ")[1])
    nrows = int((" ".join(f.readline().split())).split(" ")[1])
    xllcorner = float((" ".join(f.readline().split())).split(" ")[1])
    yllcorner = float((" ".join(f.readline().split())).split(" ")[1])
    cellsize = float((" ".join(f.readline().split())).split(" ")[1])
    NODATA_value = int((" ".join(f.readline().split())).split(" ")[1])
    header_list = [ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value]

    f.close()
    
    # Get raster data all other data of size (nrows,ncols)
    # Set up array of zeros
    data = np.loadtxt(ascii_file,datatype,skiprows=6)
    
    #reshape to raster and set NaNs
    #data = data.reshape(nrows, ncols)
    data[data==NODATA_value] = np.nan

    print "Done."
    return header_list, data  
    
#============================================================================================

def read_float(float_file):

    """
    Imports raster data from ARCMap raster_to_float 

    INPUTS:     float_file is the file descriptor created by ARCMap's
                Raster to float converter (*.flt and *.hdr)

    OUTPUTS:    header_data: a [1,6] vector containing
                [nrows, ncols, xllcorner, yllcorner, cellsize, NODATA_value]
                raster_data: A numpy array of shape [nrows,ncols]
                containing the raster data values

    Martin Hurst, August 2012

    """
    
    #Check inputs

    print "Reading raster from float..."
    
    try: float_file
    except:
        sys.exit("Process Terminated: No float_file specified.")
    if os.path.exists(float_file) is False:
        sys.exit("Process Terminated: Cannot find specified file.")
        
    #Open the header file and read in the information
    print float_file
    print float_file[:-4]
    print float_file.rstrip(".flt")
    f = open(float_file[:-4] + ".hdr", 'r')

    ncols = int(f.readline().strip().split()[1])
    nrows = int(f.readline().strip().split()[1])
    xllcorner = float(f.readline().strip().split()[1])
    yllcorner = float(f.readline().strip().split()[1])
    cellsize = float(f.readline().strip().split()[1])
    ndv = int(f.readline().strip().split()[1])
    header_list = [ncols, nrows, xllcorner, yllcorner, cellsize, ndv]

    f.close()

    #Open the float file and read in the data
    f = open(float_file, 'rb')
    read_all = f.read()
    f.close()
    data = np.fromstring(read_all, 'f')

    #Swap if bigendian data
    if sys.byteorder == 'big':
        data = data.byteswap()

    data = data.reshape(nrows, ncols)
    data[data==ndv] = np.nan
    
    print "Done."
    
    return header_list, data  
    
#============================================================================================

def write_ascii(file_name,header,raster):

    """
    Writes a raster_ascii file for reading into ARCMap given header data and a 2d data array 

    INPUTS:     header: header data for the raster
                [ncols nrows xllcorner yllcorner cellsize NODATA_value]

    OUTPUTS:    file_name: ascii file raster

    Martin Hurst, March 2010

    """
    
    #-------------------------------------------------------------------------------

    #write header and raster to file_name    
    #create file
    f = open(file_name, 'w')
    
    #create and print header info and data
    print_header = ['ncols\t'+str(header[0])+
                    '\nnrows\t'+str(header[1])+
                    '\nxllcorner\t'+ str(header[2])+
                    '\nyllcorner\t'+ str(header[3])+
                    '\ncellsize\t'+str(header[4])+
                    '\nNODATA_value\t'+str(header[5])+'\n']
    f.write(print_header[0])

    raster[np.isnan(raster)] = header[5]
    for row in raster:
        row.tofile(f, sep="\t")
        f.write("\n")
    
    f.close

#============================================================================================

def plot_raster(raster, header=None):
    if header==None:
        raster = np.flipud(raster)
        plt.imshow(raster) #, origin='lower', interpolation='nearest')
        plt.colorbar()
        plt.show()
