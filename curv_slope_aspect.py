#=============================================================================
#
#   Function to calculate topographic curvature of a landscape from a dem in
#   ascii format as generated by ARCMap's raster2ascii converter. Also gets
#   slopes and aspect values
#
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
#
#
#   Martin Hurst, April 2010
#
#===============================================================================

# Import Modules
import os, sys, numpy as np, scipy as sp, arc_raster as asc
    
def polyfit_6term(dem_ascii, prefix, kernel_size):

    #-------------------------------------------------------------------------------------------------------------------------
    """
    Curvature calculated by fitting polynomial of the form

        A[x^2] + B[y^2] + C[x*y] + D[x] + E[y] + F
    
    to moving window kernel in the gridded data for the specified  kernel/window
    size. A-F are coefficients from which topographic parameters can be calculated.
    Curvature is determined as C = 2(AD^2 + BE^2 + CDE)/(D^2+E^2)

    INPUTS          dem_ascii: Text file containing DEM in raster format
                    prefix: Prefix for text files to be generated
                    kernel_size: Size of the moving window to fit polynomial
                                                (must be odd number)
    
    OUTPUTS         curv_out_ascii: Text file containing calculated curvatures
                    slope_out_ascii: Text file containing calculated slopes
                    aspect_out_ascii: Text file containing calculated aspects
    
        Martin Hurst, April 2010

    aspect fixed july 2011!


    """
    #-------------------------------------------------------------------------------------------------------------------------

    # Print time
    print "Running crv.polyfit_6term on " + dem_ascii + " with kernel size " + str(kernel_size) + "..."
    
    # Load DEM file
    hdr, dem = asc.read_ascii(dem_ascii)
    dem_res = hdr[4]
    ndv = hdr[5]

    print ndv
    
    # Set up DEM
    [nrows,ncols] = np.shape(dem)
    curv = np.zeros( (nrows,ncols) )
    slope = np.zeros( (nrows,ncols) )
    aspect = np.zeros( (nrows,ncols) )
    mean_residuals = np.zeros( (nrows,ncols) )

    #Generate kernel
    kernel = np.ones( (kernel_size,kernel_size) )

    #scale to midpoint of kernel
    scaling = dem_res*((kernel_size-1)/2)
    print 'scaling = ', scaling
    x = np.ones( (kernel_size,kernel_size) )
    y = np.ones( (kernel_size,kernel_size) )
    m,n = np.shape(kernel)
    for i in range(0,m):
        for j in range(0,n):
            y[i][j] = i*dem_res - scaling
            x[i][j] = j*dem_res - scaling
    x = x.reshape(kernel_size**2)
    y = y.reshape(kernel_size**2)

    centre_cell = (m-1)/2

    # Pass over DEM and fit polynomial surface to elevations in kernel each time
    # Use coefficients to calculate curvature C = 2(AD^2 + BE^2 + CDE)/(D^2+E^2)
    for i in range (0,nrows):
        #print 'row is ' + str(i)
        for j in range (0,ncols):

            #can't caluclate curvature at the edges or where NODATA_values
            if (i < centre_cell) | (j < centre_cell) | (i > nrows-1-centre_cell) | (j > ncols-1-centre_cell):
                curv[i][j] = ndv
                slope[i][j] = ndv
                aspect[i][j] = ndv
                mean_residuals[i][j] = ndv

            elif (dem[i][j] == ndv or dem[i][j] == np.nan):		
                curv[i][j] = ndv
                slope[i][j] = ndv
                aspect[i][j] = ndv
                mean_residuals[i][j] = ndv

            else:
                #clip area of raster of size of kernel
                dem_kernel = dem[i-centre_cell:i+centre_cell+1,j-centre_cell:j+centre_cell+1] 

                if (ndv in dem_kernel or np.nan in dem_kernel):
                    print 'rstr[i][j] = ',rstr[i][j]
                    curv[i][j] = ndv
                    slope[i][j] = ndv
                    aspect[i][j] = ndv
                    mean_residuals[i][j] = ndv
                else:
                    #fit polynomial to dem_kernel
                    #reshape z
                    z = dem_kernel.reshape(kernel_size**2)
                    
                    # Generate Array for fitting polynomial
                    v = np.array([np.ones(kernel_size**2), y, x, y * x, y**2, x**2])
                    # Fit least squares polynomial
                    coes, residues, rank, singval = np.linalg.lstsq(v.T, z)
                    # Calculate curvature C = 2A+2B
                    curv[i][j] = 2*coes[5] + 2*coes[4]
                    # Calculate slope S = sqrt(d^2+e^2)
                    slope[i][j] = np.sqrt(coes[2]**2 + coes[1]**2)
                    # Calculate aspect A = arctan(e/d)
                    aspect[i][j] = 180 - 57.29578 * np.arctan(coes[1]/coes[2]) + 90*(coes[2]/abs(coes[2]))
                    if aspect[i][j] < 180.0:
                        aspect[i][j] = 180.0 - aspect[i][j]
                    else:
                        aspect[i][j] = 360.0 + (180.0 - aspect[i][j])
                    # get mean residuals
                    #mean_residuals[i][j] = np.mean(residues)
                    
    #------------------------------------
    #print curvature data to file
    #------------------------------------

    curv_out_ascii = prefix + "curv.txt"
    slope_out_ascii = prefix + "slope.txt"
    aspect_out_ascii = prefix + "aspect.txt"
    
    #write ascii raster function in ascii_raster.py
    asc.write_ascii(curv_out_ascii,hdr,curv)
    asc.write_ascii(slope_out_ascii,hdr,slope)
    asc.write_ascii(aspect_out_ascii,hdr,aspect)
    #asc.write_ascii_raster(hdr,mean_residuals,residuals_ascii)
    
if __name__ == "__main__":
    
    folder = "D:\\DBPR\\testing\\curv_slope_aspect\\"   #Working directory
    prefix = "db_"          #project name e.g. "db_"
    dem = "db_dem.txt"      #name of dem (ascii grid)
    scale = 13				#Size of square window (pixels) centred on pixel of interest over which to calc topographic properties
    polyfit_6term(dem,prefix,scale) #run