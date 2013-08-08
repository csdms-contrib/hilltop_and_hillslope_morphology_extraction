#===============================================================================
#
# Python functions to extract hilltops from a DEM as the intesecting margins of
# adjoining drainage basins of a given stream order.
#
# Martin Hurst, 2012
#
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
#
#===============================================================================
import sys, os, shapefile, numpy as np

# Import arcpy module
import arcpy

# Check out any necessary licenses
arcpy.CheckOutExtension("spatial")
arcpy.env.overwriteOutput=True

def hydro_prep(folder, prefix, dem, area_thresh):

	import string

	"""
	HYDRO_PREP

	Runs a series of functions necessary for the extraction fluvial data
	including fill, flow direction, flow accumulation, creating stream network
	and stream order

	Inputs: folder: directory location of DEM file
		    prefix: a file prefix (2 letters)
		    dem_file: the name of your dem file
		    area_thresh: area threshold for channel initiation in pixels

	Outputs: A series of hydrology rasters stored in folder//hydro//
		     returns max_SO to allow continuation to get basins

	Martin Hurst, February 2010
		max_SO return added May 2010
	"""

	#print "\nhydro_prep: Running hydrology preperation tool on ", dem

	# Set the Geoprocessing environment...
	#dem_file = folder + dem_file
	#print dem_file
	#gp.extent = dem

	# Allow Overwrite
	#gp.OverWriteOutput = 1

	# Set up workspace

	hydro_folder = folder + "\\hydro\\"
	junk_folder = hydro_folder + "\\junk\\"

	if os.path.exists(hydro_folder):
		pass
	else:
		os.mkdir(hydro_folder)
		os.mkdir(junk_folder)

	# Setup variables

	#dem = folder + "\\" + prefix + "dem"
	fill = hydro_folder + prefix + "fill"
	flow_dir = hydro_folder + prefix + "dir"
	flow_acc = hydro_folder + prefix + "acc"
	strm_net_ras = hydro_folder + prefix + "stn"
	strm_net_shp = hydro_folder + prefix + "stn.shp"
	strm_net_dbf = hydro_folder + prefix + "stn.dbf"
	watershed_ras = hydro_folder + prefix + "ws"
	strm_ord = hydro_folder + prefix + "so"
	strm_ord_shp = hydro_folder + prefix + "so.shp"
	max_acc_ras = junk_folder + prefix + "mxacc"
	max_acc_pnt = junk_folder + prefix + "mxaccp"
	zone4stats = junk_folder + prefix + "zone"
	dem_ascii = folder + prefix + "dem.txt"
		
	print '\tProcess: Fill...'
	arcpy.gp.Fill_sa(dem, fill, "")

	print '\tProcess: Flow Direction...'
	arcpy.gp.FlowDirection_sa(fill, flow_dir, "NORMAL", "")

	print '\tProcess: Flow Accumulation...'
	arcpy.gp.FlowAccumulation_sa(flow_dir, flow_acc, "", "FLOAT")

	# Process: Con...
	condition = "VALUE > " + str(area_thresh)
	arcpy.gp.Con_sa(flow_acc, 1, strm_net_ras, "", condition)

	print '\tProcess: Stream Order...'
	arcpy.gp.StreamOrder_sa(strm_net_ras, flow_dir, strm_ord, "STRAHLER")

	# Process: Raster to Polyline...
	arcpy.gp.RasterToPolyline_conversion(strm_ord, strm_ord_shp, "ZERO", "0", "NO_SIMPLIFY", "VALUE")

	# Process: Stream to Feature...
	arcpy.gp.StreamToFeature_sa(strm_ord, flow_dir, strm_net_shp, "NO_SIMPLIFY")


	print "\tDone!"

#==============================================================================
	
def basins_SO(folder,prefix,SO_extract,dem_res):

    """
    basins_SO

    Extracts basins form a DEM of a given Strahler order. Basin pour points defined
    by lowest end of strahler order stream segment

    INPUTS:     
                folder:     Location of hydro prop rasters (flow_dir,flow_acc,strm_ord required)
                prefix:     project identifier
                SO_extract:     Strahler order for which to extract basins
                dem_res:   DEM resolution

    OUTPUTS:
                creates folder\\basins\\ in which basins shapefiles are created
                creates folder\\basins\\junk\\ in which intermediate files are located
   
   
   Martin Hurst, Dec 2009
   
   """

    print "\n\tbasins_SO: Extracting basins at stream order = ", SO_extract
    
    
    # Set up workspace
    
    basins_folder = folder + "\\basins\\"
    basins_junk = basins_folder + "junk\\"
    hydro_folder = folder + "\\hydro\\"
    
    if os.path.exists(basins_folder):
        pass
    else:
        os.mkdir(basins_folder)
        os.mkdir(basins_junk)
        
      
    # Local variables...
    # Input rasters
    flow_dir = hydro_folder + prefix + "dir"
    flow_acc = hydro_folder + prefix + "acc"
    strm_ord = hydro_folder + prefix + "so"
        
    #Other parameters
    value_SO_extract = "VALUE = " + str(SO_extract)
        
    #Variables to create
    selected_SO = basins_junk + prefix + "SO" + str(SO_extract)
    selected_SO_shp = basins_junk + prefix + "SO" + str(SO_extract) + ".shp"
    SO_endpoints = basins_junk + prefix + "SO" + str(SO_extract) + "end.shp"
    SO_endpoints_buff = basins_junk + prefix + "SO" + str(SO_extract) + "buff.shp"
    buffered_SO = basins_junk + prefix + "SO" + str(SO_extract) + "buft.shp"
    bsn_pours = basins_junk + prefix + "bsn" + str(SO_extract) + "pours.shp"
    bsn_snaps = basins_junk + prefix + "bsn" + str(SO_extract) + "snp"
    basins = basins_junk + prefix + "SO" + str(SO_extract) + "bsns"
    basins_shp = basins_folder + prefix + "SO" + str(SO_extract) + "bsns.shp"
    
    #print "\t\t Processing: Creating buffer..."   
    #print "Process: Con..."
    arcpy.gp.Con_sa(strm_ord, 1, selected_SO, "", value_SO_extract)

    #print "Process: Stream to Feature..."
    arcpy.gp.StreamToFeature_sa(selected_SO, flow_dir, selected_SO_shp, "NO_SIMPLIFY")

    #print "Process: Feature Vertices To Points..."
    arcpy.gp.FeatureVerticesToPoints_management(selected_SO_shp, SO_endpoints, "END")

    #print "Process: Buffer..."
    arcpy.gp.Buffer_analysis(SO_endpoints, SO_endpoints_buff, str(float(dem_res)*1.5)+" Meters", "FULL", "ROUND", "NONE", "")

    #print "Process: Erase..."
    arcpy.gp.Erase_analysis(selected_SO_shp, SO_endpoints_buff, buffered_SO, "")

    #print "\t\t Processing: Creating pourpoints..."
    #print "Process: Feature Vertices To Points..."
    arcpy.gp.FeatureVerticesToPoints_management(buffered_SO, bsn_pours, "END")

    #print  "Process: Snap Pour Point..."
    arcpy.gp.SnapPourPoint_sa(bsn_pours, flow_acc, bsn_snaps, float(dem_res)/2.0, "ARCID")

    print "\t\t Processing: Extracting basins..."

    #print "Process: Watershed..."
    arcpy.gp.Watershed_sa(flow_dir, bsn_snaps, basins, "Value")

    # Process: Raster to Polygon...
    arcpy.gp.RasterToPolygon_conversion(basins, basins_shp, "NO_SIMPLIFY")
    print "\t\tDone"

#================================================================================
    
def basins_manySO(folder,prefix,dem_res):

    """
    basins_manySOs

    Iterates through stream order range and extracts basins for each
    by running basins_SO

    INPUTS:     
                folder:     Location of hydro prep rasters (flow_dir,flow_acc,strm_ord required)
                prefix:     project identifier
                dem_res:    DEM resolution

    OUTPUTS:
                creates folder\\basins\\ in which basins shapefiles are created
                creates folder\\basins\\junk\\ in which intermediate files are located
    """

    hydro_folder = folder + "\\hydro\\"
    strm_net_shp = hydro_folder + prefix + "stn.shp"
    sf = shapefile.Reader(strm_net_shp)
    records = sf.records()
    n = len(records)

    so = np.zeros(n)
    for i in range (0,n):
        so[i] = records[i][1]

    max_SO = int(np.max(so))

    print "\nbasins_manySO: Extracting basins with stream orders 1 to ", max_SO

    #iterate through stream orders getting basins
    i = 1
    while (i <= max_SO):
        basins_SO(folder,prefix,i,dem_res)
        i = i+1

#================================================================================
        
def get_hilltops(folder,prefix,buff_dist_ht,slope_thresh,slope_txt_in='none',basins_shape=None):
    """
    Define locations of hilltops from basin edges with buffer
    buffer channel by 10m to avoid interference when
    calculating curvatures and slopes

    basin_SOs, integer for maximum stream order to consider
    
    Martin Hurst, March 2010

    Updated July 2010 to buffer away from stream network

    Updated July 2011 to use intersect analysis to speed up operation
    """

    print "Running get_hilltops..."
    
    # Local Variables

    #folders
    hydro_folder = folder + "hydro\\"
    hydro_junk = hydro_folder + "junk\\"
    basins_folder = folder + "basins\\"
    basins_junk = basins_folder + "junk\\"
    hilltop_folder = folder + "hilltops\\"
    hilltop_junk = hilltop_folder + "junk\\"
    
    # check folders exist and create if not
    if os.path.exists(hilltop_folder):
        pass
    else:
        os.mkdir(hilltop_folder)

    # check folders exist and create if not
    if os.path.exists(hilltop_junk):
        pass
    else:
        os.mkdir(hilltop_junk)

    #files
    expression_string = list()

    dem_file = folder + prefix + "dem"
    strm_net_shp = hydro_folder + prefix + "stn.shp"
    raw_hilltops =  hilltop_junk + prefix + "raw_hts.shp"
    hilltops_buff1_pl = hilltop_junk + prefix + "hts_buff1.shp"
    hilltops_buff2_pl = hilltop_junk + prefix + "hts_buff2.shp"
    hts_ras = hilltop_junk + prefix + "rast_hts"
    hilltops_pl = hilltop_folder + prefix + "hilltops.shp"
    hilltop_buff_pg = hilltop_folder + prefix + "hilltops_buff.shp"
    
    slope_txt_in = folder + slope_txt_in
    slope_percent = hilltop_junk + prefix + "slope_per"
    slope_raster = hilltop_folder + prefix + "slope"
    masked_slope_raster = hilltop_junk + prefix + "slp_mask"
    masked_slope_poly = hilltop_junk + prefix + "slp_mask.shp"
    stream_net_shape = hydro_folder + prefix + "stn.shp"
    stream_net_buff = hydro_junk + prefix + "stn_buff.shp"

    if basins_shape == None:
        #get max_SO
        sf = shapefile.Reader(strm_net_shp)
        records = sf.records()
        n = len(records)

        so = np.zeros(n)
        for i in range (0,n):
            so[i] = records[i][1]

        max_SO = int(np.max(so))

        i = 1

        while i <= max_SO:
            basins = basins_folder + prefix + "so" + str(i) + "bsns.shp"
            print "\t processing basins of stream order " +str(i)
            intersect = hilltop_junk + "int_p" + str(i) + ".shp"
                
            print "\t\tlooking for intersects"

            # INTERSECT ANALYSIS - RETURNS LINES WHERE POLYS MEET
            arcpy.gp.Intersect_analysis(basins, intersect, "", "0", "LINE")

            expression_string.insert(i, intersect)
            i = i + 1

    else:
        i = 1
        intersect = hilltop_junk + "int_p" + str(i) + ".shp"
                
        print "\t\tlooking for intersects"

        # INTERSECT ANALYSIS - RETURNS LINES WHERE POLYS MEET
        arcpy.gp.Intersect_analysis(basins_shape, intersect, "", "0", "LINE")
        
        expression_string.insert(i, intersect)
       
    # MEREGES THE LINES TOGETHER
    files = list()
              
    w = shapefile.Writer()
    
    for f in expression_string:
        try:
            r = shapefile.Reader(f)
            w._shapes.extend(r.shapes())
            w.records.extend(r.records())
            w.fields = list(r.fields)
        except:
            pass

    w.save(raw_hilltops)

    #--- Filter hilltops not to include slopes > x ---#
    
    # Process: ASCII to Raster...
    if os.path.exists(slope_raster):
        print "\tslope_raster exists"
    elif os.path.exists(slope_txt_in):
        arcpy.gp.ASCIIToRaster_conversion(slope_txt_in, slope_raster, "FLOAT")
        print '\tconverted slope raster'
    else:
        print "\tCacluating Slope over 3x3m window"
        arcpy.gp.Slope_sa(dem_file, slope_raster, "PERCENT_RISE",0.01)
        

    # Process: Con...
    con = "VALUE > " + str(slope_thresh)
    arcpy.gp.Con_sa(slope_raster, "1", masked_slope_raster, "", con)

    # Process: RasterToPolygon_conversion
    arcpy.gp.RasterToPolygon_conversion(masked_slope_raster, masked_slope_poly, "NO_SIMPLIFY")

    # Process: Buffer stream network to 30m...
    arcpy.gp.Buffer_analysis(stream_net_shape, stream_net_buff, 5, "FULL", "ROUND", "NONE", "")

    # Process: Mask hilltops <Xm from a stream...
    arcpy.gp.Erase_analysis(raw_hilltops, stream_net_buff, hilltops_buff1_pl, "")

    # Process: Mask hilltops for slopes >0.4 ...
    arcpy.gp.Erase_analysis(hilltops_buff1_pl, masked_slope_poly, hilltops_buff2_pl, "")

    
    #CONVERTS TO RASTER
    print "\tdeleting duplicates"
    arcpy.gp.FeatureToRaster_conversion(hilltops_buff2_pl, "GRIDCODE", hts_ras, 0.1)
    arcpy.gp.RasterToPolyline_conversion(hts_ras, hilltops_pl, "NODATA", 30)
    

    # Process: Buffer basin margins to define hilltops...
    arcpy.gp.Buffer_analysis(hilltops_pl, hilltop_buff_pg, 0.5, "FULL", "ROUND", "NONE", "")


    print "\tDone!\n"

#===============================================================================
    
def get_cht(folder,prefix,curv_txt_in='none'):
    
    """
    Creates raster of hilltop curvature from ascii curv_file created by curv_func_v1.1.cpp
    requires get_hilltops to have been run first
    
    Martin Hurst, March 2010
    
    """

    print "Running get_cht..."
    

    # Files
    hilltop_folder = folder + "hilltops\\"
    final_hilltop_buff = hilltop_folder + prefix + "hilltops_buff.shp"
    hilltop_curv_raster = hilltop_folder + prefix + "cht"
    curv_txt_in = folder + curv_txt_in
    
    # check folders exist and create if not
    if os.path.exists(hilltop_folder):
        pass
    else:
        os.mkdir(hilltop_folder)

    curv_raster = hilltop_folder + prefix + "curv"

    if os.path.exists(curv_txt_in):
        arcpy.gp.ASCIIToRaster_conversion(curv_txt_in, curv_raster, "FLOAT")
        print '\tconverted curv raster'
    elif os.path.exists(curv_raster):
        print "\tcurv_raster exists"
    else:
        raise runtimeerror('No curvature data')
        
    # Process: Extract by Mask (2)...
    arcpy.gp.ExtractByMask_sa(curv_raster, final_hilltop_buff, hilltop_curv_raster)

    print "\tDone!\n"

#=============================================================================
    
# This allows hydro_prep to be launched from a command line
if __name__ == "__main__":
    folder = "E:\\DBPR\\testing\\"   #Working directory
    prefix = "db_"          #project name e.g. "db_"
    dem = "db_dem"          #name of dem (arc grid)
    dem_res = 0.25          #dem resolution (e.g. 1m)
    area_thresh = 10000     #area threshold (pixels) for channels

    buff_distance = 0.5             #distance to buffer stream net by
    curv_txt_in = "db_curv.txt"     #curvature ascii grid
    slope_txt_in = "db_slope.txt"   #slope ascii grid
    slope_thresh = 0.4              #mask steep hilltops
	
    hydro_prep(folder, prefix, dem, area_thresh)
    basins_manySO(folder,prefix,dem_res)
    get_hilltops(folder,prefix,buff_distance,slope_thresh,slope_txt_in)
    get_cht(folder,prefix,curv_txt_in)

