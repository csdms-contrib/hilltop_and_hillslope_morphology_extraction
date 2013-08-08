/*

read_dem_v1.3.h functions to read and write ascii and flt gridded data files in a
format compatible with ARC GIS.

Copyright (C) 2013 Martin D. Hurst

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Martin D. Hurst
mhurst@bgs.ac.uk

*/

// HEADER FILE for read_dem_v1.3
//==============================================================================
//
// Reads and writes raster files in ascii or float format as created by ARCMap
//
//	read_raster_txt : 	Reads a raster in ascii format. Data is stored as float
//						so that it can later be written to a *.flt file
//
//						Requires arguments for ncols, nrows (both ints), xmin, ymin
//						(both doubles), dem_res (float), no data value ndv (int) and
//						the data array (TNT:Array2D) into which to read the data, and
//						the char string containing the file name to read.
//
//	write_raster_txt :	Writes a raster to an ARCMap style ascii format.
//
//						Requires all header data (see read_raster_txt INPUTS) plus the
//						TNT:Array2D containing the raster data (type float) and a char
//						string containing the filename to write to.
//
//	read_raster_hdr :	Reads the header file (*.hdr) for a float
//
//						Requires arguments for ncols, nrows (both ints), xmin,
//						ymin (both doubles), dem_res (float) and no data value ndv (int)
//						into which to read the data, and the char string containing the
//						file name to read.
//
//	read_raster_flt : 	Reads the data file (*.flt) for a float
//
//						Requires nrows and ncols as obtained using read_raster_hdr,
//						a TNT:Array2D into which to read the data, and a char string
//						containing the filename to be read (*.flt)
//
//	write_raster_hdr :	Writes the header file (*.hdr) for a float
//
//						Requires ncols, nrows (both ints), xmin, ymin (both doubles),
//						dem_res (float) and no data value ndv (int), and the char
//						string containing the file name to write.
//
//	write_raster_flt :	Writes the data file (*.flt) for a float
//
//						Requires nrows and ncols as obtained using read_raster_hdr,
//						a TNT:Array2D containing the data, and a char string
//						containing the filename to be written (*.flt)
//
//
//------------------------------------------------------------------------------
// David Milodowski, February 2012
//------------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <iomanip>
#include "TNT/tnt.h"

#ifndef READ_DEM_H
#define READ_DEM_H

// ASCII tools
  void read_raster_txt(int& ncols, int& nrows, double& xmin, double& ymin, float& dem_res, int& ndv, TNT::Array2D<float>& data, char *in_file);
  void write_raster_txt(int& ncols, int& nrows, double& xmin, double& ymin, float& dem_res, int& ndv, TNT::Array2D<float>& data, char *out_file);
// Float tools
  void read_raster_hdr(int& ncols, int& nrows, double& xmin, double& ymin, float& dem_res, int& ndv, char *header_file);
  void read_raster_flt(int& ncols, int& nrows, TNT::Array2D<float>& data, char *in_file);
  void write_raster_hdr(int& ncols, int& nrows, double& xmin, double& ymin, float& dem_res, int& ndv, char *header_file);
  void write_raster_flt(int& ncols, int& nrows, TNT::Array2D<float>& data, char *out_file);

#endif

