//INCLUDE MODULES
#include <iostream>
#include <fstream>
#include <iomanip>
#include "TNT/tnt.h"

//SET NAME SPACE
using namespace std;
using namespace TNT;



//========================================================================================
//
// Reads and writes raster files in ascii or float format as created by ARCMap
//
//Copyright (C) 2013 Martin D. Hurst
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//Martin D. Hurst
//mhurst@bgs.ac.uk
//
//	read_raster_txt : 	Reads a raster in ascii format. Data is stored as float in order
//						that it can later be written to a *.flt file
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
//----------------------------------------------------------------------------------------
//
// V1.2 November 2010	Updated to allow reading and writing of *.flt format with header.
//						Functions added to read and write integer floats.
//
//----------------------------------------------------------------------------------------
//
// Martin Hurst, October 2010
//
//========================================================================================



void read_raster_txt(int& ncols, int& nrows, double& xmin, double& ymin,
						float& dem_res, int& ndv, Array2D<float>& data, char *in_file) {

	ifstream data_in(in_file);
	//Read in raster data
	string str;
	data_in >> str >> ncols >> str >> nrows
		   >> str >> xmin >> str >> ymin
		   >> str >> dem_res
		   >> str >> ndv;

	for (int i=0; i<nrows; ++i)	{
		for (int j=0; j<ncols; ++j)	{
			data_in >> data[i][j];
		}
	}
	data_in.close();
}


void write_raster_txt(int& ncols, int& nrows, double& xmin, double& ymin,
						float& dem_res, int& ndv, Array2D<float>& data, char *out_file) {

	ofstream data_out(out_file);

	if( data_out.fail() ){
		cout << "\nFATAL ERROR: unable to write to " << out_file << endl;
		exit(EXIT_FAILURE);
	}

	data_out <<  "ncols         " << ncols
			<< "\nnrows         " << nrows
			<< "\nxllcorner     " << setprecision(14) << xmin
			<< "\nyllcorner     " << setprecision(14) << ymin
			<< "\ncellsize      " << dem_res
			<< "\nNODATA_value  " << ndv << endl;

	for (int i=0; i<nrows; ++i)	{
		for (int j=0; j<ncols; ++j)	{
			data_out << setprecision(6) << data[i][j] << " ";
		}
		if (i != nrows-1) data_out << endl;
	}
	data_out.close();
}

void read_raster_hdr(int& ncols, int& nrows, double& xmin, double& ymin,
						float& dem_res, int& ndv, char *header_file) {

	ifstream ifs(header_file);
	if( ifs.fail() ){
		cout << "\nFATAL ERROR: the header file \"" << header_file << "\" doesn't exist" << endl;
		exit(EXIT_FAILURE);
	}
	else {
	string str;
	ifs >> str >> ncols >> str >> nrows
		>> str >> xmin >> str >> ymin
		>> str >> dem_res
		>> str >> ndv;
	}
	ifs.close();
}

void read_raster_flt(int& ncols, int& nrows, Array2D<float>& data, char *in_file) {

	ifstream ifs(in_file, ios::in | ios::binary);
	if( ifs.fail() ){
		cout << "\nFATAL ERROR: the data file \"" << in_file << "\" doesn't exist" << endl;
		exit(EXIT_FAILURE);
	}
	else {
		float temp;
		for (int i=0; i<nrows; ++i)	{
			for (int j=0; j<ncols; ++j)	{
				ifs.read(reinterpret_cast<char*>(&temp), sizeof(temp));
				data[i][j] = temp;
			}
		}
	}
	ifs.close();
}

void write_raster_hdr(int& ncols, int& nrows, double& xmin, double& ymin,
						float& dem_res, int& ndv, char *header_file) {
	ofstream ofs(header_file);
	string str;
	ofs <<  "ncols         " << ncols
		<< "\nnrows         " << nrows
		<< "\nxllcorner     " << setprecision(14) << xmin
		<< "\nyllcorner     " << setprecision(14) << ymin
		<< "\ncellsize      " << dem_res
		<< "\nNODATA_value  " << ndv << endl;
	ofs.close();
}

void write_raster_flt(int& ncols, int& nrows, Array2D<float>& data, char *out_file) {
	ofstream ofs(out_file, ios::out | ios::binary);
	float temp;
	for (int i=0; i<nrows; ++i)	{
		for (int j=0; j<ncols; ++j)	{
			temp = data[i][j];
			ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
		}
	}
	ofs.close();
}
