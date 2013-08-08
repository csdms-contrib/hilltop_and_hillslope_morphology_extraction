/*--------------------------------------------------------------------------------

Kinematic Trace Downstream using Aspect

takes path of steepest descent routing according to aspect from hilltop
to valley bottom and gets mean slope, hillslope length, and relief.
Algorithm created after Lea [1991] 'An aspect driven kinematic routing algorithm'

Requires a dtm, stream network and hilltop network to run, as well as curvature, slope
and aspect grids. Grids are read as *.flts

Martin Hurst. July 2011;

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

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <cstring>
#include "TNT/tnt.h"
#include "read_dem_v1.3.hpp"

using namespace std;
using namespace TNT;

int main(int nNumberofArgs,char *argv[]) {

	//Test for correct input arguments
	if (nNumberofArgs!=2) {
		cout << "FATAL ERROR: You must provide an input arguments\n\t - a dem descriptor\n\t" << endl;
		exit(EXIT_SUCCESS);
	}

	//get the file names as arguments
	char *dem_in = argv[1];

	//generate file names
	int slen = strlen(dem_in);

	char *fill_file = new char[slen+9], *fill_header = new char[slen+9];
	char *stnet_header = new char[slen+11], *stnet_file = new char[slen+11];
	char *slope_header = new char[slen+10], *slope_file = new char[slen+10];
	char *aspect_header = new char[slen+12], *aspect_file = new char[slen+12];
	char *seg_header = new char[slen+8], *seg_file = new char[slen+8];
	char *cht_header = new char[slen+9], *cht_file = new char[slen+9];

	strcpy(fill_file,dem_in), strcat(fill_file,"_fill.flt");
	strcpy(fill_header,dem_in), strcat(fill_header,"_fill.hdr");
	strcpy(stnet_header,dem_in), strcat(stnet_header,"_stnet.hdr");
	strcpy(stnet_file,dem_in), strcat(stnet_file,"_stnet.flt");
	strcpy(slope_header, dem_in), strcat(slope_header,"_slope.hdr");
	strcpy(slope_file, dem_in), strcat(slope_file,"_slope.flt");
	strcpy(aspect_header,dem_in), strcat(aspect_header,"_aspect.hdr");
	strcpy(aspect_file,dem_in), strcat(aspect_file,"_aspect.flt");
	strcpy(seg_header,dem_in), strcat(seg_header,"_seg.hdr");
	strcpy(seg_file,dem_in), strcat(seg_file,"_seg.flt");
	strcpy(cht_header,dem_in), strcat(cht_header,"_curv.hdr");
	strcpy(cht_file,dem_in), strcat(cht_file,"_curv.flt");

	//Declare parameters
	//Declare parameters
	int ncols,nrows,i,j,a,b;
	double xmin,ymin,X,Y;
	float dem_res, mean_slope, relief;
	int ndv;
	double slope_total, length, d;
	int flag, count;
	double PI = 3.14159265;
	double degs, degs_old, degs_new, theta;
	double s_local, s_edge;
	double xo, yo, xi, yi, temp_yo1, temp_yo2, temp_xo1, temp_xo2;

	// a direction flag numbered 1,2,3,4 for E,S,W,N respectively
	int dir;

	//READ IN DATA
	//read in header
	read_raster_hdr(ncols, nrows, xmin, ymin, dem_res, ndv, fill_header);

	double ymax = ymin + nrows*dem_res;
	double xmax = xmin + ncols*dem_res;

	//Declare Arrays

	Array2D<float> zeta(nrows, ncols);
	Array2D<float> stnet(nrows,ncols);
	//Array2D<float> area(nrows,ncols);
	Array2D<float> slope(nrows,ncols);
	Array2D<float> aspect(nrows,ncols);
	Array2D<float> rads(nrows,ncols);
	Array2D<float> cht(nrows,ncols);
	Array2D<float> seg(nrows,ncols);
	Array2D<float> path(nrows, ncols);
	Array2D<float> blank(nrows,ncols);

	int vec_size = 1000000;

	Array1D<double> easting(ncols);
	Array1D<double> northing(nrows);
	Array1D<double> east_vec(vec_size);
	Array1D<double> north_vec(vec_size);

	//read in arrays
	read_raster_flt(ncols, nrows, zeta, fill_file);
	//read_raster_flt(ncols, nrows, area, area_file);
	read_raster_flt(ncols, nrows, stnet, stnet_file);
	read_raster_flt(ncols, nrows, slope, slope_file);
	read_raster_flt(ncols, nrows, aspect, aspect_file);
	read_raster_flt(ncols, nrows, cht, cht_file);
	read_raster_flt(ncols, nrows, seg, seg_file);

	cout	<< 	"\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
			<<	"\n* Kinematic trace: steepest descent routing according to aspect   *"
			<<	"\n*                  from hilltop to valley bottom and collects     *"
			<<	"\n*                  gets mean slope, hillslope length, and relief. *"
			<<	"\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << endl;


	//Setup output files
	char *outfile = new char[26];
	strcpy(outfile,dem_in), strcat(outfile,"_hillslope_metrics.txt");

	ofstream ofs;
	ofs.open(outfile);
	if( ofs.fail() ){
		cout << "\nFATAL ERROR: unable to write to " << outfile << endl;
		exit(EXIT_FAILURE);
	}
	ofs << "X Y Seg Cht S R Lh A\n";

	//for creating file names for hillslope profiles
	string file_part_1 = "prof_";
	string file_part_2;
	string file_part_3 = ".txt";
	string filename;
	char filename_c[256];

	//calculate northing and easting
	for (i=0;i<nrows;++i){
			northing[i] = ymax - i*dem_res - 0.5;
	}
	for (j=0;j<ncols;++j){
			easting[j] = xmin + j*dem_res + 0.5;
	}

	ofstream ofs2;
	ofs2.open("seg_recorder.txt");
	ofs2 << "segs" << endl;

	//convert aspects to radians with east as theta = 0/2*pi
	for (i=0; i<nrows; ++i) {
		for (j=0; j<ncols; ++j) {
			//convert aspects to radians with east as theta = 0/2*pi
			if (rads[i][j] != ndv) rads[i][j] = (PI/180)*((-1*aspect[i][j]) + 90);
			blank[i][j] = ndv;
			if (seg[i][j] != ndv) ofs2 << seg[i][j] << endl;
		}
	}

	ofs2.close();

	int ht_count = 0;

	// cycle through study area, find hilltops and trace downstream
	for (i=1; i<nrows-1; ++i) {
		for (j=1; j<ncols-1; ++j) {
			cout << flush << "\t\ti = " << i << " / " << nrows << "\r";
			// ignore edge cells and non-hilltop cells
			// route initial node by aspect and get outlet coordinates
			if (cht[i][j] != ndv) {
				//cout << "\n\t\t" << seg[i][j] << " " << cht[i][j] << " " << i << " " << j << endl;
				slope_total = 0;
				length = 0;
				flag = true;
				count = 0;
				path = blank.copy();


				++ht_count;
				count = 1;
				degs = aspect[i][j];
				theta = rads[i][j];
				a = i;
				b = j;
				path[a][b] = 1;

				east_vec[0] = easting[b];
				north_vec[0] = northing[a];

				s_local = slope[a][b];
				//test direction, calculate outlet coordinates and update indicies
				// easterly
				if (degs >= 45 && degs < 135) {
					xo = 1, yo = (1+tan(theta))/2;
					d = abs(1/(2*cos(theta)));
					xi = 0, yi = yo;
					dir = 1;
					east_vec[count] = easting[b] + 0.5*dem_res;
					north_vec[count] = northing[a] + yo - 0.5*dem_res;
					++b;
				}
				//southerly
				else if (degs >= 135 && degs < 225) {
					xo = (1-(1/tan(theta)))/2, yo = 0;
					d = abs(1/(2*cos((PI/2)-theta)));
					xi = xo, yi = 1;
					dir = 2;
					east_vec[count] = easting[b] + xo - 0.5*dem_res;
					north_vec[count] = northing[a] - 0.5*dem_res;
					++a;
				}
				// westerly
				else if (degs >= 225 && degs < 315) {
					xo = 0, yo = (1-tan(theta))/2;
					d = abs(1/(2*cos(theta)));
					xi = 1,	yi = yo;
					dir = 3;
					east_vec[count] = easting[b] -0.5*dem_res;
					north_vec[count] = northing[a] + yo - 0.5*dem_res;
					--b;
				}
				//northerly
				else if (degs >= 315 || degs < 45) {
					xo = (1+(1/tan(theta)))/2, yo = 1;
					d = abs(1/(2*cos((PI/2) - theta)));
					xi = xo, yi = 0;
					dir = 4;
					east_vec[count] = easting[b] + xo - 0.5*dem_res;
					north_vec[count] = northing[a] + 0.5*dem_res;
					--a;
				}
				else {
					cout << "FATAL ERROR, Kinematic routing algorithm encountered null aspect value" << endl;
					exit(EXIT_FAILURE);
				}

				//collect slopes and totals weighted by path length
				slope_total += s_local*d;
				length += d;

				s_local = slope[a][b];


				//continue trace until a stream node is encountered
				while (flag == true) {

					path[a][b] = 1;

					degs_old = degs;
					degs_new = aspect[a][b];
					theta = rads[a][b];

					++ count;

				//Test for perimeter flow paths
					if ((dir == 1 && degs_new > 0 && degs_new < 180)
						||	(dir == 2 && degs_new > 90 && degs_new < 270)
						||	(dir == 3 && degs_new > 180 && degs_new < 360)
						||	(dir == 4 && degs_new > 270) 
						|| (dir == 4 && degs_new < 90)) {

						//DO NORMAL FLOW PATH
						//set xo, yo to 0 and 1 in turn and test for true outlet (xi || yi == 0 || 1)
						temp_yo1 = yi + (1-xi)*tan(theta); 		// xo = 1
						temp_xo1 = xi + (1-yi)*(1/tan(theta)); 	// yo = 1
						temp_yo2 = yi - xi*tan(theta);			// xo = 0
						temp_xo2 = xi - yi*(1/tan(theta));		// yo = 0

						// can't outlet at same point as inlet
						if (dir == 1) temp_yo2 = -1;
						else if (dir == 2) temp_xo1 = -1;
						else if (dir == 3) temp_yo1 = -1;
						else if (dir == 4) temp_xo2 = -1;

						s_local = slope[a][b];

						if (temp_yo1 <= 1 && temp_yo1 > 0) {
							xo = 1, yo = temp_yo1;
							d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
							xi = 0, yi = yo,
							dir = 1;
							//east_vec[count] = easting[b] + 0.5*dem_res;
							//north_vec[count] = northing[a] + yo - 0.5*dem_res;
							++b;
							if (xi== 0 && yi == 0) yi = 0.00001;
							else if (xi== 0 && yi == 1) yi = 1 - 0.00001;
						}
						else if (temp_xo2 <= 1 && temp_xo2 > 0) {
							xo = temp_xo2, yo = 0;
							d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
							xi = xo, yi = 1,
							dir = 2;
							//east_vec[count] = easting[b] + xo - 0.5*dem_res;
							//north_vec[count] = northing[a] - 0.5*dem_res;
							++a;
							if (xi== 0 && yi == 1) xi = 0.00001;
							else if (xi== 1 && yi == 1) xi = 1 - 0.00001;
						}
						else if (temp_yo2 <= 1 && temp_yo2 > 0) {
							xo = 0, yo = temp_yo2;
							d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
							xi = 1, yi = yo,
							dir = 3;
							//east_vec[count] = easting[b] -0.5*dem_res;
							//north_vec[count] = northing[a] + yo - 0.5*dem_res;
							--b;
							if (xi== 1 && yi == 0) yi = 0.00001;
							else if (xi== 1 && yi == 1) yi = 1 - 0.00001;
						}

						else if (temp_xo1 <= 1 && temp_xo1 > 0) {
							xo = temp_xo1, yo = 1;
							d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
							xi = xo, yi = 0,
							dir = 4;
							//east_vec[count] = easting[b] + xo - 0.5*dem_res;
							//north_vec[count] = northing[a] + 0.5*dem_res;
							--a;
							if (xi == 0 && yi == 0) xi = 0.00001;
							else if (xi== 1 && yi == 0) xi = 1 - 0.00001;
						}
						slope_total += s_local*d;
					}

					else {

						// ROUTE ALONG EDGES
						if (dir	== 1) {
							if 	(degs_old <= 90 || degs_new >= 270) {
								xo = 0.00001, yo = 1;
								s_edge = abs(s_local*sin(theta));
								d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
								xi = xo, yi = 1-yo;
								dir = 4;
								//east_vec[count] = easting[b] + xo - 0.5*dem_res;
								//north_vec[count] = northing[a] + 0.5*dem_res;
								--a;
							}
							else if (degs_old > 90 || degs_new < 270) {
								xo = 0.00001, yo = 0;
								s_edge = abs(s_local*sin((PI/2)-theta));
								d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
								xi = xo, yi = 1-yo;
								dir = 2;
								//east_vec[count] = easting[b] + xo - 0.5*dem_res;
								//north_vec[count] = northing[a] - 0.5*dem_res;
								++a;
							}
							else {
								cout << "Flow unable to route N or S" << endl;
								exit(EXIT_FAILURE);
							}
						}
						else if (dir == 2) {
							if 	(degs_old <= 180 || degs_new >= 0) {
								xo = 1, yo = 1-0.00001;
								s_edge = abs(s_local*sin((2/PI)-theta));
								d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
								xi = 1-xo, yi = yo;
								dir = 1;
								//east_vec[count] = easting[b] + 0.5*dem_res;
								//north_vec[count] = northing[a] + yo - 0.5*dem_res;
								++b;
							}
							else if (degs_old > 180 || degs_new < 360) {
								xo = 0, yo = 1-0.00001;
								s_edge = abs(s_local*sin(theta));
								d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
								xi = 1-xo, yi = yo;
								dir = 3;
								//east_vec[count] = easting[b] -0.5*dem_res;
								//north_vec[count] = northing[a] + yo - 0.5*dem_res;
								--b;

							}
							else {
								cout << "Flow unable to route E or W" << endl;
								exit(EXIT_FAILURE);
							}
						}
						else if (dir == 3) {
							if 	(degs_old <= 270 || degs_new >= 90) {
								xo = 1-0.00001, yo = 0;
								s_edge = abs(s_local*sin(theta));
								d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
								xi = xo, yi = 1-yo;
								dir = 2;
								//east_vec[count] = easting[b] + xo - 0.5*dem_res;
								//north_vec[count] = northing[a] - 0.5*dem_res;
								++a;
							}
							else if (degs_old > 270 || degs_new < 90) {
								xo = 1-0.00001, yo = 1;
								s_edge = abs(s_local*sin((2/PI) - theta));
								d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
								xi = xo, yi = 1- yo;
								dir = 4;
								//east_vec[count] = easting[b] + xo - 0.5*dem_res;
								//north_vec[count] = northing[a] + 0.5*dem_res;
								--a;
							}
							else {
								cout << "Flow unable to route N or S" << endl;
								exit(EXIT_FAILURE);
							}
						}
						else if (dir == 4) {
							if 	(degs_old <= 360 || degs_new >= 180) {
								xo = 0, yo = 0.00001;
								s_edge = abs(s_local*sin((PI/2) - theta));
								d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
								xi = 1-xo, yi = yo;
								dir = 3;
								//east_vec[count] = easting[b] -0.5*dem_res;
								//north_vec[count] = northing[a] + yo - 0.5*dem_res;
								--b;
							}
							else if (degs_old > 0 || degs_new < 180) {
								xo = 1, yo = 0.00001;
								s_edge = abs(s_local*sin(theta));
								d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
								xi = 1-xo, yi = yo;
								dir = 1;
								//east_vec[count] = easting[b] + 0.5*dem_res;
								//north_vec[count] = northing[a] + yo - 0.5*dem_res;
								++b;
							}
							else {
								cout << "Flow unable to route E or W" << endl;
								exit(EXIT_FAILURE);
							}
						}
						slope_total += s_edge*d;


					}
					length += d;
					degs = degs_new;

					//cout << "[a][b]: " << a << " " << b << endl;

					if (a == 0 || b == 0 ||	a == nrows-1 || b == ncols-1 || stnet[a][b] != ndv || path[a][b] == 1) flag = false;
				}

				//if trace finished at a stream, print hillslope info.
				if (stnet[a][b] != ndv) {

					// PRINT TO FILE Cht Sbar Relief Lh
					X = xmin + j*dem_res;
					Y = ymin + (nrows-i)*dem_res;
					relief = zeta[i][j] - zeta[a][b];
					length = length*dem_res;
					mean_slope = slope_total/(length/dem_res);

					ofs << X << " " << Y << " " << seg[i][j] << " "
						<< cht[i][j] << " " << mean_slope << " "
						<< relief << " " << length << " " << "/n"; //area[a][b] << "\n";

					//PRINT FILE OF PATH NODES FOR EACH HILLSLOPE VECTOR
					//stringstream s;
					//s << ht_count;
					//file_part_2 = s.str();
					//filename = file_part_1;
					//filename.append(file_part_2), filename.append(file_part_3);
					//strcpy(filename_c,filename.c_str());

					//ofstream prof_out;
					//prof_out.open(filename_c);
					//prof_out << "Easting " << "Northing" << endl;
					//prof_out.precision(10);

					//for (int c=0;c<count;++c) {
						//print
					//	prof_out 	<< east_vec[c] << " "
					//				<< north_vec[c] << endl;
					//}
					//prof_out.close();
				}
			}
		}
	}
	ofs.close();
}



