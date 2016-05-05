//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Cht_Res.cpp
//
// Driver created to generate hilltop curvature data at a range of input resolutions.
//
// Driver expects unfilled DEMs at a range of grid resolutions in the input directory
//
// Run the driver with the following arguments:
//
// path to the DEM files with a trailing slash
// DEM filename Prefix
// file extension without the dot
// Window size
//
// DEM naming convention should be <prefix>_<resolution>_DEM
//
// A usage example is:
// nice ./CurvatureRes.out /home/s0675405/DataStore/Data/GM/ GM bil 5
//
// Output data is stored in the input directory in files called <Prefix>_ChtResData_<curvature type>.txt
// and <Prefix>_<Resolution>_Hist_CHT.txt
// This data can be plotted using python code found at https://github.com/sgrieve/Resolution_Paper_Figs/
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Stuart W.D. Grieve
// University of Edinburgh
// November 2015
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDIndexRaster.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDChannel.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDIndexChannel.hpp"
#include "../LSDMostLikelyPartitionsFinder.hpp"
#include "../LSDBasin.hpp"
#include "../LSDShapeTools.hpp"

int main(int nNumberofArgs, char *argv[])
{
  //Test for correct input arguments
	if (nNumberofArgs!=5)
	{
		cout << "FATAL ERROR: wrong number of inputs. The program needs the path (with trailing slash), the filename prefix, the DEM file format, the window size in spatial units.";
		exit(EXIT_FAILURE);
	}

  //get input args
  string path = argv[1];
  string Prefix = argv[2];
  string DEM_Format = argv[3];
  int WindowSize = atoi(argv[4]);

  //surface fitting
  vector<int> raster_selection;
  raster_selection.push_back(0);
  raster_selection.push_back(1); //slope
  raster_selection.push_back(0);
	raster_selection.push_back(1); //curvature

  //set up a writer to write the output data
  ofstream WriteData;

  //create an output filename based on the dem name
  stringstream ss;

  ss << path << Prefix << "_ChtResData.txt";
  WriteData.open(ss.str().c_str());

  //write headers
  WriteData << "resoulution 2pc 25pc median mean 75pc 98pc minimum maximum" << endl;

  //array of resolutions to load
  int Resolutions[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

  // vectors to hold the stats about the fitted surface
  vector<float> Curv_vec;

	//set boundary conditions
  vector<string> BoundaryConditions(4, "No Flux");

  for (int a = 0; a < 10; ++a){

	  cout << "Processing DEM " << a+1 << " of " << "10" << endl;

		//load the DEM
    //build the string of the filename to load
		stringstream ss2;
		ss2 << path << Prefix << "_" << Resolutions[a] << "_DEM";
		LSDRaster DEM(ss2.str(), DEM_Format);

		//Fill
		float MinSlope = 0.0001;
		LSDRaster FilledDEM = DEM.fill(MinSlope);

	  int CurrentWindowSize = WindowSize;
	  vector<LSDRaster> Surfaces = FilledDEM.calculate_polyfit_surface_metrics(CurrentWindowSize, raster_selection);

		// get a flow info object
	  LSDFlowInfo FlowInfo(BoundaryConditions,FilledDEM);

	  //get stream net from channel heads
	  vector<int> sources = FlowInfo.Ingest_Channel_Heads((path+Prefix+"_CH"), "csv", 2);

	  LSDJunctionNetwork ChanNetwork(sources, FlowInfo);

		// extract hilltops
    LSDRaster hilltops = ChanNetwork.ExtractRidges(FlowInfo);
    LSDRaster Hilltops = ChanNetwork.ExtractHilltops(hilltops, Surfaces[1], 0.4);


    //get hilltop curvature using filter to remove positive curvatures
    LSDRaster cht_raster = FilledDEM.get_hilltop_curvature(Surfaces[3], Hilltops);
    LSDRaster CHT = FilledDEM.remove_positive_hilltop_curvature(cht_raster);

	  //go through the landscape and get every curvature value into a 1D vector
    Curv_vec = Flatten_Without_Nodata(CHT.get_RasterData(), CHT.get_NoDataValue());

    stringstream ss3;
    ss3 << path << Prefix << "_" << Resolutions[a] << "_Hist_CHT.txt";

    print_histogram(Curv_vec, 0.01, ss3.str());

    vector<float> Boxplot = BoxPlot(Curv_vec);

	  //write the values to the output file
	  WriteData << Resolutions[a] << " " << Boxplot[0] << " " << Boxplot[1] << " " << Boxplot[2] << " " << Boxplot[3] << " " << Boxplot[4] << " " << Boxplot[5] << " " << Boxplot[6] << " " << Boxplot[7];

		WriteData << endl;

  }
  WriteData.close();

}
