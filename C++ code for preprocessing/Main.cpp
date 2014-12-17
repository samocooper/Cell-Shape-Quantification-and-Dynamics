#include <math.h>
#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include <Windows.h>
#include <opencv2\ml\ml.hpp>
#include <opencv2\imgproc\imgproc.hpp>
#include <opencv2\highgui\highgui.hpp>

#include "fileFinder.h"
#include "structures.h"
#include "fileSorter.h"
#include "inputReader.h"
#include "fileReadout.h"
#include "plateKeyReader.h"
#include "clusterData.h"
#include "vecToMat.h"
#include "hierarchalCluster.h"
#include "hierarchalOutput.h"
#include "fileReader.h"
#include "staticCombine.h"
#include "gapStatistic.h"
#include "testReader.h"

using namespace std;
using namespace cv;

bool down = false;
int prescision = 50;

static void onMouse( int event, int x, int y, int, void* )
{	
	if( event == EVENT_LBUTTONDOWN )
	{
		down = true;
		return;
	}
}

int main()
{

	// Use two user inputted strings to identify files within folder

	string input_1 = "";
	string input_2 = "";

	cout << "Please enter a sequence unique to the file names" << '\n';
	getline(cin, input_1);
	cout << "Please enter another sequence unique to the file names" << '\n';
	getline(cin, input_2);
	vector<string> dataFiles = fileFinder(input_1,input_2);

	// Find the plate key for the files only nessacary for dynamic data

	vector<Well_name> plate_key = findPlateKey();

	//Vector to hold the data on cells from either the dynamic or static file. 
	//In the case of "dynamic files" vector<vector<Cell>> is a well (one csv file) vector<Cell> is a timepoint, each Cell has information on the field of view and the X Y position within that field.
	//In static vector<vector<Cell_static>> is a single file, vector<Cell_static> is a single well. Currently only wells with over 20 cells detected are included. 
	
		//----------------STATIC-----------------
	//vector<vector<vector<Cell_static>>> file_readout;

		//----------------DYNAMIC-----------------
	vector<vector<vector<Cell>>> file_readout;
	
	// Read the data from the files found and return a list of Cell objects
		
	int total_number_of_tcs =  dataFiles.size();
	for( int i = 0; i<total_number_of_tcs; i++)
	{
			//----------------STATIC-----------------
		//file_readout.push_back(Cell_data_static(dataFiles[i]));	

			//----------------DYNAMIC-----------------
		file_readout.push_back(Cell_data(dataFiles[i]));	
	}
	cout << '\n' << "Number of files found:  "<< file_readout.size() << '\n';

	// Number of cells per well to be used in the overall mixture of cells to be analysed
	// Here a sample population is created to train the models upon or determine the gap statistic for K means

		//--------------TEST---------------------------
	//vector<Ranked_cell> all_cells = testRead(dataFiles[0]);
	//int max_cells = 30;

		//----------------STATIC-----------------

	//vector<Ranked_cell> all_cells = staticCombine(file_readout,max_cells);
	
	// When data is dynamic the cells need to be tracked over timepoints.
	// This is done through the X,Y position and the field of view being a) the closest b) within a cutoff distance to filter out incorrectly segmented cells.
	// Only cells which have been tracked for over 3 timepoints are returned as "Tracked Cell" objects, this seems to greatly improve the quality of data.

		//----------------DYNAMIC-----------------
	
	vector<TC> sorted_readout;

	for( int i = 0; i<total_number_of_tcs; i++)
	{
		TC temp_tc;
		for(int j=0; j<total_number_of_tcs;j++)
		{
			size_t found = dataFiles[i].find(plate_key[j].well);
			if(found != std::string::npos)
			{
				temp_tc.condition = plate_key[j].condition;				
				temp_tc.well = plate_key[j].well;
			}
		}
		temp_tc.cells = dataSort(file_readout[i]);
		sorted_readout.push_back(temp_tc);
	}
		//----------------DYNAMIC-----------------

	vector<Ranked_cell> all_cells = dynamicCombine (sorted_readout);

	cout << '\n' << "Size of sample population:" << all_cells.size() << '\n';

	// Convert the the vector of "Ranked" cells into a matrix, where each row of the matrix corresponds to a ranked cell.
	// This is necessary since the OpenCV clustering algorithms use a 32 bit floating point matrix to train upon, the type is"CV_32FC1".

	Mat sample = vecToMat(all_cells);

	// Train a vector of EM models, one for each feature, upon the training sample. Currently the number of clusters is either 1 or 2 and their is a fiddle factor, I would like to change this.
	int count =0;
	float X = 0.24;
	//for(float X = 0.16; X<=0.3; X += 0.02)
	//{
	cout << "REPEAT DONE ";

	vector<EM> models = trainModels(sample,X);

	// Open a file to output data to.

	ofstream output("K = 5 final with dynamics.txt");

	// Use the trained models to assign for each feature the cluster to which it belongs i.e. the cluster_index member of Ranked_cell.

	vector<Ranked_cell> clustered_cells;
	clustered_cells = cellProfile(all_cells,models);

	int sample_size = clustered_cells.size();
	int cols = clustered_cells[0].cluster_index.size();

	// Create a 32_FC1 matrix, where each row corresponds to the cluster_indexs of a Ranked_cell.

	Mat Kmeans_clustered = int_vecToMat(clustered_cells);
	//Optionally determine the Gap statistic.
	//gapStatisticOutput(Kmeans_clustered,count,10);
	count++;
	//}
	
	Mat Kmeans_index;
	Mat Kmeans_centers;
	cv::TermCriteria Kmeans_parameters;
	Kmeans_parameters.maxCount = 100;

	for(int i=0; i<cols; i++)
	{
		cout << Clusters_for_feature[i] << ' ';
	}
	int K=5;
	
	// The returned value is the quality of clustering, i.e. the mean square distnaces from each center, this is what is plotted in the gap statistic
	
	double compact = kmeans(Kmeans_clustered,K,Kmeans_index,Kmeans_parameters,200,KMEANS_PP_CENTERS,Kmeans_centers);
	int cols1 = sorted_readout[0].cells[0].feature_values.size();

	//Determine TC effect on populations
	/*/
	for( int i = 0; i<total_number_of_tcs; i++)//file_readout[0].size()
	{
		vector<Tracked_cell> population = sorted_readout[i].cells;
		vector<Ranked_cell> ranked_cells;
		int population_size = population.size();

		for(int j=0; j<population_size; j++)
		{
			Tracked_cell temp_cell;
			temp_cell = population[j];			
			while(!temp_cell.feature_values.empty())
			{
				Ranked_cell temp_ranked;
				vector<float> features = temp_cell.feature_values.back();
				for( int k =0; k<cols1; k++)
				{
					if(k!=10)
					{
						temp_ranked.feature_values.push_back(features[k]);
					}
				}							
				ranked_cells.push_back(temp_ranked);
				temp_cell.feature_values.pop_back();	
			}
		}
		/*/
		/*/
		vector<Cell_static> well = file_readout[0][i];
		vector<Ranked_cell> ranked_cells;
		int population_size = well.size();
		for(int j=0; j<population_size; j++)
		{
		Cell_static temp_cell;
		temp_cell = well[j];
		Ranked_cell temp_ranked;
		for( int k =0; k<Number_of_features; k++)
		{
		if(k!=10)
		{
		temp_ranked.feature_values.push_back(temp_cell.Features[k]);
		}
		}
		ranked_cells.push_back(temp_ranked);
		}
		
	/*/
		/*/
		vector<Ranked_cell> profiled_cells = cellProfile(ranked_cells,models);
		int cols = profiled_cells[0].cluster_index.size();
		vector<float> TC_result (K,0);

		for(int j=0; j<population_size; j++)
		{
			float smallest = cols;
			float index;

			for(int k=0; k<K; k++)
			{
				float distance = 0;
				for(int l=0; l<cols; l++)
				{
					distance += fabs(Kmeans_centers.at<float>(k,l)-profiled_cells[j].cluster_index[l]);
				}
				if(distance < smallest)
				{
					smallest = distance;
					index = k;
				}
			}
			TC_result[index] +=1;
		}
		
		float total =0;
		for (int i=0; i<K; i++)
		{
			total += TC_result[i];
		}
		for (int i=0; i<K; i++)
		{
			TC_result[i] = TC_result[i]/total;
			cout << TC_result[i] << ',';
			output << TC_result[i] << ',';
		}
		cout << '\n';
		output << '\n';
	}
	output.close();
	/*/
	/*/
	vector<int> cluster_size (K,0);
	for(int i = 0; i<sample_size; i++)
	{
	cluster_size[Kmeans_index.at<int>(i,0)] +=1;
	}
	for(int j=0; j<K; j++)
	{
	output << cluster_size[j] << ',';
	}
	output << '\n';
	}
	output.close();
	/*/
	//Output overall population clustered
	float *p;
	int dynamics = clustered_cells[0].dynamic_values.size();
	for(int i = 0; i<sample_size; i++)
	{
	output << Kmeans_index.at<int>(i,0) << ',';
	p =  Kmeans_clustered.ptr<float>(i);
	for (int j=0; j<cols; j++)
	{
	output << p[j] << ',';
	}
	for (int j=0; j<dynamics; j++)
	{
		output << clustered_cells[i].dynamic_values[j] << ',';
	}
	output << '\n';
	}

	output.close();
	

	//hierarchalOutput(clustered_cells);
}