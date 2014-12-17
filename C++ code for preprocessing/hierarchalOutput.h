#ifndef hierarchalOutput_H
#define hierarchalOutput_H

#include <cstdlib>
#include <vector>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\imgproc\imgproc.hpp>
#include "structures.h"
#include <cmath>

void hierarchalOutput(std::vector<Ranked_cell> clustered_cells)
{
	std::vector<Compared_row> compared_cells =  clusterSample(clustered_cells);

	std::vector<Compared_row> clustered_rows = clusterReturn(compared_cells, compared_cells.size()-1);
	std::vector<Ranked_cell> output;

	int row_number = clustered_rows.size();
	int cols = clustered_cells[0].cluster_index.size();

	for(int i =0; i<row_number; i++)
	{
		output.push_back(clustered_cells[clustered_rows[i].index]);
	}
	std::ofstream output_file ("clustered sample population single cells2.txt");
	for(int i =0; i<row_number; i++)
	{
		for(int j=0; j<cols; j++)
		{
		output_file << output[i].cluster_index[j] << ',';
		}
		output_file << '\n';
	}
	std::cout << "Hurrah";	
	output_file.close();
}

#endif