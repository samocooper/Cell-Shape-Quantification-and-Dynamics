#ifndef staticCombine_H
#define staticCombine_H

#include <vector>
#include "structures.h"
#include <iostream>

std::vector<Ranked_cell> staticCombine(std::vector<std::vector<std::vector<Cell_static>>> file_read, int max_cells)
{
	std::vector<Ranked_cell> all_cells;
	std::vector<std::vector<Cell_static>> experiment = file_read.back();
	int experiment_size = experiment.size();
	for(int i=0; i<experiment_size; i++)
	{
		std::vector<Cell_static> well = experiment[i];
		int temp_size = max_cells;
		if (well.size()<max_cells)
		{
			temp_size = well.size();
		}
		for(int i=0; i<temp_size; i++)
		{

			Ranked_cell temp_ranked;
			Cell_static temp_cell = well[i];
			for(int j=0; j <Number_of_features; j++)
			{
				temp_ranked.feature_values.push_back(temp_cell.Features[j]);
			}
			all_cells.push_back(temp_ranked);
		}
		std::cout << "@";
	}
	return all_cells;
}

#endif