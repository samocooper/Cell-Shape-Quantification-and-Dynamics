#ifndef sortAndRank_H
#define sortAndRank_H

#include <math.h>
#include <fstream>
#include <string.h>
#include <cstdlib>
#include <vector>
#include "structures.h"
#include <algorithm>

bool greatest(Ranked_cell i, Ranked_cell j)
{
	return (i.feature_values[Current_feature]<j.feature_values[Current_feature]);
}

std::vector<Ranked_cell> sortRanked(std::vector<Ranked_cell> input_vector)
{
	int cols = input_vector[0].feature_values.size();
	Current_feature = 0;

	for(int h=0; h<cols; h++)
	{
	std::sort(input_vector.begin(),input_vector.end(),greatest);
	for(int i =0; i<input_vector.size(); i++)
	{
		input_vector[i].feature_ranks.push_back(i);
	}
	std::cout << '&';
	Current_feature += 1;
	}
	return input_vector;
}

bool greatest_rank(Ranked_cell i, Ranked_cell j)
{
	return (i.feature_ranks[Current_feature]<j.feature_ranks[Current_feature]);
}

std::vector<std::vector<double>> sortAndCorrelate(std::vector<Ranked_cell> input_vector)
{
	int cols = input_vector[0].feature_ranks.size();
	double size = input_vector.size();
	Current_feature = 0;
	std::vector<std::vector<double>> heatmap;
	for(int h=0; h<cols; h++)
	{
	std::sort(input_vector.begin(),input_vector.end(),greatest);	
	std::vector<double> spearmans (cols,0);

	for(int i =0; i<input_vector.size(); i++)
	{
		int temp_val = input_vector[i].feature_ranks[Current_feature];

		for(int j=0; j<cols; j++)
		{
			float temp = temp_val-input_vector[i].feature_ranks[j];
			spearmans[j] += (double)(temp*temp);
		}
	}
	for(int j=0; j<cols; j++)
	{
		spearmans[j] = 1-(6*spearmans[j]/(size*(size*size)-1));
	}
	heatmap.push_back(spearmans);
	std::cout << '!';
	Current_feature += 1;
	}
	return heatmap;
}
#endif