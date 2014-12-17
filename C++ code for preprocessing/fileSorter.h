#ifndef fileSorter_H
#define fileSorter_H

#include <math.h>
#include <fstream>
#include <string.h>
#include <cstdlib>
#include <vector>
#include "structures.h"

std::vector<Tracked_cell> dataSort(std::vector<std::vector<Cell>> Well)
{
	signed int short well_size = Well.size();
	std::vector<Tracked_cell> Feature_result;
	std::vector<std::vector<std::vector<Cell>>> Time_Field_list;

	//Intitialise the vector with the number of fields of view used in the experiment and the number of timepoints specified

	for(int i = 0; i<9; i++) // i<number_of_fields
	{
		std::vector<std::vector<Cell>> Field_temp ;			
		for(int j =0; j < well_size; j++)
		{
			std::vector<Cell> At_time;
			Field_temp.push_back(At_time);
		}
		Time_Field_list.push_back(Field_temp);
	}

	// Use the field variable of the cell object to assign the cell object to it's field

	for(int i =0; i<well_size; i++)
	{
		std::vector<Cell> Time_point = Well[i];
		while(!Time_point.empty())
		{
			Cell Temp_cell = Time_point.back();
			Time_point.pop_back();
			Time_Field_list[Temp_cell.field-1][i].push_back(Temp_cell);
		}
	}	
	signed int short time_field_list_size = Time_Field_list.size();
	for(int i=0; i<time_field_list_size; i++)
	{

		// Use the X and Y positions to track the cells 
		signed int short field_list = Time_Field_list[i].size()-1;
		for(int j=0; j < field_list; j++)
		{
			signed int short list_1 = Time_Field_list[i][j+1].size();
			for(int k=1; k< list_1; k++)
			{
				int x1 = Time_Field_list[i][j+1][k].X_pos;
				int y1 = Time_Field_list[i][j+1][k].Y_pos;	
				Time_Field_list[i][j+1][k].field = -1;
				float distance = 100000;
				int list_2 = Time_Field_list[i][j].size();
				for(int l=0; l<list_2; l++)
				{
					int X_dist = x1-Time_Field_list[i][j][l].X_pos;
					int Y_dist = y1-Time_Field_list[i][j][l].Y_pos;
					int temp_dist = X_dist*X_dist +Y_dist*Y_dist;
					if(distance > temp_dist)
					{
						distance = (float)temp_dist;
						Time_Field_list[i][j+1][k].field = l;
					}						
				}

				for(int l=0; l<k; l++)
				{								
					if(Time_Field_list[i][j+1][k].field == Time_Field_list[i][j+1][l].field && Time_Field_list[i][j+1][k].field >0)
					{
						float temp_dist2 = Time_Field_list[i][j+1][l].X_pos*Time_Field_list[i][j][Time_Field_list[i][j+1][k].field].X_pos + Time_Field_list[i][j+1][l].Y_pos*Time_Field_list[i][j][Time_Field_list[i][j+1][k].field].Y_pos;
						if (distance > temp_dist2)
						{
							Time_Field_list[i][j+1][k].field = -1;
						}	
						if (distance < temp_dist2)
						{
							Time_Field_list[i][j+1][l].field = -1;
						}
					}
				}
				if(Time_Field_list[i][j+1][k].field >-1)
				{
					Time_Field_list[i][j+1][k].moved = sqrt(distance);
					for(int m=0; m<Number_of_features; m++)
					{
						Time_Field_list[i][j+1][k].Features_dynamic[m] = Time_Field_list[i][j+1][k].Features[m] - Time_Field_list[i][j][Time_Field_list[i][j+1][k].field].Features[m];
					}
				}
			}
		}

		// Use the X and Y positions to track the cells 

		for(int j = field_list-1; j>=0; j--)
		{
			signed int short list_1 = Time_Field_list[i][j].size();
			for(int l=0; l<list_1; l++)
			{
				int next_cell = l;
				Tracked_cell Temp_tracked = {};
				int Time_position = j;					
				int count = 0;				
				int previous_cell = -1;
				int previous_movement[Number_of_features] = {};				

				if(Time_Field_list[i][Time_position][l].tracked == false)
				{
					while(Time_position >= 0 && next_cell > -1 && next_cell<Time_Field_list[i][Time_position].size())
					{

						// Break loop if a cell has moved beyond 20 pixels - this suggests a different cell is being tracekd

						if(Time_Field_list[i][Time_position][next_cell].moved>20)
						{
							break;
						}						

						Time_Field_list[i][Time_position][next_cell].tracked = true;
						Temp_tracked.X_values.push_back(Time_Field_list[i][Time_position][next_cell].X_pos);
						Temp_tracked.Y_values.push_back(Time_Field_list[i][Time_position][next_cell].Y_pos);
						Temp_tracked.moved.push_back(Time_Field_list[i][Time_position][next_cell].moved);
						
						std::vector<float> features;
						std::vector<float> dynamics;
						for(int m=0; m<Number_of_features; m++)
						{							
							features.push_back(Time_Field_list[i][Time_position][next_cell].Features[m]);				
							dynamics.push_back(fabs(Time_Field_list[i][Time_position][next_cell].Features_dynamic[m]));
						}
						Temp_tracked.dynamic_values.push_back(dynamics);						
						Temp_tracked.feature_values.push_back(features);

						Temp_tracked.frames_tracked += 1;
						previous_cell = next_cell;
						count ++;
						next_cell = Time_Field_list[i][Time_position][next_cell].field;	
						Time_position -= 1;						
					}
				}
				if(Temp_tracked.frames_tracked > 2)
				{
					Feature_result.push_back(Temp_tracked);
				}
				
			}			
		}
	}
	std::cout << Feature_result.size() << ' ';
	return Feature_result;
}
// Here you can manually determine which features should be used within the clustering proccess
std::vector<Ranked_cell> dynamicCombine(std::vector<TC> sorted_readout)
{
	std::vector<Ranked_cell> all_cells;
	int max_frames = 1; // must be less than cut off point for number of frames stored
	int temp_size = sorted_readout.size();
	for(int g =0; g<temp_size; g++)
	{
		int max_cells = 10;
		TC temp_tc = sorted_readout[g];
		if(temp_tc.cells.size()<max_cells)
		{
			max_cells = temp_tc.cells.size();
		}
		for(int h =0; h< max_cells; h++)
		{
			Tracked_cell temp_cell = temp_tc.cells[h];			
			for(int i =0; i<max_frames; i++)
			{
				Ranked_cell temp_ranked;

				for( int j =0; j<Number_of_features; j++)
				{
					if(j!=10)
					{
						temp_ranked.feature_values.push_back(temp_cell.feature_values[i][j]);
						temp_ranked.dynamic_values.push_back(temp_cell.dynamic_values[i][j]);						
					}
				}
				temp_ranked.dynamic_values.push_back((float)temp_cell.moved[i]);
				all_cells.push_back(temp_ranked);
			}
		}
		std::cout << '@';
	}
	return all_cells;
}
#endif