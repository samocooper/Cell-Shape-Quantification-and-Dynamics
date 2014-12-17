#ifndef fileReader_H
#define fileReader_H

#include <math.h>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>

#include "structures.h"

std::vector<std::vector<Cell_static>> Cell_data_static (std::string fileName)
{

	// open file with name stored in list
	
	std::ifstream PC_data;
	PC_data.open(fileName);
	std::vector<std::vector<Cell_static>> well_list;

	if(PC_data.is_open())
	{
		std::cout << '#';

		PC_data.ignore(6000,'\n');	

		bool comma = false;
		unsigned short int j =0;

		std::string temp_well_name = "";
		std::string prev_well_name = "";
		int temp_well_index = 0;
		
		float temp_feature[Number_of_features] = {};
		float temp_place[Number_of_features];
		bool decimal[Number_of_features] = {};	

		for(unsigned char i =0; i<Number_of_features; i++)
		{
			 temp_place[i] = 10;
		}

		std::vector<Cell_static> time_list;
		
		// convert csv file into an array of Cell objects

		while(PC_data.good())
		{

			char c = PC_data.get();
			if(j == 6)
			{
				if(c != ',')
				{
					temp_well_name.push_back(c);
				}
			}
			if(j >= 15 && c!= ',' && c!='\n' && j < 61)
			{
				if(c == '.')
				{
					decimal[j-15] = true;
				}
				if(!decimal[j-15] && c != '.')
				{
					temp_feature[j-15] = temp_feature[j-15]*10;
					temp_feature[j-15] +=  c - '0';
				}
				if(decimal[j-15] && c != '.')
				{
					temp_feature[j-15] += ((float)(c - '0'))/temp_place[j-15];
					temp_place[j-15] = temp_place[j-15]*10;
				}
			}
			if (j<6)
			{
				PC_data.ignore(500,',');
				j += 1;
			}
			if( c == ',')
			{
				j += 1;
			}					
			if(c == '\n')
			{
				Cell_static temp_cell;
				for(int k =0; k<Number_of_features; k++)
				{
				temp_cell.Features[k] = temp_feature[k];				
				temp_feature[k] =0;
				temp_place[k] = 10;
				decimal[k] = false;
				}
				temp_cell.well_name = temp_well_name;
				
				if(temp_well_name.compare(prev_well_name) !=0)
				{
					prev_well_name = temp_well_name;
					temp_well_index +=1;	
					if(time_list.size()>20)
					{
					well_list.push_back(time_list);
					}
					time_list.clear();
				}
				temp_cell.well_index = temp_well_index;				
				time_list.push_back(temp_cell);
				temp_well_name = "";
				j =0;
			}
		}		
		well_list.push_back(time_list);
	}
	PC_data.close();
	return well_list;
}

#endif