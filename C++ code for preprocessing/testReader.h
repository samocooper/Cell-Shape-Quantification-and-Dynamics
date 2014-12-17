#ifndef testReader_H
#define testReader_H

#include <math.h>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>

#include "structures.h"

std::vector<Ranked_cell> testRead (std::string fileName)
{

	// open file with name stored in list
	
	std::ifstream PC_data;
	PC_data.open(fileName);
	std::vector<Ranked_cell> well_list;

	if(PC_data.is_open())
	{
		std::cout << '#';

		PC_data.ignore(6000,'\n');	

		bool comma = false;
		unsigned short int j =0;
		
		float temp_feature[Number_of_features] = {};
		float temp_place[Number_of_features];
		bool decimal[Number_of_features] = {};	

		for(unsigned char i =0; i<Number_of_features; i++)
		{
			 temp_place[i] = 10;
		}
		// convert csv file into an array of Cell objects

		while(PC_data.good())
		{			
			char c = PC_data.get();

			if(j >= 2 && c!= ',' && c!='\n')
			{
				if(c == '.')
				{
					decimal[j-2] = true;
				}
				if(!decimal[j-2] && c != '.')
				{
					temp_feature[j-2] = temp_feature[j-2]*10;
					temp_feature[j-2] +=  c - '0';
				}
				if(decimal[j-2] && c != '.')
				{
					temp_feature[j-2] += ((float)(c - '0'))/temp_place[j-2];
					temp_place[j-2] = temp_place[j-2]*10;
				}
			}			
			if( c == ',')
			{
				j += 1;
			}
			if(c == '\n')
			{
				Ranked_cell temp_cell;
				for(int k =0; k<Number_of_features; k++)
				{
					temp_cell.feature_values.push_back(temp_feature[k]);	
					std::cout << ' ' <<temp_feature[k];
					temp_feature[k] =0;
					temp_place[k] = 10;
					decimal[k] = false;
				}		
				well_list.push_back(temp_cell);
				j =0;
			}			
		}
	}
	PC_data.close();
	return well_list;
}

#endif