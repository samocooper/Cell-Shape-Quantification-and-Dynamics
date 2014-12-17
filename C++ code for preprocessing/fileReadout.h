#ifndef fileReadout_H
#define fileReadout_H

#include <math.h>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include "structures.h"

std::vector<std::vector<Cell>> Cell_data (std::string fileName)
{

	// open file with name stored in list
	
	std::ifstream PC_data;
	PC_data.open(fileName);

	// create a vetor for each well: size = number of timepoints the sample was imaged for

	unsigned short int number_of_timepoints = 10; // add one	
	std::vector<std::vector<Cell>> well_list(number_of_timepoints);

	if(PC_data.is_open())
	{
		std::cout << '#';

		// ignore opening line containing headings
		PC_data.ignore(800,'\n');

		// counts the number of commas that have passed i.e. stor data after each comma
				
		bool comma = false;
		unsigned short int j =0;

		// intitialise variables to store numbers as they are read

		unsigned short int X_position = 0;
		unsigned short int Y_position = 0;
		unsigned short int T_position =0;
		unsigned short int F_position = 0;	
		unsigned short int prev_timepoint = 0;
		float temp_feature[Number_of_features] = {};
		float temp_place[Number_of_features];
		bool decimal[Number_of_features] = {};	

		for(unsigned char i =0; i<Number_of_features; i++)
		{
			 temp_place[i] = 10;
		}

		std::vector<Cell> time_list;

		bool first = false;

		// convert csv file into an array of Cell objects

		while(PC_data.good())
		{

			char c = PC_data.get();

			if(j == 9)
			{
				if( c!= ',')
				{
					F_position = F_position*10;
					F_position += c - '0';	
				}					
			}
			if(j == 11)
			{
				if( c!= ',')
				{
					T_position = T_position*10;
					T_position += c - '0';				
				}
				if( c == ',' && prev_timepoint != T_position && T_position < number_of_timepoints)
				{
					prev_timepoint = T_position;
					well_list[T_position] = time_list;
					first = true;
					time_list.clear();
				}
			}
			if(j == 13 && c!= ',')
			{
				X_position = X_position*10;
				X_position +=  c - '0';					
			}
			if(j == 14 && c!= ',')
			{
				Y_position = Y_position*10;
				Y_position +=  c - '0';					
			}
			if(j >= 16 && c!= ',' && c!='\n')
			{
				if(c == '.')
				{
					decimal[j-16] = true;
				}
				if(!decimal[j-16] && c != '.')
				{
					temp_feature[j-16] = temp_feature[j-16]*10;
					temp_feature[j-16] +=  c - '0';
				}
				if(decimal[j-16] && c != '.')
				{
					temp_feature[j-16] += ((float)(c - '0'))/temp_place[j-16];
					temp_place[j-16] = temp_place[j-16]*10;
				}
			}
			if (j<9)
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
				Cell temp_cell = {};
				
				temp_cell.X_pos = X_position;
				temp_cell.Y_pos = Y_position;
				temp_cell.field = F_position;

				for(int k =0; k<Number_of_features; k++)
				{
				temp_cell.Features[k] = temp_feature[k];				
				temp_feature[k] =0;
				temp_place[k] = 10;
				decimal[k] = false;
				}				

				time_list.push_back(temp_cell);

				if(temp_cell.field > 9)
				{
					time_list.pop_back();
				}
				X_position = 0;
				Y_position = 0;
				T_position =0;			
				F_position =0;

				j =0;							
			}
		}			
	}	
	PC_data.close();
	return well_list;
}

#endif