#ifndef plateKeyReader_H
#define plateKeyReader_H

#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include "structures.h"
#include "fileFinder.h"

std::vector<Well_name> plateKeyRead(std::string fileName)
{
	std::ifstream fileReadout;
	fileReadout.open(fileName);	
	std::vector<Well_name> plate_key;
	Well_name temp_name;
	temp_name.condition="";
	temp_name.well="";	
	bool is_well = true;
	while(fileReadout.good())
	{
		char c = fileReadout.get();				
		if(c!=',' && is_well && c!= '\n' && c !=' ')
		{
			temp_name.well.push_back(c);
		}
		if(c!=',' && !is_well && c!= '\n' && c !=' ')
		{
			temp_name.condition.push_back(c);
		}
		if(c == ',')
		{
			is_well = false;
		}
		if(c == '\n')
		{
			is_well = true;
			plate_key.push_back(temp_name);
			temp_name.condition="";
			temp_name.well="";
		}
	}
	return plate_key;
}
std::vector<Well_name> findPlateKey (void)
{
	std::string input_3 = "";	
	std::string input_4 = "";

	std::cout << "Please enter a sequence unique to the plate key" << '\n';
	std::getline(std::cin, input_3);
	std::cout << "Please enter another sequence unique to the plate key" << '\n';
	std::getline(std::cin, input_4);
	
	std::vector<std::string> temp_plate_key = fileFinder(input_3,input_4);
	std::vector<Well_name> plate_key = plateKeyRead(temp_plate_key[0]);
	for(int i=0; i<plate_key.size(); i++)
	{
		std::cout << plate_key[i].condition << ' ' << plate_key[i].well << '\n';
	}
	return plate_key;
}

#endif