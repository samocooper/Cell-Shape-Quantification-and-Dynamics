#ifndef fileFinder_H
#define fileFinder_H
 

#include <string.h>
#include <windows.h>
#include <cstdlib>
#include <vector>

std::vector<std::string> fileFinder(std::string identifier , std::string identifier2)
{
	WIN32_FIND_DATA FindData;
	HANDLE hFind;

	std::vector<std::string> fileList;
	hFind = FindFirstFile("columbus/*.csv", &FindData);
	
	do{
		std::string fileName = "columbus/";
		fileName.append(FindData.cFileName);
		size_t found = fileName.find(identifier);
		size_t found2 = fileName.find(identifier2);
		if(found != std::string::npos && found2 != std::string::npos)
		{
			fileList.push_back(fileName);
		}
	}while (FindNextFile(hFind, &FindData));

	return fileList;
}
 
#endif