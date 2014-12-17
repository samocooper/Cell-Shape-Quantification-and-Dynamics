#ifndef hierarchalCluster_H
#define hierarchalCluster_H

#include <cstdlib>
#include <vector>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\imgproc\imgproc.hpp>
#include "structures.h"
#include <cmath>

std::vector <Compared_row> clusterSample ( std::vector<Ranked_cell> input_cells )
{
	int rows = input_cells.size();
	int cols = input_cells[0].cluster_index.size();
	std::vector<Compared_row> compared;

	for(int i =0; i<rows; i++)
	{
		Compared_row compared_row;
		Ranked_cell temp_cell = input_cells[i];
		for(int j = 0; j<rows; j++)
		{
			Ranked_cell compare_cell = input_cells[j];
			uchar value = 0;
			for(int k=0; k<cols; k++)
			{
				if(compare_cell.cluster_index[k] == temp_cell.cluster_index[k])
				{
					value++;
				}
			}
			compared_row.value.push_back(value);
		}		
		compared_row.integrated = false;
		compared_row.is_last = true;
		compared_row.index = i;
		compared.push_back(compared_row);
	}
	
	for(int h=0; h<rows-1; h++)
	{
		float max =0;
		int maxi;
		int maxj;
		for(int i =0; i<compared.size(); i++)
		{
		
			Compared_row compared_row = compared[i];		
			if(!compared_row.integrated)
			{
				for(int j =0; j<rows; j++)
				{
					uchar temp_val = compared_row.value[j];
					if(temp_val >= max && i != j)
					{
						max = temp_val;
						maxi = i;
						maxj = j;
					}
				}
			}
		}
		std::cout << max;
		for(int i =0; i<rows; i++)
		{
			if(compared[maxj].value[i] == max)
			{
				compared[maxj].value[i] = 0;				
			}
			if(compared[maxi].value[i] == max)
			{
				compared[maxi].value[i] = 0;
			}
		}
		while(compared[maxj].integrated)
		{
			maxj = compared[maxj].next;
		}
		compared[maxi].integrated = true;
		compared[maxj].integrated = true;
		Compared_row combined;

		compared[maxi].next = compared.size();
		compared[maxj].next = compared.size();

		combined.left = maxi;
		combined.right = maxj;

		for(int i =0; i<rows; i++)
		{		
			combined.value.push_back(compared[maxj].value[i]);
			if(compared[maxi].value[i]>combined.value[i])
			{
				combined.value[i] = compared[maxi].value[i];
			}
			if(i == maxi || i == maxj || compared[maxi].value[i] ==0 || compared[maxj].value[i] ==0)
			{
				combined.value[i] =0;
			}
		}
		combined.is_last = false;
		combined.integrated = false;
		compared.push_back(combined);
	}
	return compared;
}
std::vector<Compared_row> clusterReturn(std::vector<Compared_row> input, int index)
{
	std::vector<Compared_row> list;
	Compared_row node = input[index];
	if(node.is_last)
	{
		list.push_back(node);
	}
	if(!node.is_last)
	{
		std::vector<Compared_row> leftlist = clusterReturn(input, node.left);
	int lsize = leftlist.size();
	for(int i=0; i<lsize; i++)
	{
		list.push_back(leftlist[i]);
	}
	std::vector<Compared_row> rightlist =clusterReturn(input, node.right);
	int rsize=rightlist.size();
	for(int i=0; i<rsize; i++)
	{
		list.push_back(rightlist[i]);
	}
	}
	return list;
}
#endif
			

