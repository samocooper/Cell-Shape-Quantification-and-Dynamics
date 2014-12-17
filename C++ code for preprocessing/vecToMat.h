#ifndef vecToMat_H
#define vecToMat_H

#include <cstdlib>
#include <vector>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\imgproc\imgproc.hpp>
#include "structures.h"

cv::Mat vecToMat (std::vector<Ranked_cell> input_cells)
{
	int size = input_cells.size();
	int features = input_cells[0].feature_values.size();
	cv::Mat input_mat (size,features,CV_32FC1);

	float *p;
	for(int i = 0; i<size; i++)
	{
		p =  input_mat.ptr<float>(i);
		for (int j=0; j<features; j++)
		{
			p[j] =  input_cells[i].feature_values[j];
		}
	}
	return input_mat;
}

cv::Mat int_vecToMat (std::vector<Ranked_cell> input_cells)
{
	int size = input_cells.size();
	int features = input_cells[0].cluster_index.size();
	cv::Mat input_mat (size,features,CV_32FC1);

	float *p;
	for(int i = 0; i<size; i++)
	{
		p =  input_mat.ptr<float>(i);
		for (int j=0; j<features; j++)
		{
			p[j] =  (float)input_cells[i].cluster_index[j];
		}
	}
	return input_mat;
}


#endif