#ifndef gapStatistic_H
#define gapStatistic_H

#include <cstdlib>
#include <vector>
#include <math.h>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\imgproc\imgproc.hpp>
#include "structures.h"
#include <fstream>
#include <string>
cv::Mat generateTest (int rows, int cols)
{
	cv::Mat test_mat = cv::Mat::zeros(rows,cols,CV_32FC1);
	cv::RNG rng;

	for(int i=0; i<cols; i++)
	{
		int max_val = Clusters_for_feature[i]+1;
		for(int j =0; j<rows; j++)
		{
			int x = (int)rng.uniform(0,max_val);
			test_mat.at<float>(j,i) += x;
		}
	}
	return test_mat;
}

void gapStatisticOutput (cv::Mat Kmeans_clustered,uchar repeats, int max_clusters)
{
	int rows = Kmeans_clustered.rows;
	int cols = Kmeans_clustered.cols;
	std::string name = "Kmeans clustering qualities Melanoma gtpase screen 1.txt";
	std:: cout << "NAME "<< name;
	name[28] += repeats;
	std::ofstream output(name);
		for(int K =1; K<max_clusters; K++)
		{
			std::cout << "!";
			cv::Mat Kmeans_centers;
			cv::TermCriteria Kmeans_parameters;
			Kmeans_parameters.maxCount = 100;
			cv::Mat Kmeans_index;

			cv::Mat test_mat = generateTest(rows,cols);
			cv::Mat test_index;
			cv::Mat test_centers;

			double compact = kmeans(Kmeans_clustered,K,Kmeans_index,Kmeans_parameters,400,cv::KMEANS_PP_CENTERS,Kmeans_centers);
			double compact_test = kmeans(test_mat,K,test_index,Kmeans_parameters,400,cv::KMEANS_PP_CENTERS,test_centers);

			std::vector<float> within_cluster_distance (K,0);
			for(int l=0; l<rows; l++)
			{

				float smallest = cols;
				int cluster =0;
				for (int m =0; m<K; m++)
				{	
					float distance = 0;
					for(int n=0; n<cols; n++)
					{
						float temp = fabs(Kmeans_centers.at<float>(m,n)-Kmeans_clustered.at<float>(l,n));
						
						distance += temp*temp;
					}
					if (distance <= smallest)
					{
						smallest = distance;
						cluster = m;
					}
				}				
				within_cluster_distance[cluster] += smallest;
			}
			

			double statistic =0;
			for(int l=0; l<K; l++)
			{
				statistic += (double)within_cluster_distance[l]/(double)rows;
			}

			std::vector<float> test_within_cluster_distance (K,0);

			for(int l=0; l<rows; l++)
			{

				float smallest = cols*cols;
				int cluster = 0;
				for (int m =0; m<K; m++)
				{
					float distance = 0;
					for(int n=0; n<cols; n++)
					{
						float temp = fabs(test_centers.at<float>(m,n)-test_mat.at<float>(l,n));
						distance += temp*temp;
					}	
					if (distance <= smallest)
					{
						smallest = distance;
						cluster = m;
					}
				}
				
				//Test output 
				test_within_cluster_distance[cluster] += smallest;
			}

			double test_statistic =0;
			for(int l=0; l<K; l++)
			{
				test_statistic += (double)test_within_cluster_distance[l]/(double)rows;
			}

			statistic = statistic/2;
			test_statistic = test_statistic/2;

			std::cout << statistic << ',' << test_statistic << ',' << test_statistic-statistic << '\n';
			output << statistic << ',' << test_statistic <<  ',' << test_statistic-statistic <<'\n';
	}
	output.close();
}
#endif