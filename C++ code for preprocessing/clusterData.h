#ifndef clusterData_H
#define clusterData_H

#include <cstdlib>
#include <vector>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\imgproc\imgproc.hpp>
#include <opencv2\ml\ml.hpp>
#include "structures.h"
#include <cmath>


std::vector<cv::EM> trainModels (cv::Mat cluster_input,float beta)
{
	cv::TermCriteria parameters;
	parameters.maxCount = 100;
	parameters.type = 3;

	cv::TermCriteria parameters_dud;
	parameters.maxCount =0;
	parameters.type = 3;

	std::vector<std::vector<cv::EM>> all_final_models;
	std::vector<std::vector<double>> sum_of_logs;

	int cols = cluster_input.cols;
	int rows = cluster_input.rows;
	int number_of_clusters = 6;
	for(int g =1; g <=number_of_clusters; g++)
	{		
		std::vector<cv::EM> final_models;
		std::vector<std::vector<cv::Mat>> average_covs;
		std::vector<cv::Mat> average_means;
		std::vector<cv::Mat> average_weights;
		cv::Mat log_likelihoods(rows,1,CV_64FC1);		
		std::vector<double> mean_of_logs;
		// Initialize average matrices and final model

		for(int h =0; h<cols; h++)
		{
			double temp_double = 0;
			mean_of_logs.push_back(temp_double);
			cv::EM temp_model (g,1,parameters_dud);
			cv::Mat dud = cv::Mat::zeros(g,1,CV_64FC1);
			temp_model.train(dud);
			final_models.push_back(temp_model);
			cv::Mat temp_mean = cv::Mat::zeros(g,1,CV_64FC1);
			cv::Mat temp_weight = cv::Mat::zeros(1,g,CV_64FC1);
			average_means.push_back(temp_mean);
			average_weights.push_back(temp_weight);

			std::vector<cv::Mat> temp_covs;
			for (int i =0; i<g; i++)
			{
				cv::Mat temp_cov = cv::Mat::zeros(1,1,CV_64FC1);
				temp_covs.push_back(temp_cov);
			}
			average_covs.push_back(temp_covs);
		}

		// Cross validate
		float alpha = 0.1;
		for(int h=0; h<10; h++)
		{
			// Remove 10% of the data

			cv::Mat sample_input (0,cols,CV_32FC1);
			int rows = cluster_input.rows;
			int number_selected = rows - rows/10;
			int  row = h*rows/10;
			std::vector<cv::EM> reduced_models;
			for(int k =0; k<number_selected;k++)
			{
				if(row + k < rows && row+k >=0)
				{
					sample_input.push_back(cluster_input.row(row+k));
				}
				if(row + k >= rows &&(row+k-rows) >=0 )
				{
					sample_input.push_back(cluster_input.row(row+k-rows));
				}
			}
			//Train the models

			for(int i=0; i<cols; i++)
			{
				double logs =0;
				cv::EM model (g,1,parameters);
				model.train(sample_input.col(i),log_likelihoods);
				for(int j = 0; j < number_selected; j++ )
				{					
					if(log_likelihoods.at<double>(j) + std::log(sample_input.at<float>(j,i))>-100000)
					{
					logs += log_likelihoods.at<double>(j) + log(sample_input.at<float>(j,i));
					}
				}
				logs = logs - (beta*g*fabs(logs));
				mean_of_logs[i] += 0.1*logs;
				reduced_models.push_back(model);
			}
			//Sum the models

			for(int i=0; i<cols; i++)
			{

				cv::Mat temp_mean = reduced_models[i].get<cv::Mat>("means");
				cv::Mat temp_weight = reduced_models[i].get<cv::Mat>("weights");			
				std::vector<cv::Mat> temp_covs = reduced_models[i].get<std::vector<cv::Mat>>("covs");

				for (int j = 0; j<g; j++)
				{
					float greatest = 0;
					int index =0;

					for (int k = 0; k<g; k++)
					{
						float temp = temp_mean.at<double> (k,0);

						if(temp >greatest)
						{
							greatest = temp;
							index = k;
						}
					}
					cv::addWeighted(average_means[i].row(j),1,temp_mean.row(index),alpha,0,average_means[i].row(j));					
					cv::addWeighted(average_weights[i].col(j),1,temp_weight.col(index),alpha,0,average_weights[i].col(j));
					for(int j=0; j<g; j++)
					{
						cv::addWeighted(average_covs[i][j],1,temp_covs[index],alpha,0,average_covs[i][j]);
					}
					temp_mean.at<double> (index,0) =0;
				}
			}
		}
		//EM HACK ----- copy averaged models to final model
		for(int h=0; h<cols; h++)
		{
			cv::Mat temp_means = final_models[h].get<cv::Mat>("means");
			cv::Mat temp_weights = final_models[h].get<cv::Mat>("weights");
			std::vector<cv::Mat> temp_covs = final_models[h].get<std::vector<cv::Mat>>("covs");

			average_means[h].copyTo(temp_means);
			for(int j =0; j<g; j++)
			{
				cv::Mat temp_mat;
				temp_covs.push_back(temp_mat);
				average_covs[h][j].copyTo(temp_covs[j]);
			}
			average_weights[h].copyTo(temp_weights);
		}
		sum_of_logs.push_back(mean_of_logs);
		all_final_models.push_back(final_models);
	}
	std::vector<cv::EM> chosen_models;
	for(int g=0; g<cols; g++)
	{	
		int index = 0;
		for(int h = 0; h < number_of_clusters; h++)
		{
			std::cout<< sum_of_logs[h][g] << ' ';
			if(sum_of_logs[h][g]>sum_of_logs[index][g])
			{
				index = h;
			}
		}
		chosen_models.push_back(all_final_models[index][g]);
		Clusters_for_feature[g] = index;
		std::cout << index << '\n';
	}
	return chosen_models;
}
std::vector<Ranked_cell> cellProfile (std::vector<Ranked_cell> sample, std::vector<cv::EM> models)
{
	int cols = sample[0].feature_values.size();
	int rows = sample.size();

	for(int i=0; i<rows; i++)
	{
		for(int j=0; j<cols; j++)
		{
			cv::Mat feature( 1, 1, CV_32FC1);
			feature.at<float>(0,0) = sample[i].feature_values[j];
			cv::Vec2d response = models[j].predict(feature);
			int cluster = (int)response[1];
			sample[i].cluster_index.push_back(cluster);
		}
	}
	return sample;
}
#endif