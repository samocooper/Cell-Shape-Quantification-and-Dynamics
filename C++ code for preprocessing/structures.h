#ifndef structures_H
#define structures_H

#include <vector>
#include <opencv2\highgui\highgui.hpp>
#include <string>

extern const unsigned char Number_of_features = 16;
extern unsigned int Current_feature = 0;
extern unsigned int Clusters_for_feature[100] = {};

struct Cell_static
{	
	std::string well_name;
	int well_index;
	float Features[Number_of_features];
};
struct Cell
{
	short int X_pos;
	short int Y_pos;
	short int field;
	float moved;
	float Features[Number_of_features];
	float Features_dynamic[Number_of_features];
	bool tracked;
};
struct Tracked_cell
{
	int frames_tracked;
	std::vector<int> X_values;
	std::vector<int> Y_values;
	std::vector<float> moved;
	std::vector<std::vector<float>> feature_values;
	std::vector<std::vector<float>> dynamic_values;
	int condition;
};
struct Ranked_cell
{
	std::vector<float> feature_values;
	std::vector<float> dynamic_values;
	std::vector<int> feature_ranks;
	std::vector<float> cluster_index;
};
struct Well_name
{
	std::string well;
	std::string condition;
};
struct TC
{
	std::string condition;	
	std::string well;
	std::vector<Tracked_cell> cells;
};
struct Compared_row
{	
	int next;
	bool is_last;
	bool integrated;
	std::vector<uchar> value;
	int left;
	int right;
	std::string name;
	int index;
};


#endif