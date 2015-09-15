#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
using namespace cv;

#include <iostream>
using namespace std;

struct cloud_pt
{
	Point3f pt;
	double luminace;
};

struct sum_map
{
	int index;
	vector<cloud_pt> pts;
};

class map_manager
{
private:
	vector<cloud_pt> pts;

public:
	map_manager();
	void update(vector<cloud_pt>  pts);
	void clear();
	~map_manager();
};
