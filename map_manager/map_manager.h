#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
using namespace cv;

#include <iostream>
#include <fstream>
using namespace std;

struct cloud_pt
{
	int x, y, z;
	double luminance;
};


class map_manager
{
private:
	int max_sdist;
	ofstream ofs;
	ifstream ifs;
	cloud_pt * cpts;
	bool * cpt_flags;
	int buf_sz;
	Point3i position;

public:
	map_manager(int buf_sz);
	void update(cloud_pt * new_pts, int num_cpts);
	void clear();
	void load(Point3i pos);
	void save();
	void check(Point3i pos);
	~map_manager();
};

void split_xyz(char * src, int buf_sz, Point3i &pt);
