#include "map_manager.h"

map_manager::map_manager()
{	
}

void update(vector<cloud_pt> new_pts)
{
	vector<cloud_pt>::iterator new_pts_it, pts_it;
	new_pts_it = new_pts.begin();
	pts_it = pts.end();
	
	for(;new_pts_it != new_pts.end(); ++new_pts_it, ++pts_it)
		*pts_it = *new_pts_it;
	

	
}

void map_manager::clear()
{
}

map_manager::~map_manager()
{
	clear();
}
