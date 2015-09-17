#include "map_manager.h"

map_manager::map_manager(int buf_sz)
{
	ofs.open("map.txt");
	ifs.open("map.txt");
	this->buf_sz = buf_sz;
	cpts = new cloud_pt[buf_sz];

	cpt_flags = new bool[buf_sz];
	for(int i = 0; i < buf_sz; ++i)
	{
		cpt_flags[i] = true;
	}

	max_sdist = 10000;
}

void map_manager::update(cloud_pt * new_cpts, int num_cpts)
{
	for(int i = 0; i < num_cpts; ++i)
	{
		for(int j = 0; j < buf_sz; ++j)
		{
			if(cpts[j].x == new_cpts[i].x 
			   && cpts[j].y == new_cpts[i].y 
			   && cpts[j].z == new_cpts[i].z)
				cpts[j] = new_cpts[i];
		}
	}
}

void map_manager::check(Point3i pos)
{
	int sdist;
	for(int i = 0; i < buf_sz; ++i)
	{
		sdist = 0;
		sdist += pow(cpts[i].x - pos.x, 2);
		sdist += pow(cpts[i].y - pos.y, 2);
		sdist += pow(cpts[i].z - pos.z, 2);
		
		if(sdist > max_sdist)
		{
			cpt_flags[i] = false;
		}
	}
}

void map_manager::save()
{
	char buf[46];
	for(int i = 0; i < buf_sz; ++i)
	{
		while(!ifs.eof())
		{
			ifs >> buf;
		
			Point3i pt;
			split_xyz(buf, 46, pt);
			if(cpts[i].x == pt.x &&
			   cpts[i].y == pt.y &&
			   cpts[i].z == pt.z)
			{
			
				break;
			}
			
			if(ifs.eof())
			{
				
			}
		}
	}
}

void map_manager::load(Point3i pos)
{
	check(pos);

	char buf[46];
	for(int i = 0; i < buf_sz; ++i)
	{
		while(cpt_flags[i])
		{
			ifs >> buf;
			
			Point3i pt;
			split_xyz(buf, 46, pt);
	
			if(cpts[i].x == pt.x  
			   && cpts[i].y == pt.y && cpts[i].z == pt.z)
			{
				continue;
			}
;
			int sdist = 0;
			sdist += pow(pt.x - pos.x, 2);
			sdist += pow(pt.y - pos.y, 2);
			sdist += pow(pt.z - pos.z, 2);

			if(sdist < max_sdist)
			{
				cpts[i].x = pt.x;
				cpts[i].y = pt.y;
				cpts[i].z = pt.z;
			}

			if(ifs.eof())
				break;

		}
	}

}

void map_manager::clear()
{
}

map_manager::~map_manager()
{
	clear();
}

void split_xyz(char * src, int buf_sz, Point3i &pt)
{
	char * buf = new char[buf_sz/4];
	int num_found = 0, start = 0;
	for(int i = 0; i < buf_sz; ++i)
	{
		if(src[i] == ' ')
		{
			memcpy(buf, src + start, i-start);
			start = i + 1;
			if(num_found == 0)
			{
				pt.x = atoi(buf);
			}
			else if(num_found == 1)
			{
				pt.y = atoi(buf);
			}
			else
			{
				pt.z = atoi(buf);
			}
		}
	}
}
