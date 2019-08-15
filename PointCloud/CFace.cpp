#include "stdafx.h"
#include "CFace.h"


CFace::CFace(const pcl::PointXYZ &point1, const pcl::PointXYZ &point2, const pcl::PointXYZ &point3)
{
	_point1 = point1;
	_point2 = point2;
	_point3 = point3;
}


CFace::~CFace()
{
}
