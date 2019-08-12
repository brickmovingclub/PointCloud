#include "stdafx.h"
#include "CFace.h"


CFace::CFace(const struct Point &point1, const struct Point &point2, const struct Point &point3)
{
	_point1 = point1;
	_point2 = point2;
	_point3 = point3;
}


CFace::~CFace()
{
}
