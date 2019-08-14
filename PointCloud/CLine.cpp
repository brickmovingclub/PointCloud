#include "stdafx.h"
#include "CLine.h"


CLine::CLine(struct Point pointStart, struct Point pointEnd)
{
	_pointStart = pointStart;
	_pointEnd = pointEnd;

}

CLine::CLine()
{
}

CLine::~CLine()
{
}


float CLine::LineLength(CLine line)
{
	float p1_x = line._pointStart._x;
	float p1_y = line._pointStart._y;
	float p1_z = line._pointStart._z;
	float p2_x = line._pointEnd._x;
	float p2_y = line._pointEnd._y;
	float p2_z = line._pointEnd._z;
	return sqrt(pow(p1_x - p2_x, 2) + pow(p1_y - p2_y, 2) + pow(p1_z - p2_z, 2));
}


float CLine::LineLength_Point(Point p1, Point p2)
{
	// TODO: 在此处添加实现代码.
	return sqrt(pow(p1._x - p2._x, 2) + pow(p1._y - p2._y, 2) + pow(p1._z - p2._z, 2));
}
