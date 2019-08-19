#include "stdafx.h"
#include "CLine.h"


CLine::CLine(const pcl::PointXYZ &pointStart, const pcl::PointXYZ &pointEnd)
{
	_pointStart = pointStart;
	_pointEnd = pointEnd;
}

CLine::CLine(const Point &pointStart, const Point &pointEnd)
{
	_pointStart.x = pointStart._x;	_pointStart.y = pointStart._y;	_pointStart.z = pointStart._z;
	_pointEnd.x = pointStart._x;	_pointEnd.y = pointStart._y;	_pointEnd.z = pointStart._z;
}
CLine::CLine()
{

}

CLine::CLine(const CLine &obj)
{
	this->_pointStart = obj._pointStart;
	this->_pointEnd = obj._pointEnd;
}

CLine::~CLine()
{
}


bool CLine::operator ==(const CLine& a)const
{
	return ((this->_pointStart.x == a._pointStart.x)
		&& (this->_pointStart.y == a._pointStart.y)
		&& (this->_pointStart.z == a._pointStart.z)
		&& (this->_pointEnd.x == a._pointEnd.x)
		&& (this->_pointEnd.y == a._pointEnd.y)
		&& (this->_pointEnd.z == a._pointEnd.z));
}

float CLine::LineLength(CLine line)
{
	float p1_x = line._pointStart.x;
	float p1_y = line._pointStart.y;
	float p1_z = line._pointStart.z;
	float p2_x = line._pointEnd.x;
	float p2_y = line._pointEnd.y;
	float p2_z = line._pointEnd.z;
	return sqrt(pow(p1_x - p2_x, 2) + pow(p1_y - p2_y, 2) + pow(p1_z - p2_z, 2));
}


float CLine::LineLength_Point(Point p1, Point p2)
{
	// TODO: 在此处添加实现代码.
	return sqrt(pow(p1._x - p2._x, 2) + pow(p1._y - p2._y, 2) + pow(p1._z - p2._z, 2));
}


Point CLine::getPointStart()
{
	// TODO: 在此处添加实现代码.
	return Point(_pointStart.x, _pointStart.y, _pointStart.z);
}


Point CLine::getPointEnd()
{
	// TODO: 在此处添加实现代码.
	return Point(_pointEnd.x, _pointEnd.y, _pointEnd.z);
}

pcl::PointXYZ CLine::GetPCLPointStart()
{
	return _pointStart;
}
pcl::PointXYZ CLine::GetPCLPointEnd()
{
	return _pointEnd;
}