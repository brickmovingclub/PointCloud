#include "stdafx.h"
#include "CFace.h"


CFace::CFace(const pcl::PointXYZ &point1, const pcl::PointXYZ &point2, const pcl::PointXYZ &point3)
{
	_point1._x = point1.x;	_point1._y = point1.y;	_point1._z = point1.z;

	_point2._x = point2.x;	_point2._y = point2.y;	_point2._z = point2.z;
	_point3._x = point3.x;	_point3._y = point3.y;	_point3._z = point3.z;

	
}
CFace::CFace(const Point &point1, const Point &point2, const Point &point3)
{
	_point1 = point1;
	_point2 = point2;
	_point3 = point3;
}


CFace::~CFace()
{
}


Point CFace::GetPoint1()const
{
	return _point1;
}
Point CFace::GetPoint2()const
{
	return _point2;
}
Point CFace::GetPoint3()const
{
	return _point3;
}

void  CFace::operator =(const CFace &p)
{
	this->_point1._x = p._point1._x;
	this->_point1._y = p._point1._y;
	this->_point1._z = p._point1._z;

	this->_point2._x = p._point2._x;
	this->_point2._y = p._point2._y;
	this->_point2._z = p._point2._z;

	this->_point3._x = p._point3._x;
	this->_point3._y = p._point3._y;
	this->_point3._z = p._point3._z;
}