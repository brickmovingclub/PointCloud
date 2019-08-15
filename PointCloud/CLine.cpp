#include "stdafx.h"
#include "CLine.h"


CLine::CLine(const pcl::PointXYZ &pointStart, const pcl::PointXYZ &pointEnd)
{
	_pointStart = pointStart;
	_pointEnd = pointEnd;
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
