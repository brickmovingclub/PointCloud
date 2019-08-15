#pragma once


//	œﬂ¿‡
class CLine
{
public:
	CLine();
	CLine(const CLine &obj);
	CLine(const pcl::PointXYZ &pointStart, const pcl::PointXYZ &pointEnd);
	~CLine();
private:
	pcl::PointXYZ _pointStart;
	pcl::PointXYZ _pointEnd;
};

