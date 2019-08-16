#pragma once


//	œﬂ¿‡
class CLine
{
public:
	CLine();
	CLine(const CLine &obj);
	CLine(const pcl::PointXYZ &pointStart, const pcl::PointXYZ &pointEnd);
	~CLine();
	bool CLine::operator ==(const CLine& a)const;

	
	pcl::PointXYZ _pointStart;
	pcl::PointXYZ _pointEnd;

	float LineLength(CLine line);
	float LineLength_Point(Point p1, Point p2);
	Point getPointStart();
	Point getPointEnd();
private:
};

