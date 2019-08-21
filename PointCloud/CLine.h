#pragma once


//	ÏßÀà
class CLine
{
public:
	CLine();
	CLine(const CLine &obj);
	CLine(const pcl::PointXYZ &pointStart, const pcl::PointXYZ &pointEnd);
	CLine(const Point &pointStart, const Point &pointEnd);
	~CLine();
	bool CLine::operator ==(const CLine& a)const;

	void operator= (const CLine &line)
	{
		this->_pointStart = line._pointStart;
		this->_pointEnd = line._pointEnd;
	}
	pcl::PointXYZ _pointStart;
	pcl::PointXYZ _pointEnd;

	float LineLength(CLine line);
	static float LineLength_Point(Point p1, Point p2);
	pcl::PointXYZ GetPCLPointStart();
	pcl::PointXYZ GetPCLPointEnd();
	Point getPointStart();
	Point getPointEnd();
private:
};

