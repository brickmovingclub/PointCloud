#pragma once

//	中頭
class CFace
{
public:
	CFace(const pcl::PointXYZ &point1, const pcl::PointXYZ &point2, const pcl::PointXYZ &point3);
	~CFace();

	Point GetPoint1()const;
	Point GetPoint2()const;
	Point GetPoint3()const;

	void  CFace::operator =(const CFace &p);
private:
	pcl::PointXYZ _point1;
	pcl::PointXYZ _point2;
	pcl::PointXYZ _point3;
};

// 利鯉爆中 std::list<CFace> ST