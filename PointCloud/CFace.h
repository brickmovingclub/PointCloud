#pragma once

//	��Ƭ
class CFace
{
public:
	CFace(const pcl::PointXYZ &point1, const pcl::PointXYZ &point2, const pcl::PointXYZ &point3);
	~CFace();
private:
	pcl::PointXYZ _point1;
	pcl::PointXYZ _point2;
	pcl::PointXYZ _point3;
};

// �������� std::list<CFace> ST