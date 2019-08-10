#pragma once
class Common
{
public:
	Common();
	~Common();
	
	/**************通用功能函数****************/
	static void PrintString(const string &str);
	// 判断三点是否在同一直线上以及领域点是否在三点构成的圆球中
	bool Condition_a_b(pcl::PointXYZ pi, pcl::PointXYZ pj, pcl::PointXYZ pk, std::vector<pcl::PointXYZ> near_pi);
};

