#pragma once
class Common
{
public:
	Common();
	~Common();
	
	/**************ͨ�ù��ܺ���****************/
	static void PrintString(const string &str);
	// �ж������Ƿ���ͬһֱ�����Լ�������Ƿ������㹹�ɵ�Բ����
	bool Condition_a_b(pcl::PointXYZ pi, pcl::PointXYZ pj, pcl::PointXYZ pk, std::vector<pcl::PointXYZ> near_pi);
};

