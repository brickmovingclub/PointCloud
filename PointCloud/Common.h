#pragma once
class Common
{
public:
	Common();
	~Common();
	/**************ͨ�ù��ܺ���****************/

	static void PrintString(const string &str);
//<<<<<<< HEAD
	//	�����еĵ���������Ƭͬ��
	static bool OnTheSameSide(const CVector &normal,const Point &origin,const std::list<Point> &nearPoints);

	//	�������Ĵ�С�Ƚ�
	#ifndef ABS
	#define ABS(x) ((x)<0?-(x):(x))//���xС��0������-x�����򷵻�x
	#endif
	template<class A>
	static bool fuzzyCompare1D(A a, A b)//�Ƚ�a�Ƿ�С��b
	{
		return ABS(a - b) < std::numeric_limits<A>::epsilon();
	}
	
	//	��ռ�ƽ��ķ�����
	static bool CalNormalVector(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3,float &dx, float &dy, float &dz);

	//	���ݿ���ķ�����������������ķ���
	static bool GetTwoLineIntersection(float _a1, float _b1, float _c1, float _a2, float _b2, float _c2, float &x, float &y);

//=======
	// �ж������Ƿ���ͬһֱ�����Լ�������Ƿ������㹹�ɵ�Բ����
	bool Condition_a_b(pcl::PointXYZ pi, pcl::PointXYZ pj, pcl::PointXYZ pk, std::vector<pcl::PointXYZ> near_pi);
//>>>>>>> origin/dev_hhy
	void findCandidatePoints(Point pi, Point pj, Point pk);
};

