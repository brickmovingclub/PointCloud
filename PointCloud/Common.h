#pragma once
class Common
{
public:
	Common();
	~Common();
	/**************ͨ�ù��ܺ���****************/

//<<<<<<< HEAD
	//	�����еĵ���������Ƭͬ��(normal Ϊ��������origin Ϊ������Ƭ������һ��)
	static bool OnTheSameSide(const CVector &normal,const pcl::PointXYZ &origin,const std::vector<std::pair<double, pcl::PointXYZ>> &nearPoints);
	//static bool OnTheSameSide(const pcl::PointXYZ &pi, const pcl::PointXYZ &pj, const pcl::PointXYZ &pk, const std::vector<std::pair<double,Point>> &nearPoints);

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
	static bool Condition_a_b(pcl::PointXYZ pi, pcl::PointXYZ pj, pcl::PointXYZ pk, std::vector<std::pair<double, pcl::PointXYZ>> &near_pi);

	//	�뾶����㼯����
	static void NearRadiusSearch(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,const pcl::PointXYZ &point, const float &radius, std::vector<std::pair< double, pcl::PointXYZ>> &nearPoint);

	//	k ��������
	static void KNearSearch(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, const pcl::PointXYZ &searchPoint, const int &k, std::vector<std::pair< double, pcl::PointXYZ>> &KnearPoint);

	//	��������
	static void VoxelSearch(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, const pcl::PointXYZ &searchPoint, std::vector<std::pair< double, pcl::PointXYZ>> &KnearPoint);


	//	pcl ��������
	static void PCLDrawLine(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, pcl::visualization::PCLVisualizer::Ptr viewer, std::list<CLine> &activeList);
//>>>>>>> origin/dev_hhy
//<<<<<<< HEAD
	static std::vector<int> findCandidatePoints(pcl::PointCloud<pcl::PointXYZ>::Ptr _cloud, Point pi, Point pj, Point pk, std::vector<bool> flag,std::vector<CLine> ActiveE, CLine CurrentE);

	static void  GetNormal(pcl::PointXYZ &p1, pcl::PointXYZ &p2, pcl::PointXYZ &p3, CVector &vector);

//=======
	//ѡ���ѡ�㼯
	std::vector<int> findCandidatePoints(pcl::PointCloud<pcl::PointXYZ>::Ptr _cloud, Point pi, Point pj, Point pk, std::map<Point, bool> flag,std::vector<CLine> ActiveE, CLine CurrentE);
	
	//ѡ����ѵ�
	Point FindBestPoint(pcl::PointCloud<pcl::PointXYZ>::Ptr _cloud, Point pi, Point pj, Point pk, std::vector<int> near_pm, CLine CurrentE, std::list<CFace> ST, std::vector<CLine> InnerE, std::map<Point, bool> flag);
	
	// Ѱ�ҵ���ڽ�������,���ж��Ƿ񹲱�
	bool findNearFace_Point(CFace curFace, Point pc, std::list<CFace> ST);
	
	// �ж������������Ƿ��ཻ
	bool IntersectTriangle(Point pi, Point pj, Point pc, std::list<CFace> ST);
	
	// �ж������߶��Ƿ��ཻ
	bool IntersectionLine(Point pi, Point pj, Point pc, Point pk);
	
	// �ж�ǰ�������㹹�ɵ��������Ƿ�������ĸ���
	bool TriangleIncludeSubpoint(Point pi, Point pj, Point bestP, Point subpoint);
	
	// �õ���������ɵ������ε����
	float TriangleArea(Point A, Point B, Point C);
	
	// ���»����
	void UpdateActiveList(std::vector<CLine> ActiveE, CLine CurrentE, Point bestP, std::vector<CLine> InnerE, std::vector<Point> FreeP, std::vector<Point> ActiveP, std::map<Point, bool> flag, std::vector<CFace> ST);
	// �ж���ѵ���ӵ�λ������
	int BestPositionType(std::vector<CLine> ActiveE, CLine CurrentE, Point bestP, std::vector<Point> FreeP);
	/*************************���»����*************************************/
	// ��ѵ������ɵ�
	void UpdateMode(std::vector<CLine> ActiveE, CLine CurrentE, Point bestP, std::vector<CLine> InnerE, std::vector<Point> FreeP, std::vector<Point> ActiveP);
	// ��ѵ�λ�ڻ������Ϊ��ǰ���ǰ���ڱߵĶ˵�
	void UpdateMode1(std::vector<CLine> ActiveE, CLine CurrentE, Point bestP, std::vector<CLine> InnerE, std::vector<Point> ActiveP, std::map<Point, bool> flag);
	// ��ѵ�λ�ڻ������Ϊ��ǰ��ߺ����ڱߵĶ˵�
	void UpdateMode2(std::vector<CLine> ActiveE, CLine CurrentE, Point bestP, std::vector<CLine> InnerE, std::vector<Point> ActiveP, std::map<Point, bool> flag);
	//��ѵ�λ�ڻ�������뵱ǰ���û�����ڹ�ϵ
	void UpdateMode3(std::vector<CLine> ActiveE, CLine CurrentE, Point bestP, std::vector<CLine> InnerE, std::vector<Point> ActiveP, std::map<Point, bool> flag, std::vector<CFace> ST);
//>>>>>>> origin/dev_hhy
};

