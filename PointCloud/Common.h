#pragma once
class Common
{
public:
	Common();
	~Common();
	/**************通用功能函数****************/

//<<<<<<< HEAD
	//	领域中的点在三角面片同侧(normal 为法向量，origin 为三角面片上任意一点)
	static bool OnTheSameSide(const CVector &normal,const pcl::PointXYZ &origin,const std::vector<std::pair<double, pcl::PointXYZ>> &nearPoints);
	//static bool OnTheSameSide(const pcl::PointXYZ &pi, const pcl::PointXYZ &pj, const pcl::PointXYZ &pk, const std::vector<std::pair<double,Point>> &nearPoints);

	//	浮点数的大小比较
	#ifndef ABS
	#define ABS(x) ((x)<0?-(x):(x))//如果x小于0，返回-x，否则返回x
	#endif
	template<class A>
	static bool fuzzyCompare1D(A a, A b)//比较a是否小于b
	{
		return ABS(a - b) < std::numeric_limits<A>::epsilon();
	}
	
	//	求空间平面的法向量
	static bool CalNormalVector(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3,float &dx, float &dy, float &dz);
	static void CalNormalVector(const Point &p1, const Point &o2, const Point &p3, CVector &vector);
	//	根据克莱姆法则求解带两个参数的方程
	static bool GetTwoLineIntersection(float _a1, float _b1, float _c1, float _a2, float _b2, float _c2, float &x, float &y);

//=======
	// 判断三点是否在同一直线上以及领域点是否在三点构成的圆球中
	static bool Condition_a_b(pcl::PointXYZ pi, pcl::PointXYZ pj, pcl::PointXYZ pk, std::vector<std::pair<double, pcl::PointXYZ>> &near_pi);

	//	半径领域点集搜索
	static void NearRadiusSearch(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,const pcl::PointXYZ &point, const float &radius, std::vector<std::pair< double, pcl::PointXYZ>> &nearPoint);

	//	k 近邻搜索
	static void KNearSearch(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, const pcl::PointXYZ &searchPoint, const int &k, std::vector<std::pair< double, pcl::PointXYZ>> &KnearPoint);

	//	体素搜索
	static void VoxelSearch(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, const pcl::PointXYZ &searchPoint, std::vector<std::pair< double, pcl::PointXYZ>> &KnearPoint);


	//	pcl 绘制线条
	static void PCLDrawLine(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, pcl::visualization::PCLVisualizer::Ptr viewer, std::vector<CLine> &activeList);
//>>>>>>> origin/dev_hhy
//<<<<<<< HEAD

//=======
	//选择候选点集
	static void  findCandidatePoints(pcl::PointCloud<pcl::PointXYZ>::Ptr _cloud, Point pi, Point pj, Point pk, std::map<Point, bool> flag, std::vector<CLine> ActiveE, CLine CurrentE, std::vector<std::pair< double, pcl::PointXYZ>> &result);
	
	//选择最佳点
	static Point FindBestPoint(pcl::PointCloud<pcl::PointXYZ>::Ptr _cloud, Point pi, Point pj, Point pk, std::vector<std::pair< double, pcl::PointXYZ>> &result, CLine CurrentE, std::list<CFace> ST, std::vector<CLine> InnerE, std::map<Point, bool> flag);
	
	// 寻找点的邻接三角形,并判断是否共边
	static bool findNearFace_Point(CFace curFace, Point pc, std::list<CFace> ST);
	
	// 判断两个三角形是否相交
	static bool IntersectTriangle(Point pi, Point pj, Point pc, std::list<CFace> ST);
	
	// 判断两条线段是否相交
	static bool IntersectionLine(Point pi, Point pj, Point pc, Point pk);
	
	// 判断前面三个点构成的三角形是否包含第四个点
	static bool TriangleIncludeSubpoint(Point pi, Point pj, Point bestP, Point subpoint);
	
	// 得到三个点组成的三角形的面积
	static float TriangleArea(Point A, Point B, Point C);
	
	// 更新活动链表
	void UpdateActiveList(std::vector<CLine> ActiveE, CLine CurrentE, Point bestP, std::vector<CLine> InnerE, std::vector<Point> FreeP, std::vector<Point> ActiveP, std::map<Point, bool> flag, std::vector<CFace> ST);
	// 判断最佳点添加的位置类型
	int BestPositionType(std::vector<CLine> ActiveE, CLine CurrentE, Point bestP, std::vector<Point> FreeP);
	/*************************更新活动链表*************************************/
	// 最佳点是自由点
	void UpdateMode(std::vector<CLine> ActiveE, CLine CurrentE, Point bestP, std::vector<CLine> InnerE, std::vector<Point> FreeP, std::vector<Point> ActiveP);
	// 最佳点位于活动边上且为当前活动边前相邻边的端点
	void UpdateMode1(std::vector<CLine> ActiveE, CLine CurrentE, Point bestP, std::vector<CLine> InnerE, std::vector<Point> ActiveP, std::map<Point, bool> flag);
	// 最佳点位于活动边上且为当前活动边后相邻边的端点
	void UpdateMode2(std::vector<CLine> ActiveE, CLine CurrentE, Point bestP, std::vector<CLine> InnerE, std::vector<Point> ActiveP, std::map<Point, bool> flag);
	//最佳点位于活动边上且与当前活动边没有相邻关系
	void UpdateMode3(std::vector<CLine> ActiveE, CLine CurrentE, Point bestP, std::vector<CLine> InnerE, std::vector<Point> ActiveP, std::map<Point, bool> flag, std::vector<CFace> ST);
//>>>>>>> origin/dev_hhy
};

