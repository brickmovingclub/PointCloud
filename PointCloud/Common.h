#pragma once
class Common
{
public:
	Common();
	~Common();
	
	/**************通用功能函数****************/
	static void PrintString(const string &str);
	//	三点共线
	static bool Collineation(const Point &point1, const Point &point2, const Point &point3);
	//	经过pi,pj,pk三点的圆球内不包括Near(Pi)中的领域点
	static bool NoInclude();
	//	领域中的点在三角面片同侧
	static bool OnTheSameSide(const CVector &normal,const Point &origin,const std::list<Point> &nearPoints);

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

	//	根据克莱姆法则求解带两个参数的方程
	static bool GetTwoLineIntersection(float _a1, float _b1, float _c1, float _a2, float _b2, float _c2, float &x, float &y);

};

