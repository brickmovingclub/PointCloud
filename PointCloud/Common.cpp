#include "stdafx.h"

#include "CVector.h"

#include "Common.h"


Common::Common()
{
}


Common::~Common()
{
}


void Common::PrintString(const string &str)
{

}


//<<<<<<< HEAD
bool Common::NoInclude()
{

	return false;
}

bool Common::Collineation(const Point &point1, const Point &point2, const Point &point3)
{
	return false;
}

bool Common::OnTheSameSide(const CVector &normal, const Point &origin, const std::list<Point> &nearPoints)
{
	int i = 0;
	double temp = 0.0f;
	vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
	plane->SetOrigin(origin.x, origin.y, origin.z);
	plane->SetNormal(normal.GetX(), normal.GetY(), normal.GetZ());

	for (auto iter : nearPoints)
	{
		temp = plane->EvaluateFunction(iter.x, iter.y, iter.z);
		if ((temp  - 0.0f)> 0)
			i++;
		else if ((temp - 0.0f) < 0)
			i--;
	}
	return (i == nearPoints.size() ? true :false);
}

bool Common::GetTwoLineIntersection(float _a1, float _b1, float _c1, float _a2, float _b2, float _c2, float &x, float &y)
{
	//_a1x+_b1y=_c1;---(1)
	//_a2x+_b2y=_c2;---(2)
	//
	if (_c1 == 0 && _c2 == 0)
	{
		//2个未知数2个方程组成的齐次方程组求解
		//系数矩阵B
		//| _a1 _b1 |
		//| _a2 _b2 |
		float DB = _a1 * _b2 - _a2 * _b1;
		if (DB != 0)//有唯一零解
		{
			x = 0;
			y = 0;
			return true;
		}
		else//有无数解
		{
			x = 0;
			y = 0;
			return false;
		}
	}
	else
	{
		//2个未知数2个方程组成的非齐次方程组求解
		//系数矩阵B
		//| _a1 _b1 |
		//| _a2 _b2 |
		//
		float DB = _a1 * _b2 - _a2 * _b1;
		if (DB != 0)//有唯一解
		{
			float dD1 = _c1 * _b2 - _c2 * _b1;
			float dD2 = _a1 * _c2 - _a2 * _c1;
			x = dD1 / DB;
			y = dD2 / DB;
			return true;
		}
		else//有无数解或者无解
		{
			x = 0;
			y = 0;
			return false;
		}
	}
	return false;
}

bool Common::CalNormalVector(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3,float &dx, float &dy, float &dz)
{
	//float a1,a2,b1,b2,c1,c2;
	//
	float _a1 = x2 - x1; float _b1 = y2 - y1; float _c1 = z2 - z1;
	float _a2 = x3 - x1; float _b2 = y3 - y1; float _c2 = z3 - z1;
	float _a3 = x3 - x2; float _b3 = y3 - y2; float _c3 = z3 - z2;
	//_a1x+_b1y+_c1z=0;
	//_a2x+_b2y+_c2z=0;
	//_a3x+_b3y+_c3z=0;
	//3个未知数3个方程组成的齐次方程组求解
	//系数矩阵A
	//| _a1 _b1 _c1 |
	//| _a2 _b2 _c2 |
	//| _a3 _b3 _c3 |
	//如果行列式A的值不等于0，则有唯一解且为零解
	float DA = _a1 * _b2*_c3 + _b1 * _c2*_a3 + _a2 * _b3*_c1 - _a3 * _b2*_c1 - _a1 * _b3*_c2 - _a2 * _b1*_c3;
	if (DA != 0)
	{
		dx = 0.0f;
		dy = 0.0f;
		dz = 0.0f;
		return false;
	}
	//---------------------------------------------//
	//如果行列式A的值等于0，则有非零解
	//非零解即x!=0时有解或者y!=0时有解或者z!=0时有解
	float x = 0.0f, y = 0.0f, z = 0.0f;
	//若z!=0时有解,取z=-1
	//_a1x+_b1y=_c1;---(1)
	//_a2x+_b2y=_c2;---(2)
	//_a3x+_b3y=_c3;---(3)
	//任取2个方程即可，在此取(1)(2)
	x = 0.0f; y = 0.0f;
	bool flag3 = GetTwoLineIntersection(_a1, _b1, _c1, _a2, _b2, _c2, x, y);
	if (flag3)//假设成立
	{
		dx = -x;
		dy = -y;
		dz = 1.0f;
		return true;
	}
	//假设不成立，继续试验另一个假设
	//若x!=0时有解取x=-1，平面中两条直线求交点问题
	//_b1y+_c1z=_a1;---(1)
	//_b2y+_c2z=_a2;---(2)
	//_b3y+_c3z=_a3;---(3)
	//任取2个方程即可，在此取(1)(2)
	y = 0.0f; z = 0.0f;
	bool flag1 = GetTwoLineIntersection(_b1, _c1, _a1, _b2, _c2, _a2, y, z);
	if (flag1)//假设成立
	{
		dx = 1.0f;
		dy = -y;
		dz = -z;
		return true;
	}
	//假设不成立，继续试验另一个假设
	//若y!=0时有解取y=-1，平面中两条直线求交点问题
	//_a1x+_c1z=_b1;---(1)
	//_a2x+_c2z=_b2;---(2)
	//_a3x+_c3z=_b3;---(3)
	//任取2个方程即可，在此取(1)(2)
	x = 0.0f; z = 0.0f;
	bool flag2 = GetTwoLineIntersection(_a1, _c1, _b1, _a2, _c2, _b2, x, z);
	if (flag2)//假设成立
	{
		dx = -x;
		dy = 1.0f;
		dz = -z;
		return true;
	}

	//所有假设都不成立，求解失败
	return false;
}


//=======
// 判断三点是否在同一直线上以及领域点是否在三点构成的圆球中
bool Common::Condition_a_b(pcl::PointXYZ pi, pcl::PointXYZ pj, pcl::PointXYZ pk, std::vector<pcl::PointXYZ> near_pi)
{
	// TODO: 在此处添加实现代码.	
	double a1, b1, c1, d1;
	double a2, b2, c2, d2;
	double a3, b3, c3, d3;

	double x1 = pi.x, y1 = pi.y, z1 = pi.z;
	double x2 = pj.x, y2 = pj.y, z2 = pj.z;
	double x3 = pk.x, y3 = pk.y, z3 = pk.z;

	//判断三点是否共线
	double s = 0.5 * (x1 * y2 + x2 * y3 + x3 * y1 - x1 * y3 - x2 * y1 - x3 * y2);
	if (s == 0)
		return false;

	//求圆心
	a1 = (y1 * z2 - y2 * z1 - y1 * z3 + y3 * z1 + y2 * z3 - y3 * z2);
	b1 = -(x1 * z2 - x2 * z1 - x1 * z3 + x3 * z1 + x2 * z3 - x3 * z2);
	c1 = (x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2);
	d1 = -(x1 * y2 * z3 - x1 * y3 * z2 - x2 * y1 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - x3 * y2 * z1);

	a2 = 2 * (x2 - x1);
	b2 = 2 * (y2 - y1);
	c2 = 2 * (z2 - z1);
	d2 = x1 * x1 + y1 * y1 + z1 * z1 - x2 * x2 - y2 * y2 - z2 * z2;

	a3 = 2 * (x3 - x1);
	b3 = 2 * (y3 - y1);
	c3 = 2 * (z3 - z1);
	d3 = x1 * x1 + y1 * y1 + z1 * z1 - x3 * x3 - y3 * y3 - z3 * z3;

	//圆心坐标
	pcl::PointXYZ centerpoint;
	centerpoint.x = -(b1*c2*d3 - b1 * c3*d2 - b2 * c1*d3 + b2 * c3*d1 + b3 * c1*d2 - b3 * c2*d1)
		/ (a1*b2*c3 - a1 * b3*c2 - a2 * b1*c3 + a2 * b3*c1 + a3 * b1*c2 - a3 * b2*c1);
	centerpoint.y = (a1*c2*d3 - a1 * c3*d2 - a2 * c1*d3 + a2 * c3*d1 + a3 * c1*d2 - a3 * c2*d1)
		/ (a1*b2*c3 - a1 * b3*c2 - a2 * b1*c3 + a2 * b3*c1 + a3 * b1*c2 - a3 * b2*c1);
	centerpoint.z = -(a1*b2*d3 - a1 * b3*d2 - a2 * b1*d3 + a2 * b3*d1 + a3 * b1*d2 - a3 * b2*d1)
		/ (a1*b2*c3 - a1 * b3*c2 - a2 * b1*c3 + a2 * b3*c1 + a3 * b1*c2 - a3 * b2*c1);

	//半径
	double dR = sqrt(pow(centerpoint.x - pi.x, 2) + pow(centerpoint.y - pi.y, 2) + pow(centerpoint.z - pi.z, 2));

	//判断领域点是否在三点构成的圆球中
	double distance;
	for (auto it = near_pi.begin(); it != near_pi.end(); it++)
	{
		distance = sqrt(pow(centerpoint.x - it->x, 2) + pow(centerpoint.y - it->y, 2) + pow(centerpoint.z - it->z, 2));
		if (distance < dR)
			return false;
	}
	return true;
//>>>>>>> origin/dev_hhy
}
