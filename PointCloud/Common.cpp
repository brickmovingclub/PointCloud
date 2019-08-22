#include "stdafx.h"

#include "CLine.h"

#include "CVector.h"

#include "CFace.h"

#include "Common.h"


Common::Common()
{
}


Common::~Common()
{
}



//	领域中的点在三角面片同侧
bool Common::OnTheSameSide(const CVector &normal, const pcl::PointXYZ &origin, const std::vector<std::pair<double, pcl::PointXYZ>> &nearPoints)
{
	int i = 0,j = 0;
	double temp = 0.0f;
	vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
	plane->SetOrigin(origin.x, origin.y, origin.z);
	plane->SetNormal(normal.GetX(), normal.GetY(), normal.GetZ());

	for (auto iter : nearPoints)
	{
		temp = plane->EvaluateFunction(iter.second.x, iter.second.y, iter.second.z);
		if ((temp - 0.0f) > 0)
			i++;
		else if ((temp - 0.0f) < 0)
			i--;
		else
		{
			j++;
			//std::cout << 123 << std::endl;
		}
	}
	return (((i +j ) == nearPoints.size())? true : false);
}
/*
bool Common::OnTheSameSide(const pcl::PointXYZ &pi, const pcl::PointXYZ &pj, const pcl::PointXYZ &pk, const std::vector<std::pair<double, Point>> &nearPoints)
{
	//	求空间平面的法向量
	float dx, dy, dz;
	//float a1,a2,b1,b2,c1,c2;
	//
	float _a1 = pj.x - pi.x; float _b1 = pj.y - pi.y; float _c1 = pj.z - pi.z;
	float _a2 = pk.x - pi.x; float _b2 = pk.y - pi.y; float _c2 = pk.z - pi.z;
	float _a3 = pk.x - pj.x; float _b3 = pk.y - pj.y; float _c3 = pk.z - pj.z;
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
*/
//=======

//	根据克莱姆法则求解带两个参数的方程
//>>>>>>> origin/dev_hhy
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

//<<<<<<< HEAD
//bool Common::CalNormalVector(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3, float &dx, float &dy, float &dz)
//=======

//	求空间平面的法向量
void Common::CalNormalVector(const Point &p1, const Point &p2, const Point &p3, CVector &vector)
{
	CVector a(p2._x - p1._x, p2._y - p1._y, p2._z - p1._z);
	CVector b(p3._x - p1._x, p3._y - p1._y, p3._z - p1._z);

	float na = a.GetY() * b.GetZ() - a.GetZ() * b.GetY();
	float nb = a.GetZ() * b.GetX() - a.GetX() * b.GetZ();
	float nc = a.GetX() * b.GetY() - a.GetY() * b.GetX();
	//float dx, dy, dz;
	//bool flag  = CalNormalVector(p1._x,p1._y,p1._z,p2._x,p2._y,p2._z,p3._x,p3._y,p3._z,dx,dy,dz);
	vector.SetVector(na, nb, nc);
}

float Common::Round(const float &src, const int &bits)
{
	float f = src;
	if (bits > 0)
	{
		stringstream ss;
		ss << fixed << setprecision(bits) << f;
		ss >> f;
	}
	return f;
}
double Common::CalRadius(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud)
{
	pcl::PointXYZ min;//用于存放三个轴的最小值
	pcl::PointXYZ max;//用于存放三个轴的最大值
	float dx, dy, dz;
	double temp = 1.0f;
	pcl::getMinMax3D(*cloud, min, max);
	dx = (max.x - min.x);
	dy = max.y - min.y;
	dz = max.z - min.z;
	if (dx != 0)
		temp *= dx;
	if (dy != 0.0f)
		temp *= dy;
	if (dz != 0.0f)
		temp *= dz;
	temp = (double)(temp / cloud->points.size());
	return (double)sqrt(temp);	//	搜索半径r
}

bool Common::CalNormalVector(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3,float &dx, float &dy, float &dz)
//>>>>>>> origin/dev_hhy
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
bool Common::Condition_a_b(pcl::PointXYZ pi, pcl::PointXYZ pj, pcl::PointXYZ pk, std::vector<std::pair<double,pcl::PointXYZ>> &near_pi)
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
		distance = sqrt(pow(centerpoint.x - it->second.x, 2) + pow(centerpoint.y - it->second.y, 2) + pow(centerpoint.z - it->second.z, 2));
		if (distance < dR)
			return false;		
	}
	return true;
	//>>>>>>> origin/dev_hhy
}

void Common::NearRadiusSearch(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, pcl::PointCloud<pcl::RGB>::Ptr &pPointsRGB,const pcl::PointXYZ &point, const float &radius, std::vector<std::pair< double,pcl::PointXYZ>> &nearPoint)
{
	//  // Neighbors within radius search

	//pcl::PointXYZ searchPoint;
	//searchPoint.x = point._x;
	//searchPoint.y = point._y;
	//searchPoint.z = point._z;

	float resolution = 5.0f;

	pcl::octree::OctreePointCloudSearch<pcl::PointXYZ> octree(resolution);

	octree.setInputCloud(cloud);
	octree.addPointsFromInputCloud();


	std::vector<int> pointIdxRadiusSearch;
	std::vector<float> pointRadiusSquaredDistance;

	if (octree.radiusSearch(point, radius, pointIdxRadiusSearch, pointRadiusSquaredDistance) > 0)
	{
		int size = pointIdxRadiusSearch.size();
		for (size_t i = 0; i < pointIdxRadiusSearch.size(); ++i)
		{
			nearPoint.push_back(std::make_pair(pointRadiusSquaredDistance[i], cloud->points[pointIdxRadiusSearch[i]]));
			pPointsRGB->points[pointIdxRadiusSearch[i]].b = 255;		//所有领域点标记为蓝色
			pPointsRGB->points[pointIdxRadiusSearch[i]].r = 0;
				

			//std::cout << "    " << cloud->points[pointIdxRadiusSearch[i]].x
			//<< " " << cloud->points[pointIdxRadiusSearch[i]].y
		//	<< " " << cloud->points[pointIdxRadiusSearch[i]].z
			//<< " (squared distance: " << pointRadiusSquaredDistance[i] << ")" << std::endl;
		}
			
	}


}

void Common::KNearSearch(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,const pcl::PointXYZ &searchPoint,const int &k, std::vector<std::pair< double, pcl::PointXYZ>> &KnearPoint)
{
	

	float resolution = 5.0f;

	pcl::octree::OctreePointCloudSearch<pcl::PointXYZ> octree(resolution);

	octree.setInputCloud(cloud);
	octree.addPointsFromInputCloud();

	std::vector<int> pointIdxNKNSearch;
	std::vector<float> pointNKNSquaredDistance;


	if (octree.nearestKSearch(searchPoint, k, pointIdxNKNSearch, pointNKNSquaredDistance) > 0)
	{
		int size = pointIdxNKNSearch.size();

		for (size_t i = 0; i < pointIdxNKNSearch.size(); ++i)
		{
			KnearPoint.push_back(std::make_pair(pointNKNSquaredDistance[i], cloud->points[pointIdxNKNSearch[i]]));

			std::cout << "    " << cloud->points[pointIdxNKNSearch[i]].x
				<< " " << cloud->points[pointIdxNKNSearch[i]].y
				<< " " << cloud->points[pointIdxNKNSearch[i]].z
				<< " (squared distance: " << pointNKNSquaredDistance[i] << ")" << std::endl;
		}
			
	}
}

void Common::VoxelSearch(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, const pcl::PointXYZ &searchPoint, std::vector<std::pair< double, pcl::PointXYZ>> &KnearPoint)
{

	float resolution = 5.0f;

	pcl::octree::OctreePointCloudSearch<pcl::PointXYZ> octree(resolution);

	octree.setInputCloud(cloud);
	octree.addPointsFromInputCloud();

	std::vector<int> pointIdxVec;

	if (octree.voxelSearch(searchPoint, pointIdxVec))
	{
		std::cout << "Neighbors within voxel search at (" << searchPoint.x
			<< " " << searchPoint.y
			<< " " << searchPoint.z << ")"
			<< std::endl;

		int size = pointIdxVec.size();

		for (size_t i = 0; i < pointIdxVec.size(); ++i)
		{
			KnearPoint.push_back(std::make_pair(0, cloud->points[pointIdxVec[i]]));
			std::cout << "    " << cloud->points[pointIdxVec[i]].x
				<< " " << cloud->points[pointIdxVec[i]].y
				<< " " << cloud->points[pointIdxVec[i]].z << std::endl;
		}
			
	}
}
void Common::PCLDrawLine(QString &faceName, pcl::visualization::PCLVisualizer::Ptr viewer, const CFace &face,const QColor &color )
{
	int i = 0;
	string temp;
	QString name = faceName + QString::number(i++);
	temp = name.toStdString();
	viewer->addLine<pcl::PointXYZ>(face.GetPCLPoint1(), face.GetPCLPoint2(), color.red(), color.green(), color.blue(), temp.c_str()); //红色线段,线的名字叫做"line1
	name = faceName + QString::number(i++);
	temp = name.toStdString();
	viewer->addLine<pcl::PointXYZ>(face.GetPCLPoint2(), face.GetPCLPoint3(), color.red(), color.green(), color.blue(), temp.c_str()); //红色线段,线的名字叫做"line1
	name = faceName + QString::number(i++);
	temp = name.toStdString();
	viewer->addLine<pcl::PointXYZ>(face.GetPCLPoint3(), face.GetPCLPoint1(), color.red(), color.green(), color.blue(), temp.c_str()); //红色线段,线的名字叫做"line1

}


void Common::PCLDrawLine(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, pcl::visualization::PCLVisualizer::Ptr viewer, std::list<CFace> &ST)
{
	//cloud = getpoint();//实时获取点云
	pcl::PointXYZ  minPt, maxPt;
	//viewer->removeAllShapes();
	pcl::getMinMax3D(*cloud, minPt, maxPt);

	//pcl::PointXYZ origin(0, 0, 0);

	//viewer->setBackgroundColor(125, 125, 125); //背景色
	//pcl::ModelCoefficients plane;
	int i = 0;
	QString name;
	string temp;
	for (auto iter : ST)
	{
		name = QString("line") + QString::number(i++);
		temp = name.toStdString();
		viewer->addLine<pcl::PointXYZ>(iter.GetPCLPoint1(), iter.GetPCLPoint2(), 255, 0, 0,temp.c_str()); //红色线段,线的名字叫做"line1
		name = QString("line") + QString::number(i++);
		temp = name.toStdString();
		viewer->addLine<pcl::PointXYZ>(iter.GetPCLPoint2(), iter.GetPCLPoint3(), 255, 0, 0, temp.c_str()); //红色线段,线的名字叫做"line1
		name = QString("line") + QString::number(i++);
		temp = name.toStdString();
		viewer->addLine<pcl::PointXYZ>(iter.GetPCLPoint3(), iter.GetPCLPoint1(), 255, 0, 0, temp.c_str()); //红色线段,线的名字叫做"line1


	}
	//viewer->addPointCloud<pcl::PointXYZ>(cloud, "cloud");

	viewer->spinOnce(100);
}

Point Common::GetOtherPoint(const Point &pointi, const Point &pointj, std::list<CFace> &ST)
{
	int i = 0;
	for (auto iter : ST)
	{
		i = 0;
		if (pointi == iter.GetPoint1() || pointi == iter.GetPoint2() || pointi == iter.GetPoint3())
			i++;
		else
			continue;
		if (pointj == iter.GetPoint1() || pointj == iter.GetPoint2() || pointj == iter.GetPoint3())
			i++;
		else
			continue;
		if (i == 2)
		{
			if (pointi != iter.GetPoint1() && pointj != iter.GetPoint1())
				return iter.GetPoint1();
			else if (pointi != iter.GetPoint2() && pointj != iter.GetPoint2())
				return iter.GetPoint2();
			else if (pointi != iter.GetPoint3() && pointj != iter.GetPoint3())
				return iter.GetPoint3();
			
		}
	}
	return Point();
}

// 去除冗余点
void Common::RemoveRedundantPoints(Point pi, Point pj, Point bestP, std::map<Point, bool> &flag, std::vector<std::pair< double, pcl::PointXYZ>> &result)
{
	// TODO: 在此处添加实现代码.
	//检测是否有冗余点
	CVector vec_new;
	CalNormalVector(pi, pj, bestP, vec_new);

	float a = vec_new.GetX();
	float b = vec_new.GetY();
	float c = vec_new.GetZ();
	float d = -a * pi._x - b * pi._y - c * pi._z;

	for (auto it = result.begin(); it != result.end(); it++)
	{
		if ((it->second.x == bestP._x)
			&& (it->second.y == bestP._y)
			&& (it->second.z == bestP._z))
			continue;
		Point  spatialpoint;;
		spatialpoint._x = it->second.x;
		spatialpoint._y = it->second.y;
		spatialpoint._z = it->second.z;
		Point subpoint; //投影点
		subpoint._x = ((b * b + c * c) * spatialpoint._x - a * (b * spatialpoint._y + c * spatialpoint._z + d)) / (a * a + b * b + c * c);
		subpoint._y = b / a * (subpoint._x - spatialpoint._x) + spatialpoint._y;
		subpoint._z = c / a * (subpoint._x - spatialpoint._x) + spatialpoint._z;

		//判断投影点是否在三角形内
		if (TriangleIncludeSubpoint(pi, pj, bestP, subpoint))
			flag[spatialpoint] = true; //当前点设置为自由点/排除点
	}
}
//选择候选点集
void  Common::findCandidatePoints(pcl::PointCloud<pcl::PointXYZ>::Ptr _cloud, pcl::PointCloud<pcl::RGB>::Ptr &pPointsRGB, Point pi, Point pj, Point pk, std::map<Point, bool> flag, std::vector<CLine> ActiveE, CLine CurrentE, std::vector<std::pair< double, pcl::PointXYZ>> &result)
{
	//得到三边边长，判断三角形类型
	CLine l;
	float line_ij = l.LineLength_Point(pi, pj);
	float line_ik = l.LineLength_Point(pi, pk);
	float line_jk = l.LineLength_Point(pj, pk);
	float a, b, c;
	c = max(line_ij, line_ik);
	c = max(c, line_jk);
	if (line_ij == c)
	{
		a = line_ik;
		b = line_jk;
	}
	else if (line_ik == c)
	{
		a = line_ij;
		b = line_jk;
	}
	else
	{
		a = line_ij;
		b = line_ik;
	}
	//求领域半径r
	float r;
	if (c * c > a * a + b * b) //钝角三角形
		r = c;
	else
		r = (a + b + c) / 3;

	pcl::PointXYZ pm;
	pm.x = (pi._x + pj._x) / 2;
	pm.y = (pi._y + pj._y) / 2;
	pm.z = (pi._z + pj._z) / 2;

	

	//计算点pm的r范围内的领域点集
	//std::vector<int> near_pm; 
	
	Common::NearRadiusSearch(_cloud, pPointsRGB, pm, 2 * r , result);

	//精简领域点集
	//int size = near_pm.size();
	//std::cout << "-**************************************" << std::endl;
	for (auto it = result.begin(); it != result.end();)
	{
		//std::cout << it->second.x << "  " << it->second.y << "  " << it->second.z << std::endl;
		//删除固定点与排除点,以及当前活动边上的两个点
		Point point(it->second.x, it->second.y, it->second.z);
		if (flag[point] || point == pi || point == pj || point == pk)
			it = result.erase(it);
		else
			++it;
	}
	int resultSize = result.size();

	//找到pi,pj相邻的活动点(可能存在问题，当活动边并没有顺序时，此可算法失效)
	int i = 0;
	for (; i < ActiveE.size();)
	{
		if (ActiveE[i] == CurrentE)
			break;
		else
			i++;
	}
	Point pa = ActiveE[(i - 1 + ActiveE.size()) % ActiveE.size()].getPointStart();
	Point pb = ActiveE[(i + 1 + ActiveE.size()) % ActiveE.size()].getPointEnd();
	CVector vec; //pi->pj->pk的法向量
	Common::CalNormalVector(pi, pj, pk, vec);

	//边角度约束简化
	float angle = 0.0f;

	CVector p_ij(pj._x - pi._x, pj._y - pi._y, pj._z - pi._z);
	CVector p_ia(pa._x - pi._x, pa._y - pi._y, pa._z - pi._z);
	float angle_ij_ia_cos = CVector::vectorInnerProduct(p_ij, p_ia) / (p_ij.vectorMag(p_ij) * p_ia.vectorMag(p_ia));

	CVector f = CVector::Multiplication(p_ij, p_ia); //ij与ia所在平面的法向量
	float f_vec = (CVector::vectorInnerProduct(f, vec) / (f.vectorMag(f) * vec.vectorMag(vec)));

	if (f_vec >= 0)
		angle = 360 - acos(angle_ij_ia_cos);
	else
		angle = acos(angle_ij_ia_cos);


	float angle1 = 0.0f;

	CVector p_ji(pi._x - pj._x, pi._y - pj._y, pi._z - pj._z);
	CVector p_jb(pb._x - pj._x, pb._y - pj._y, pb._z - pj._z);
	float angle_ji_jb_cos = CVector::vectorInnerProduct(p_ji, p_jb) / (p_ji.vectorMag(p_ji) * p_jb.vectorMag(p_jb));

	CVector f1 = CVector::Multiplication(p_ji, p_jb); //ib与ji所在平面的法向量
	float f1_vec = CVector::vectorInnerProduct(f1, vec) / (f1.vectorMag(f1) * vec.vectorMag(vec));

	if (f1_vec <= 0)
		angle1 = 360 - acos(angle_ji_jb_cos);
	else
		angle1 = acos(angle_ji_jb_cos);

	float angleA = min(angle, 135.0f);
	float angleB = min(angle1, 135.0f);

	for (auto it = result.begin(); it != result.end();)
	{
		//int resultSize = result.size(); 
		Point candidateP;
		candidateP._x = it->second.x;
		candidateP._y = it->second.y;
		candidateP._z = it->second.z;

		CVector p_icandidateP(candidateP._x - pi._x, candidateP._y - pi._y, candidateP._z - pi._z);
		float angle_ij_icandidateP_cos = CVector::vectorInnerProduct(p_ij, p_icandidateP) / (p_ij.vectorMag(p_ij) * p_icandidateP.vectorMag(p_icandidateP));
		angle = acos(angle_ij_icandidateP_cos);

		CVector p_jcandidateP(candidateP._x - pj._x, candidateP._y - pj._y, candidateP._z - pj._z);
		float angle_ji_jcandidateP_cos = CVector::vectorInnerProduct(p_ji, p_jcandidateP) / (p_ji.vectorMag(p_ji) * p_jcandidateP.vectorMag(p_jcandidateP));
		angle1 = acos(angle_ji_jcandidateP_cos);

		if (angle >= angleA || angle1 >= angleB)
			it = result.erase(it);
		else
			it++;
	}

	//面角度约束简化 
	for (auto iter = result.begin(); iter != result.end();)
	{
		Point candidateP;
		candidateP._x = iter->second.x;
		candidateP._y = iter->second.y;
		candidateP._z = iter->second.z;

		CVector newVec;
		Common::CalNormalVector(pi, candidateP, pj, newVec);
		float angle_cos = CVector::vectorInnerProduct(vec, newVec) / (vec.vectorMag(vec) * newVec.vectorMag(newVec));
		if (angle_cos > 1)
			angle_cos = 1;
		else if (angle_cos < -1)
			angle_cos = -1;

		if (acos(angle_cos) > 120.0f)
			iter = result.erase(iter);
		else
			iter++;
	}

}


//选择最佳点
void Common::FindBestPoint(Point pi, Point pj, Point pk, std::vector<std::pair< double, pcl::PointXYZ>> &result, CLine CurrentE, std::list<CFace> &ST, std::vector<CLine> InnerE, std::map<Point, bool> &flag, std::vector<Point> &bestPoints)
{
	// TODO: 在此处添加实现代码.
	CVector vec; //当前活动边所在三角形的法向量
	Common::CalNormalVector(pi, pj, pk, vec);

	const float PI = 3.1415926;
	float angle_min, angle_max; //候选三角形的角度的最大值与最小值
	float a = CLine::LineLength_Point(pi, pj);
	float b = CLine::LineLength_Point(pj, pk);
	float c = CLine::LineLength_Point(pk, pi);
	float maxSide = max(a, b);
	maxSide = max(maxSide, c);
	if (a == maxSide)
		angle_max = acos((pow(b, 2) + pow(c, 2) - pow(a, 2)) / (2 * b * c)) * 180 / PI;
	else if (b == maxSide)
		angle_max = acos((pow(a, 2) + pow(c, 2) - pow(b, 2)) / (2 * a * c)) * 180 / PI;
	else
		angle_max = acos((pow(b, 2) + pow(a, 2) - pow(c, 2)) / (2 * b * a)) * 180 / PI;
	float minSide = min(a, b);
	minSide = min(minSide, c);
	if (a == minSide)
		angle_min = acos((pow(b, 2) + pow(c, 2) - pow(a, 2)) / (2 * b * c)) * 180 / PI;
	else if (b == minSide)
		angle_min = acos((pow(a, 2) + pow(c, 2) - pow(b, 2)) / (2 * a * c)) * 180 / PI;
	else
		angle_min = acos((pow(b, 2) + pow(a, 2) - pow(c, 2)) / (2 * b * a)) * 180 / PI;

	std::vector<std::pair<float, Point>> joinCosts;
	float angle_cur;//当前备选三角片的最大内角,大边对大角
	//为候选点添加代价
	for (auto it = result.begin(); it != result.end(); it++)
	{
		Point candidateP;
		candidateP._x = it->second.x;
		candidateP._y = it->second.y;
		candidateP._z = it->second.z;
		float joinCost; //节点candidateP的代价
		float cost_angle1, cost_angle2;
		//候选三角形的法向量
		CVector newVec;
		CalNormalVector(pi, candidateP, pj, newVec);
		float angle_cos; 
		float temp = vec.vectorInnerProduct(vec, newVec);
		float temp1 = (vec.vectorMag(vec) * vec.vectorMag(newVec));
		angle_cos = vec.vectorInnerProduct(vec, newVec) / (vec.vectorMag(vec) * vec.vectorMag(newVec));
		cost_angle1 = sqrt(1 - angle_cos * angle_cos);

		a = CLine::LineLength_Point(pi, candidateP);
		b = CLine::LineLength_Point(candidateP, pj);
		c = CLine::LineLength_Point(pj, pi);
		maxSide = max(a, b);
		maxSide = max(maxSide, c);		
		if (a == maxSide)
			angle_cur = acos((pow(b, 2) + pow(c, 2) - pow(a, 2)) / (2 * b * c)) * 180 / PI;
		else if (b == maxSide)
			angle_cur = acos((pow(a, 2) + pow(c, 2) - pow(b, 2)) / (2 * a * c)) * 180 / PI;
		else
			angle_cur = acos((pow(b, 2) + pow(a, 2) - pow(c, 2)) / (2 * b * a)) * 180 / PI;
		
		cost_angle2 = abs((angle_cur - angle_min) / (angle_max - angle_min));
		joinCost = cost_angle1 + cost_angle2;
		joinCosts.push_back(pair<float, Point>(joinCost, Point(it->second.x,it->second.y,it->second.z)));
	}
	
	std::sort(joinCosts.begin(), joinCosts.end(), [&](const std::pair<float, Point> &pair1, const std::pair<float, Point> &pair2) {return (((pair1.first < pair2.first) ? true : false)); });
	
	//候选三角片质量检测
	Point bestP;
	for (auto it = joinCosts.begin(); it != joinCosts.end(); it++)
	{
		Point pc;
		pc._x = it->second._x;
		pc._y = it->second._y;
		pc._z = it->second._z;

		//检测非流行点
		CFace curFace(pi, pj, pc);
		if (findNearFace_Point(curFace, pc, ST))
			continue;

		//检测非流行边
		CLine p_ic(pi, pc);
		vector<CLine>::iterator its = find(InnerE.begin(), InnerE.end(), p_ic);
		if (its != InnerE.end())
			continue;
		CLine p_cj(pc, pj);
		its = find(InnerE.begin(), InnerE.end(), p_cj);
		if (its != InnerE.end())
			continue;

		//检测自交三角形
		if (IntersectTriangle(pi, pj, pc, ST))
			continue;
		else
		{
			bestP = pc;
			bestPoints.push_back(bestP);
		}			
	}

	////检测是否有冗余点
	//CVector vec_new;
	//CalNormalVector(pi, pj, bestP, vec_new);

	//a = vec_new.GetX();
	//b = vec_new.GetY();
	//c = vec_new.GetZ();
	//float d = -a * pi._x - b * pi._y - c * pi._z;

	//for (auto it = result.begin(); it != result.end(); it++)
	//{
	//	if ((it->second.x == bestP._x) 
	//	&& (it->second.y == bestP._y)
	//	&& (it->second.z == bestP._z))
	//		continue;
	//	Point  spatialpoint;;
	//	spatialpoint._x = it->second.x;
	//	spatialpoint._y = it->second.y;
	//	spatialpoint._z = it->second.z;
	//	Point subpoint; //投影点
	//	subpoint._x = ((b * b + c * c) * spatialpoint._x - a * (b * spatialpoint._y + c * spatialpoint._z + d)) / (a * a + b * b + c * c);
	//	subpoint._y = b / a * (subpoint._x - spatialpoint._x) + spatialpoint._y;
	//	subpoint._z = c / a * (subpoint._x - spatialpoint._x) + spatialpoint._z;

	//	//判断投影点是否在三角形内
	//	if (TriangleIncludeSubpoint(pi, pj, bestP, subpoint))
	//		flag[spatialpoint] = true; //当前点设置为自由点/排除点
	//}
}


// 寻找点的所有邻接三角形,并判断是否共边
bool Common::findNearFace_Point(CFace curFace, Point pc, std::list<CFace> ST)
{
	// TODO: 在此处添加实现代码.
	std::set<Point> points;
	points.insert(curFace.GetPoint1());
	points.insert(curFace.GetPoint2());
	points.insert(curFace.GetPoint3());

	int count = 0; // 记录点pc的邻接三角形的个数
	for (auto it = ST.begin(); it != ST.end(); it++)
	{
		//是点pc的邻接三角形		
		if (((it->GetPoint1()._x == pc._x) && (it->GetPoint1()._y == pc._y) && (it->GetPoint1()._z == pc._z)) 
			|| ((it->GetPoint2()._x == pc._x) && (it->GetPoint2()._y == pc._y) && (it->GetPoint2()._z == pc._z))
			|| ((it->GetPoint3()._x == pc._x) && (it->GetPoint3()._y == pc._y) && (it->GetPoint3()._z == pc._z)))
		{
			count++;
			points.insert(it->GetPoint1());
			points.insert(it->GetPoint2());
			points.insert(it->GetPoint3());
			if (points.size() <= 4)
				return false;
		}			
	}

	if (count == 0)
		return false;

	return true;
}


// 三角形pi_pj_pc与pi，pj构成的边的一环邻接三角形，这两个三角形是否相交
bool Common::IntersectTriangle(Point pi, Point pj, Point pc, std::list<CFace> ST)
{
	// TODO: 在此处添加实现代码.
	std::vector<CLine> lines;
	for (auto it = ST.begin(); it != ST.end(); it++)
	{
		if (it->GetPoint1() == pi)
		{
			CLine line(it->GetPoint2(), it->GetPoint3());
			lines.push_back(line);
		}
		else if (it->GetPoint2() == pi)
		{
			CLine line(it->GetPoint3(), it->GetPoint1());
			lines.push_back(line);
		}
		else if (it->GetPoint3() == pi)
		{
			CLine line(it->GetPoint1(), it->GetPoint2());
			lines.push_back(line);
		}

		if (it->GetPoint1() == pj)
		{
			CLine line(it->GetPoint2(), it->GetPoint3());
			lines.push_back(line);
		}
		else if (it->GetPoint2() == pj)
		{
			CLine line(it->GetPoint3(), it->GetPoint1());
			lines.push_back(line);
		}
		else if (it->GetPoint3() == pj)
		{
			CLine line(it->GetPoint1(), it->GetPoint2());
			lines.push_back(line);
		}
	}

	for (auto it = lines.begin(); it != lines.end(); it++)
	{
		if (it->getPointStart() == pj || it->getPointEnd() == pj)
			continue;
		else
		{
			if (IntersectionLine(it->getPointStart(), it->getPointEnd(), pc, pi)
				|| IntersectionLine(it->getPointStart(), it->getPointEnd(), pc, pj))
				return true;
		}
	}
	return false;
}


// 判断两条线段是否相交
bool Common::LineInterset(CLine &line1, CLine &line2)
{
	// 检测两直线相交（向量AB x AC 与AB x AD 符号相反）
	Point line1S, line1E, line2S, line2E;
	line1S = line1.getPointStart();
	line1E = line1.getPointEnd();
	CVector vectorAB(line1S._x - line1E._x, line1S._y - line1E._y, line1S._z - line1E._z);
	CVector vectorAC(line1S._x - line2S._x, line1S._y - line2S._y, line1S._z - line2S._z);
	CVector  vectorAD(line1S._x - line2E._x, line1S._y - line2E._y, line1S._z - line2E._z);

	CVector value1 = vectorAB * vectorAC;
	CVector value2 = vectorAB * vectorAD;
	float v1 = value1.GetX() + value1.GetY() + value1.GetZ();
	float v2 = value2.GetX() + value2.GetY() + value2.GetZ();

	if ((v1 > 0 && v2 > 0) || (v1 < 0 && v2 < 0))
		return true;
	return false;

}


// 判断两条线段是否相交
bool Common::IntersectionLine(Point pi, Point pj, Point pc, Point pk)
{
	// TODO: 在此处添加实现代码.
	CVector p_ic(pc._x - pi._x, pc._y - pi._y, pc._z - pi._z);
	CVector p_ik(pk._x - pi._x, pk._y - pi._y, pk._z - pi._z);
	CVector p_jc(pc._x - pj._x, pc._y - pj._y, pc._z - pj._z);
	CVector p_jk(pk._x - pj._x, pk._y - pj._y, pk._z - pj._z);

	CVector p_ci(pi._x - pc._x, pi._y - pc._y, pi._z - pc._z);
	CVector p_ki(pi._x - pk._x, pi._y - pk._y, pi._z - pk._z);
	CVector p_cj(pj._x - pc._x, pj._y - pc._y, pj._z - pc._z);
	CVector p_kj(pj._x - pk._x, pj._y - pk._y, pj._z - pk._z);

	float sin_value = p_ic.MultiplicationCross(p_ik) * p_jc.MultiplicationCross(p_jk);
	float sin_value1 = p_ci.MultiplicationCross(p_cj) * p_ki.MultiplicationCross(p_kj);
	if (sin_value > 0 && sin_value1 > 0)
		return true;
	return false;
}


// 判断前面三个点构成的三角形是否包含第四个点
bool Common::TriangleIncludeSubpoint(Point pi, Point pj, Point bestP, Point subpoint)
{
	// TODO: 在此处添加实现代码.
	float area = TriangleArea(pi, pj, bestP);
	float area_new = TriangleArea(pi, pj, subpoint) + TriangleArea(pi, bestP, subpoint) + TriangleArea(pj, bestP, subpoint);
	if (area_new > area)
		return false;
	else
		return true;
}


// 得到三角形ABC的面积
float Common::TriangleArea(Point A, Point B, Point C)
{
	// TODO: 在此处添加实现代码.
	float s = ((B._y - A._y) * (C._z - A._z) + (B._z - A._z) * (C._x - A._x) + (B._x - A._x) * (C._y - A._y)) -
		((C._y - A._y) * (B._z - A._z) + (C._z - A._z) * (B._x - A._x) + (C._x - A._x) * (B._y - A._y));
	return 0.5 * s;
}


// 更新活动链表
void Common::UpdateActiveList(std::vector<CLine> &ActiveE, CLine CurrentE, Point bestP, std::vector<CLine> &InnerE, std::vector<Point> &FreeP, std::vector<Point> &ActiveP, std::map<Point, bool> &flag, std::list<CFace> &ST, std::vector<std::pair< double, pcl::PointXYZ>> &result)
{
	// 判断最佳点添加的位置
	CFace f(CurrentE.getPointStart(), CurrentE.getPointEnd(), bestP);
	ST.push_back(f); //添加三角面片
	AddFixedPoint(flag, f, result);
	int type = BestPositionType(ActiveE, CurrentE, bestP, FreeP);
	//int type = BestPositionType(ActiveE, CurrentE, bestP, FreeP, ActiveP, InnerE,flag);

	switch (type)
	{
	case 0:
		UpdateMode(ActiveE, CurrentE, bestP, InnerE, FreeP, ActiveP); 
		break;
	case 1:
		UpdateMode1(ActiveE, CurrentE, bestP, InnerE, ActiveP, flag);
		break;
	case 2:
		UpdateMode2(ActiveE, CurrentE, bestP, InnerE, ActiveP, flag);
		break;
	case 3:
		UpdateMode3(ActiveE, CurrentE, bestP, InnerE, ActiveP, flag, ST, result);
		break;
	default:
		break;
	}
}


// 判断最佳点添加的位置类型
/*
0 -- 最佳点是自由点
1 -- 最佳点位于活动边上且为当前活动边前相邻边的端点
2 -- 最佳点位于活动边上且为当前活动边后相邻边的端点
3 -- 最佳点位于活动边上且与当前活动边没有相邻关系
*/
int Common::BestPositionType(std::vector<CLine> ActiveE, CLine CurrentE, Point bestP, std::vector<Point> FreeP)
{
	// TODO: 在此处添加实现代码.
	std::vector<Point>::iterator it = find(FreeP.begin(), FreeP.end(), bestP);
	if (it != FreeP.end()) //最佳点是自由点
		return 0;
	std::vector<CLine>::iterator its = find(ActiveE.begin(), ActiveE.end(), CurrentE);
	std::vector<CLine>::iterator front;
	std::vector<CLine>::iterator behind;
	if (its == ActiveE.begin())
	{
		front = ActiveE.end();
		behind = its++;
	}
	else if (its == ActiveE.end())
	{
		front = its--;
		behind = ActiveE.begin();
	}
	else
	{
		front = (its--);
		its++;
		behind = (++its);
	}
	if (front->getPointStart() == bestP)
		return 1;
	if (front->getPointEnd() == bestP)
		return 2;
	return 3;
}

int Common::BestPositionType(std::vector<CLine> &ActiveE, CLine &CurrentE, const Point &bestP, std::vector<Point> &FreeP, std::vector<Point> &ActiveP, std::vector<CLine> &InnerE, std::map<Point, bool> &flag)
{
	int num = 0;
	int j = 0;
	std::vector<Point>::iterator it;
	it = std::find(FreeP.begin(), FreeP.end(), bestP);
	if (it != FreeP.end()) //最佳点是自由点
		return 0;
	it = std::find(ActiveP.begin(), ActiveP.end(), bestP);
	if (it != ActiveP.end() ){
		
		
		auto iter = std::find(ActiveE.begin(), ActiveE.end(), CLine(bestP, CurrentE.getPointStart()));
		if (iter != ActiveE.end())
		{
			j++;
			num = 1;
		}
		iter = std::find(ActiveE.begin(), ActiveE.end(), CLine(CurrentE.getPointEnd(), bestP));
		if (iter != ActiveE.end())
		{
			j++;
			num = 2;
		}
	}

	if (j == 1)		//	最佳点在活动边上，情况二
	{
		if (num == 1)
		{
			//	更新活动边
			auto it = std::find(ActiveE.begin(), ActiveE.end(), CurrentE);
			it = ActiveE.erase(it);
			it = std::find(ActiveE.begin(), ActiveE.end(), CLine(bestP, CurrentE.getPointStart()));
			it = ActiveE.erase(it);
			ActiveE.insert(it, CLine(bestP, CurrentE.getPointEnd()));
			//	更新活动点
			ActiveP.erase(std::find(ActiveP.begin(), ActiveP.end(), CurrentE.getPointStart()));
			//	更新固定点
			auto iter = flag.find(CurrentE.getPointStart());
			iter->second = true;
			//	更新固定边
			InnerE.push_back(CurrentE);
			InnerE.push_back(CLine(bestP, CurrentE.getPointStart()));


		}
		else
		{
			auto it = std::find(ActiveE.begin(), ActiveE.end(), CurrentE);
			it = ActiveE.erase(it);
			it = std::find(ActiveE.begin(), ActiveE.end(), CLine(CurrentE.getPointEnd(), bestP));
			it = ActiveE.erase(it);
			ActiveE.insert(it, CLine(CurrentE.getPointStart(), bestP));
			//	更新活动点
			ActiveP.erase(std::find(ActiveP.begin(), ActiveP.end(), CurrentE.getPointEnd()));
			//	更新固定点
			auto iter = flag.find(CurrentE.getPointEnd());
			iter->second = true;
			//	更新固定边
			InnerE.push_back(CurrentE);
			InnerE.push_back(CLine(CurrentE.getPointEnd(), bestP));
		}

		return 7;
	}
	else if (j == 0)	//	孔洞情况
		return 3;
	else
	{
		return -1;			//	错误情况发生

	}
}


/*************************更新活动链表*************************************/
// 最佳点是自由点
void Common::UpdateMode(std::vector<CLine> &ActiveE, CLine CurrentE, Point bestP, std::vector<CLine> &InnerE, std::vector<Point> &FreeP, std::vector<Point> &ActiveP)
{
	// TODO: 在此处添加实现代码.
	InnerE.push_back(CurrentE);
	std::vector<Point>::iterator its = find(FreeP.begin(), FreeP.end(), bestP);
	FreeP.erase(its);
	ActiveP.push_back(bestP);

	std::vector<CLine>::iterator it = find(ActiveE.begin(), ActiveE.end(), CurrentE);
	CLine line1(CurrentE.getPointStart(), bestP);
	CLine line2(bestP, CurrentE.getPointEnd());
	it = ActiveE.erase(it);
	it = ActiveE.insert(it, line2);
	it = ActiveE.insert(it, line1);
}


// 最佳点位于活动边上且为当前活动边前相邻边的端点
void Common::UpdateMode1(std::vector<CLine> &ActiveE, CLine CurrentE, Point bestP, std::vector<CLine> &InnerE, std::vector<Point> &ActiveP, std::map<Point, bool> &flag)
{
	// TODO: 在此处添加实现代码.
	std::vector<CLine>::iterator it = find(ActiveE.begin(), ActiveE.end(), CurrentE);
	InnerE.push_back(*it);
	if (it == ActiveE.begin())
	{
		InnerE.push_back(*ActiveE.end());
		ActiveE.erase(ActiveE.end());
	}

	std::vector<Point>::iterator its = find(ActiveP.begin(), ActiveP.end(), CurrentE.getPointStart());
	ActiveP.erase(its);

	CLine line(bestP, CurrentE.getPointEnd());
	ActiveE.insert(it, line);
	ActiveE.erase(it);
}


// 最佳点位于活动边上且为当前活动边后相邻边的端点
void Common::UpdateMode2(std::vector<CLine> &ActiveE, CLine CurrentE, Point bestP, std::vector<CLine> &InnerE, std::vector<Point> &ActiveP, std::map<Point, bool> &flag)
{
	// TODO: 在此处添加实现代码.
	std::vector<CLine>::iterator it = find(ActiveE.begin(), ActiveE.end(), CurrentE);
	InnerE.push_back(*it);
	if (it == ActiveE.end())
	{
		InnerE.push_back(*ActiveE.begin());
		ActiveE.erase(ActiveE.begin());
	}

	std::vector<Point>::iterator its = find(ActiveP.begin(), ActiveP.end(), CurrentE.getPointEnd());
	ActiveP.erase(its);

	CLine line(CurrentE.getPointStart(), bestP);
	ActiveE.insert(it, line);
	ActiveE.erase(it);
}


//最佳点位于活动边上且与当前活动边没有相邻关系
void Common::UpdateMode3(std::vector<CLine> &ActiveE, CLine CurrentE, Point bestP, std::vector<CLine> &InnerE, std::vector<Point> &ActiveP, std::map<Point, bool> &flag, std::list<CFace> &ST, std::vector<std::pair< double, pcl::PointXYZ>> &result)
{
	// TODO: 在此处添加实现代码.
	std::vector<CLine>::iterator it = find(ActiveE.begin(), ActiveE.end(), CurrentE);
	InnerE.push_back(*it);
	std::vector<CLine>::iterator its = it;
	int holeSide_count = 1;
	if (its == ActiveE.end())
		its = ActiveE.begin();
	else
		its++;
	while (its->getPointEnd() != bestP)
	{
		InnerE.push_back(*its);
		std::vector<Point>::iterator iter = find(ActiveP.begin(), ActiveP.end(), its->getPointEnd());
		iter = ActiveP.erase(iter);
		flag[its->getPointEnd()] = true;
		if (holeSide_count > 1) //可以构造三角形
		{
			CFace f(its->getPointStart(), its->getPointEnd(), CurrentE.getPointEnd());
			ST.push_back(f);
			Common::AddFixedPoint(flag, f, result);
			CLine l(its->getPointEnd(), CurrentE.getPointEnd());
			InnerE.push_back(l);
		}
		holeSide_count++;
		if (its == ActiveE.end())
			its = ActiveE.begin();
		else
			its++;
	}
	if (holeSide_count > 1)
	{
		CFace ff(its->getPointStart(), bestP, CurrentE.getPointEnd());
		ST.push_back(ff);
		Common::AddFixedPoint(flag, ff, result);
	}
	else
		InnerE.push_back(*its);

	its = ActiveE.erase(it, its);
	its = ActiveE.erase(its);
	CLine l(CurrentE.getPointStart(), bestP);
	ActiveE.insert(its, l);
}

//>>>>>>> origin/dev_hhy

//去除排除点
void Common::AddFixedPoint(std::map<Point, bool> &flag, CFace face, std::vector<std::pair< double, pcl::PointXYZ>> &result)
{
	// TODO: 在此处添加实现代码.
	Point p1 = face.GetPoint1();
	Point p2 = face.GetPoint2();
	Point p3 = face.GetPoint3();
	for (auto it = result.begin(); it != result.end(); it++)
	{
		Point p(it->second.x, it->second.y, it->second.z);
		//判断点是否在面片face中
		if (p != p1 && p != p2 && p != p3 && TriangleIncludeSubpoint(p1, p2, p3, p))
			flag[p] = true;
	}
}
