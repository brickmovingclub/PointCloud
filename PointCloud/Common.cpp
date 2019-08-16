#include "stdafx.h"
#include "CLine.h"

#include "CVector.h"

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
	int i = 0;
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
	}
	return (i == nearPoints.size() ? true : false);
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

void Common::NearRadiusSearch(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, const pcl::PointXYZ &point, const float &radius, std::vector<std::pair< double,pcl::PointXYZ>> &nearPoint)
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
			std::cout << "    " << cloud->points[pointIdxRadiusSearch[i]].x
			<< " " << cloud->points[pointIdxRadiusSearch[i]].y
			<< " " << cloud->points[pointIdxRadiusSearch[i]].z
			<< " (squared distance: " << pointRadiusSquaredDistance[i] << ")" << std::endl;
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

void Common::PCLDrawLine(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, pcl::visualization::PCLVisualizer::Ptr viewer)
{
	//cloud = getpoint();//实时获取点云
	pcl::PointXYZ  minPt, maxPt;
	//viewer->removeAllShapes();
	pcl::getMinMax3D(*cloud, minPt, maxPt);

	pcl::PointXYZ origin(0, 0, 0);

	viewer->setBackgroundColor(125, 125, 125); //背景色

	viewer->addLine<pcl::PointXYZ>(origin, minPt, 255, 0, 0, "line1"); //红色线段,线的名字叫做"line1

	viewer->spinOnce(100);
}


std::vector<int> Common::findCandidatePoints(pcl::PointCloud<pcl::PointXYZ>::Ptr _cloud, Point pi, Point pj, Point pk, std::vector<bool> flag,std::vector<CLine> ActiveE, CLine CurrentE)
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

	Point pm;
	pm._x = (pi._x + pj._x) / 2;
	pm._y = (pi._y + pj._y) / 2;
	pm._z = (pi._z + pj._z) / 2;

	//计算点pm的r范围内的领域点集
	std::vector<int> near_pm; 


	//精简领域点集
	int size = near_pm.size();
	for (int i = 0; i < size; i++)
	{
		//删除固定点与排除点
		if (flag[near_pm[i]])
			near_pm.erase(near_pm.begin() + i);
	}
	if(near_pm.empty()) //边pi-pj是边界边
		return near_pm;

	//找到pi,pj相邻的活动点
	std::vector<CLine>::iterator it = std::find(ActiveE.begin(), ActiveE.end(), CurrentE);
	it--;
	Point pa = it->getPointStart();
	it++; it++;
	Point pb = it->getPointEnd();

	//边角度约束简化
	CVector p_ia(pa._x - pi._x, pa._y - pi._y, pa._z - pi._z);
	CVector p_ij(pj._x - pi._x, pj._y - pi._y, pj._z - pi._z);
	CVector p_jb(pb._x - pj._x, pb._y - pj._y, pb._z - pj._z);
	CVector p_ji(pi._x - pj._x, pi._y - pj._y, pi._z - pj._z);

	float angle_ia_ij = p_ia.vectorInnerProduct(p_ia, p_ij) / (p_ia.vectorMag(p_ia) * p_ia.vectorMag(p_ij));
	float angle_jb_ji = p_jb.vectorInnerProduct(p_jb, p_ji) / (p_jb.vectorMag(p_jb) * p_jb.vectorMag(p_ji));
	
	float angleA = max(angle_ia_ij, (float)pow(2, 0.5) / 2);
	float angleB = max(angle_jb_ji, (float)pow(2, 0.5) / 2);

	float angle = 0.0f;
	size = near_pm.size();
	for (int i = 0; i < size; i++)
	{
		Point candidateP;
		candidateP._x = _cloud->points[near_pm[i]].x;
		candidateP._y = _cloud->points[near_pm[i]].y;
		candidateP._z = _cloud->points[near_pm[i]].z;

		CVector p_icandidateP(candidateP._x - pi._x, candidateP._y - pi._y, candidateP._z - pi._z);
		CVector p_jcandidateP(candidateP._x - pj._x, candidateP._y - pj._y, candidateP._z - pj._z);

		angle = p_ia.vectorInnerProduct(p_ia, p_icandidateP) / (p_ia.vectorMag(p_ia) * p_ia.vectorMag(p_icandidateP));
		if (angle < angleA)
		{
			near_pm.erase(near_pm.begin() + i);
			continue;
		}
		else
		{
			angle = p_jb.vectorInnerProduct(p_jb, p_jcandidateP) / (p_jb.vectorMag(p_jb) * p_jb.vectorMag(p_jcandidateP));
			if (angle < angleB)
				near_pm.erase(near_pm.begin() + i);
		}		
	}
	if (near_pm.empty()) //边pi-pj是边界边
		return near_pm;

	//面角度约束简化
	CVector vec;
	vec = vec.GetNormal(pi, pj, pk);
	size = near_pm.size();
	for (int i = 0; i < size; i++)
	{
		Point candidateP;
		candidateP._x = _cloud->points[near_pm[i]].x;
		candidateP._y = _cloud->points[near_pm[i]].y;
		candidateP._z = _cloud->points[near_pm[i]].z;

		CVector newVec = vec.GetNormal(pi, candidateP, pj);
		angle = vec.vectorInnerProduct(vec, newVec) / (vec.vectorMag(vec) * vec.vectorMag(newVec));
		if (angle < pow(3, 0.5) / 2)
			near_pm.erase(near_pm.begin() + i);
	}
	return near_pm; //若为空，则pi-pj是边界边
}
