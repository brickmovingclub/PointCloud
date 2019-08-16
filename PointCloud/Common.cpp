#include "stdafx.h"

#include "CVector.h"

#include "Common.h"


Common::Common()
{
}


Common::~Common()
{
}


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
	//	��ռ�ƽ��ķ�����
	float dx, dy, dz;
	//float a1,a2,b1,b2,c1,c2;
	//
	float _a1 = pj.x - pi.x; float _b1 = pj.y - pi.y; float _c1 = pj.z - pi.z;
	float _a2 = pk.x - pi.x; float _b2 = pk.y - pi.y; float _c2 = pk.z - pi.z;
	float _a3 = pk.x - pj.x; float _b3 = pk.y - pj.y; float _c3 = pk.z - pj.z;
	//_a1x+_b1y+_c1z=0;
	//_a2x+_b2y+_c2z=0;
	//_a3x+_b3y+_c3z=0;
	//3��δ֪��3��������ɵ���η��������
	//ϵ������A
	//| _a1 _b1 _c1 |
	//| _a2 _b2 _c2 |
	//| _a3 _b3 _c3 |
	//�������ʽA��ֵ������0������Ψһ����Ϊ���
	float DA = _a1 * _b2*_c3 + _b1 * _c2*_a3 + _a2 * _b3*_c1 - _a3 * _b2*_c1 - _a1 * _b3*_c2 - _a2 * _b1*_c3;
	if (DA != 0)
	{
		dx = 0.0f;
		dy = 0.0f;
		dz = 0.0f;
		return false;
	}
	//---------------------------------------------//
	//�������ʽA��ֵ����0�����з����
	//����⼴x!=0ʱ�н����y!=0ʱ�н����z!=0ʱ�н�
	float x = 0.0f, y = 0.0f, z = 0.0f;
	//��z!=0ʱ�н�,ȡz=-1
	//_a1x+_b1y=_c1;---(1)
	//_a2x+_b2y=_c2;---(2)
	//_a3x+_b3y=_c3;---(3)
	//��ȡ2�����̼��ɣ��ڴ�ȡ(1)(2)
	x = 0.0f; y = 0.0f;
	bool flag3 = GetTwoLineIntersection(_a1, _b1, _c1, _a2, _b2, _c2, x, y);
	if (flag3)//�������
	{
		dx = -x;
		dy = -y;
		dz = 1.0f;
		return true;
	}
	//���費����������������һ������
	//��x!=0ʱ�н�ȡx=-1��ƽ��������ֱ���󽻵�����
	//_b1y+_c1z=_a1;---(1)
	//_b2y+_c2z=_a2;---(2)
	//_b3y+_c3z=_a3;---(3)
	//��ȡ2�����̼��ɣ��ڴ�ȡ(1)(2)
	y = 0.0f; z = 0.0f;
	bool flag1 = GetTwoLineIntersection(_b1, _c1, _a1, _b2, _c2, _a2, y, z);
	if (flag1)//�������
	{
		dx = 1.0f;
		dy = -y;
		dz = -z;
		return true;
	}
	//���費����������������һ������
	//��y!=0ʱ�н�ȡy=-1��ƽ��������ֱ���󽻵�����
	//_a1x+_c1z=_b1;---(1)
	//_a2x+_c2z=_b2;---(2)
	//_a3x+_c3z=_b3;---(3)
	//��ȡ2�����̼��ɣ��ڴ�ȡ(1)(2)
	x = 0.0f; z = 0.0f;
	bool flag2 = GetTwoLineIntersection(_a1, _c1, _b1, _a2, _c2, _b2, x, z);
	if (flag2)//�������
	{
		dx = -x;
		dy = 1.0f;
		dz = -z;
		return true;
	}

	//���м��趼�����������ʧ��
	return false;
}
*/
bool Common::GetTwoLineIntersection(float _a1, float _b1, float _c1, float _a2, float _b2, float _c2, float &x, float &y)
{
	//_a1x+_b1y=_c1;---(1)
	//_a2x+_b2y=_c2;---(2)
	//
	if (_c1 == 0 && _c2 == 0)
	{
		//2��δ֪��2��������ɵ���η��������
		//ϵ������B
		//| _a1 _b1 |
		//| _a2 _b2 |
		float DB = _a1 * _b2 - _a2 * _b1;
		if (DB != 0)//��Ψһ���
		{
			x = 0;
			y = 0;
			return true;
		}
		else//��������
		{
			x = 0;
			y = 0;
			return false;
		}
	}
	else
	{
		//2��δ֪��2��������ɵķ���η��������
		//ϵ������B
		//| _a1 _b1 |
		//| _a2 _b2 |
		//
		float DB = _a1 * _b2 - _a2 * _b1;
		if (DB != 0)//��Ψһ��
		{
			float dD1 = _c1 * _b2 - _c2 * _b1;
			float dD2 = _a1 * _c2 - _a2 * _c1;
			x = dD1 / DB;
			y = dD2 / DB;
			return true;
		}
		else//������������޽�
		{
			x = 0;
			y = 0;
			return false;
		}
	}
	return false;
}

bool Common::CalNormalVector(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3, float &dx, float &dy, float &dz)
{
	//float a1,a2,b1,b2,c1,c2;
	//
	float _a1 = x2 - x1; float _b1 = y2 - y1; float _c1 = z2 - z1;
	float _a2 = x3 - x1; float _b2 = y3 - y1; float _c2 = z3 - z1;
	float _a3 = x3 - x2; float _b3 = y3 - y2; float _c3 = z3 - z2;
	//_a1x+_b1y+_c1z=0;
	//_a2x+_b2y+_c2z=0;
	//_a3x+_b3y+_c3z=0;
	//3��δ֪��3��������ɵ���η��������
	//ϵ������A
	//| _a1 _b1 _c1 |
	//| _a2 _b2 _c2 |
	//| _a3 _b3 _c3 |
	//�������ʽA��ֵ������0������Ψһ����Ϊ���
	float DA = _a1 * _b2*_c3 + _b1 * _c2*_a3 + _a2 * _b3*_c1 - _a3 * _b2*_c1 - _a1 * _b3*_c2 - _a2 * _b1*_c3;
	if (DA != 0)
	{
		dx = 0.0f;
		dy = 0.0f;
		dz = 0.0f;
		return false;
	}
	//---------------------------------------------//
	//�������ʽA��ֵ����0�����з����
	//����⼴x!=0ʱ�н����y!=0ʱ�н����z!=0ʱ�н�
	float x = 0.0f, y = 0.0f, z = 0.0f;
	//��z!=0ʱ�н�,ȡz=-1
	//_a1x+_b1y=_c1;---(1)
	//_a2x+_b2y=_c2;---(2)
	//_a3x+_b3y=_c3;---(3)
	//��ȡ2�����̼��ɣ��ڴ�ȡ(1)(2)
	x = 0.0f; y = 0.0f;
	bool flag3 = GetTwoLineIntersection(_a1, _b1, _c1, _a2, _b2, _c2, x, y);
	if (flag3)//�������
	{
		dx = -x;
		dy = -y;
		dz = 1.0f;
		return true;
	}
	//���費����������������һ������
	//��x!=0ʱ�н�ȡx=-1��ƽ��������ֱ���󽻵�����
	//_b1y+_c1z=_a1;---(1)
	//_b2y+_c2z=_a2;---(2)
	//_b3y+_c3z=_a3;---(3)
	//��ȡ2�����̼��ɣ��ڴ�ȡ(1)(2)
	y = 0.0f; z = 0.0f;
	bool flag1 = GetTwoLineIntersection(_b1, _c1, _a1, _b2, _c2, _a2, y, z);
	if (flag1)//�������
	{
		dx = 1.0f;
		dy = -y;
		dz = -z;
		return true;
	}
	//���費����������������һ������
	//��y!=0ʱ�н�ȡy=-1��ƽ��������ֱ���󽻵�����
	//_a1x+_c1z=_b1;---(1)
	//_a2x+_c2z=_b2;---(2)
	//_a3x+_c3z=_b3;---(3)
	//��ȡ2�����̼��ɣ��ڴ�ȡ(1)(2)
	x = 0.0f; z = 0.0f;
	bool flag2 = GetTwoLineIntersection(_a1, _c1, _b1, _a2, _c2, _b2, x, z);
	if (flag2)//�������
	{
		dx = -x;
		dy = 1.0f;
		dz = -z;
		return true;
	}

	//���м��趼�����������ʧ��
	return false;
}


//=======
// �ж������Ƿ���ͬһֱ�����Լ�������Ƿ������㹹�ɵ�Բ����
bool Common::Condition_a_b(pcl::PointXYZ pi, pcl::PointXYZ pj, pcl::PointXYZ pk, std::vector<std::pair<double,pcl::PointXYZ>> &near_pi)
{
	// TODO: �ڴ˴����ʵ�ִ���.	
	double a1, b1, c1, d1;
	double a2, b2, c2, d2;
	double a3, b3, c3, d3;

	double x1 = pi.x, y1 = pi.y, z1 = pi.z;
	double x2 = pj.x, y2 = pj.y, z2 = pj.z;
	double x3 = pk.x, y3 = pk.y, z3 = pk.z;

	//�ж������Ƿ���
	double s = 0.5 * (x1 * y2 + x2 * y3 + x3 * y1 - x1 * y3 - x2 * y1 - x3 * y2);
	if (s == 0)
		return false;

	//��Բ��
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

	//Բ������
	pcl::PointXYZ centerpoint;
	centerpoint.x = -(b1*c2*d3 - b1 * c3*d2 - b2 * c1*d3 + b2 * c3*d1 + b3 * c1*d2 - b3 * c2*d1)
		/ (a1*b2*c3 - a1 * b3*c2 - a2 * b1*c3 + a2 * b3*c1 + a3 * b1*c2 - a3 * b2*c1);
	centerpoint.y = (a1*c2*d3 - a1 * c3*d2 - a2 * c1*d3 + a2 * c3*d1 + a3 * c1*d2 - a3 * c2*d1)
		/ (a1*b2*c3 - a1 * b3*c2 - a2 * b1*c3 + a2 * b3*c1 + a3 * b1*c2 - a3 * b2*c1);
	centerpoint.z = -(a1*b2*d3 - a1 * b3*d2 - a2 * b1*d3 + a2 * b3*d1 + a3 * b1*d2 - a3 * b2*d1)
		/ (a1*b2*c3 - a1 * b3*c2 - a2 * b1*c3 + a2 * b3*c1 + a3 * b1*c2 - a3 * b2*c1);

	//�뾶
	double dR = sqrt(pow(centerpoint.x - pi.x, 2) + pow(centerpoint.y - pi.y, 2) + pow(centerpoint.z - pi.z, 2));

	//�ж�������Ƿ������㹹�ɵ�Բ����
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
	//cloud = getpoint();//ʵʱ��ȡ����
	pcl::PointXYZ  minPt, maxPt;
	//viewer->removeAllShapes();
	pcl::getMinMax3D(*cloud, minPt, maxPt);

	pcl::PointXYZ origin(0, 0, 0);

	viewer->setBackgroundColor(125, 125, 125); //����ɫ

	viewer->addLine<pcl::PointXYZ>(origin, minPt, 255, 0, 0, "line1"); //��ɫ�߶�,�ߵ����ֽ���"line1

	viewer->spinOnce(100);
}
