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


void Common::PrintString(const string &str)
{

}


//<<<<<<< HEAD
//	�����еĵ���������Ƭͬ��
bool Common::OnTheSameSide(const CVector &normal, const Point &origin, const std::list<Point> &nearPoints)
{
	int i = 0;
	double temp = 0.0f;
	vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
	plane->SetOrigin(origin._x, origin._y, origin._z);
	plane->SetNormal(normal.GetX(), normal.GetY(), normal.GetZ());

	for (auto iter : nearPoints)
	{
		temp = plane->EvaluateFunction(iter._x, iter._y, iter._z);
		if ((temp  - 0.0f)> 0)
			i++;
		else if ((temp - 0.0f) < 0)
			i--;
	}
	return (i == nearPoints.size() ? true :false);
}


//	���ݿ���ķ�����������������ķ���
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


//	��ռ�ƽ��ķ�����
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
bool Common::Condition_a_b(pcl::PointXYZ pi, pcl::PointXYZ pj, pcl::PointXYZ pk, std::vector<pcl::PointXYZ> near_pi)
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
		distance = sqrt(pow(centerpoint.x - it->x, 2) + pow(centerpoint.y - it->y, 2) + pow(centerpoint.z - it->z, 2));
		if (distance < dR)
			return false;
	}
	return true;
//>>>>>>> origin/dev_hhy
}


//ѡ���ѡ�㼯
std::vector<int> Common::findCandidatePoints(pcl::PointCloud<pcl::PointXYZ>::Ptr _cloud, Point pi, Point pj, Point pk, std::map<Point, bool> flag, std::vector<CLine> ActiveE, CLine CurrentE)
{
	//�õ����߱߳����ж�����������
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
	//������뾶r
	float r;
	if (c * c > a * a + b * b) //�۽�������
		r = c;
	else
		r = (a + b + c) / 3;

	Point pm;
	pm._x = (pi._x + pj._x) / 2;
	pm._y = (pi._y + pj._y) / 2;
	pm._z = (pi._z + pj._z) / 2;

	//�����pm��r��Χ�ڵ�����㼯
	std::vector<int> near_pm; 


	//��������㼯
	int size = near_pm.size();
	for (int i = 0; i < size; i++)
	{
		Point p;
		p._x = _cloud->points[near_pm[i]].x;
		p._y = _cloud->points[near_pm[i]].y;
		p._z = _cloud->points[near_pm[i]].z;
		//ɾ���̶������ų���
		if (flag.count(p) > 0 && flag[p])
			near_pm.erase(near_pm.begin() + i);
	}
	if(near_pm.empty()) //��pi-pj�Ǳ߽��
		return near_pm;

	//�ҵ�pi,pj���ڵĻ��
	std::vector<CLine>::iterator it = std::find(ActiveE.begin(), ActiveE.end(), CurrentE);
	it--;
	Point pa = it->getPointStart();
	it++; it++;
	Point pb = it->getPointEnd();

	//�߽Ƕ�Լ����
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
	if (near_pm.empty()) //��pi-pj�Ǳ߽��
		return near_pm;

	//��Ƕ�Լ����
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
	return near_pm; //��Ϊ�գ���pi-pj�Ǳ߽��
}


//ѡ����ѵ�
Point Common::FindBestPoint(pcl::PointCloud<pcl::PointXYZ>::Ptr _cloud, Point pi, Point pj, Point pk, std::vector<int> near_pm, CLine CurrentE, std::list<CFace> ST, std::vector<CLine> InnerE, std::map<Point, bool> flag)
{
	// TODO: �ڴ˴����ʵ�ִ���.
	CVector vec; //��ǰ������������εķ�����
	vec = vec.GetNormal(pi, pj, pk);
	const float PI = 3.1415926;
	float angle_min, angle_max; //��ѡ�����εĽǶȵ����ֵ����Сֵ
	CLine l;
	float a = l.LineLength_Point(pi, pj);
	float b = l.LineLength_Point(pj, pk);
	float c = l.LineLength_Point(pk, pi);
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

	std::map<float, int> joinCosts;//��ӵĴ���
	float angle_cur;//��ǰ��ѡ����Ƭ������ڽ�,��߶Դ��
	//Ϊ��ѡ����Ӵ���
	for (auto it = near_pm.begin(); it != near_pm.end(); it++)
	{
		Point candidateP;
		candidateP._x = _cloud->points[*it].x;
		candidateP._y = _cloud->points[*it].y;
		candidateP._z = _cloud->points[*it].z;
		float joinCost; //�ڵ�candidateP�Ĵ���
		float cost_angle1, cost_angle2;
		//��ѡ�����εķ�����
		CVector newVec = vec.GetNormal(pi, candidateP, pj);
		float angle_cos; 
		angle_cos = vec.vectorInnerProduct(vec, newVec) / (vec.vectorMag(vec) * vec.vectorMag(newVec));
		cost_angle1 = sqrt(1 - angle_cos * angle_cos);

		a = l.LineLength_Point(pi, candidateP);
		b = l.LineLength_Point(candidateP, pj);
		c = l.LineLength_Point(pj, pi);
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
		joinCosts.insert(pair<float, int>(joinCost, *it));
	}
	
	//��ѡ����Ƭ�������
	Point bestP;
	for (auto it = joinCosts.begin(); it != joinCosts.end(); it++)
	{
		Point pc;
		pc._x = _cloud->points[it->second].x;
		pc._y = _cloud->points[it->second].y;
		pc._z = _cloud->points[it->second].z;	

		//�������е�
		CFace curFace(pi, pj, pc);
		if (findNearFace_Point(curFace, pc, ST))
			continue;

		//�������б�
		CLine p_ic(pi, pc);
		vector<CLine>::iterator its = find(InnerE.begin(), InnerE.end(), p_ic);
		if (its != InnerE.end())
			continue;
		CLine p_cj(pc, pj);
		its = find(InnerE.begin(), InnerE.end(), p_cj);
		if (its != InnerE.end())
			continue;

		//����Խ�������
		if (IntersectTriangle(pi, pj, pc, ST))
			continue;
		else
		{
			bestP = pc;
			break;
		}			
	}

	//����Ƿ��������
	CVector vec_new;
	vec_new = vec_new.GetNormal(pi, pj, bestP);
	a = vec_new.GetX();
	b = vec_new.GetY();
	c = vec_new.GetZ();
	float d = -a * pi._x - b * pi._y - c * pi._z;

	for (auto it = near_pm.begin(); it != near_pm.end(); it++)
	{
		if ((_cloud->points[*it].x == bestP._x) 
		&& (_cloud->points[*it].y == bestP._y)
		&& (_cloud->points[*it].z == bestP._z))
			continue;
		Point  spatialpoint;;
		spatialpoint._x = _cloud->points[*it].x;
		spatialpoint._y = _cloud->points[*it].y;
		spatialpoint._z = _cloud->points[*it].z;
		Point subpoint; //ͶӰ��
		subpoint._x = ((b * b + c * c) * spatialpoint._x - a * (b * spatialpoint._y + c * spatialpoint._z + d)) / (a * a + b * b + c * c);
		subpoint._y = b / a * (subpoint._x - spatialpoint._x) + spatialpoint._y;
		subpoint._z = c / a * (subpoint._x - spatialpoint._x) + spatialpoint._z;

		//�ж�ͶӰ���Ƿ�����������
		if (TriangleIncludeSubpoint(pi, pj, bestP, subpoint))
			flag[spatialpoint] = true; //��ǰ������Ϊ���ɵ�/�ų���
	}
	CFace f(pi, pj, bestP);
	ST.push_back(f); //���������Ƭ
	return bestP;
}


// Ѱ�ҵ���ڽ�������,���ж��Ƿ񹲱�
bool Common::findNearFace_Point(CFace curFace, Point pc, std::list<CFace> ST)
{
	// TODO: �ڴ˴����ʵ�ִ���.
	std::set<Point> points;
	points.insert(curFace.GetPoint1());
	points.insert(curFace.GetPoint2());
	points.insert(curFace.GetPoint3());

	for (auto it = ST.begin(); it != ST.end(); it++)
	{
		//�ǵ�pc���ڽ�������		
		if (((it->GetPoint1()._x == pc._x) && (it->GetPoint1()._y == pc._y) && (it->GetPoint1()._z == pc._z)) 
			|| ((it->GetPoint2()._x == pc._x) && (it->GetPoint2()._y == pc._y) && (it->GetPoint2()._z == pc._z))
			|| ((it->GetPoint3()._x == pc._x) && (it->GetPoint3()._y == pc._y) && (it->GetPoint3()._z == pc._z)))
		{
			points.insert(it->GetPoint1());
			points.insert(it->GetPoint2());
			points.insert(it->GetPoint3());
			if (points.size() <= 4)
				return false;
		}			
	}
	return true;
}


// ������pi_pj_pc��pi_pj���ɵıߵ��ڽ������Σ��������������Ƿ��ཻ
bool Common::IntersectTriangle(Point pi, Point pj, Point pc, std::list<CFace> ST)
{
	// TODO: �ڴ˴����ʵ�ִ���.
	for (auto it = ST.begin(); it != ST.end(); it++)
	{
		if ((it->GetPoint1()._x == pi._x) && (it->GetPoint1()._y == pi._y) && (it->GetPoint1()._z == pi._z))
		{
			if ((it->GetPoint2()._x == pj._x) && (it->GetPoint2()._y == pj._y) && (it->GetPoint2()._z == pj._z))
			{
				Point pk = it->GetPoint3();
				if (IntersectionLine(pi, pc, pj, pk) 
				 || IntersectionLine(pj, pc, pi, pk))
					return true;
			}
			else if ((it->GetPoint3()._x == pj._x) && (it->GetPoint3()._y == pj._y) && (it->GetPoint3()._z == pj._z))
			{
				Point pk = it->GetPoint2();
				if (IntersectionLine(pi, pc, pj, pk)
					|| IntersectionLine(pj, pc, pi, pk))
					return true;
			}
		}
		else if ((it->GetPoint2()._x == pi._x) && (it->GetPoint2()._y == pi._y) && (it->GetPoint2()._z == pi._z))
		{
			if ((it->GetPoint1()._x == pj._x) && (it->GetPoint1()._y == pj._y) && (it->GetPoint1()._z == pj._z))
			{
				Point pk = it->GetPoint3();
				if (IntersectionLine(pi, pc, pj, pk)
					|| IntersectionLine(pj, pc, pi, pk))
					return true;
			}
			else if ((it->GetPoint3()._x == pj._x) && (it->GetPoint3()._y == pj._y) && (it->GetPoint3()._z == pj._z))
			{
				Point pk = it->GetPoint1();
				if (IntersectionLine(pi, pc, pj, pk)
					|| IntersectionLine(pj, pc, pi, pk))
					return true;
			}
		}
		else if ((it->GetPoint3()._x == pi._x) && (it->GetPoint3()._y == pi._y) && (it->GetPoint3()._z == pi._z))
		{
			if ((it->GetPoint2()._x == pj._x) && (it->GetPoint2()._y == pj._y) && (it->GetPoint2()._z == pj._z))
			{
				Point pk = it->GetPoint1();
				if (IntersectionLine(pi, pc, pj, pk)
					|| IntersectionLine(pj, pc, pi, pk))
					return true;
			}
			else if ((it->GetPoint1()._x == pj._x) && (it->GetPoint1()._y == pj._y) && (it->GetPoint1()._z == pj._z))
			{
				Point pk = it->GetPoint2();
				if (IntersectionLine(pi, pc, pj, pk)
					|| IntersectionLine(pj, pc, pi, pk))
					return true;
			}
		}
	}
	return false;
}


// �ж������߶��Ƿ��ཻ
bool Common::IntersectionLine(Point pi, Point pj, Point pc, Point pk)
{
	// TODO: �ڴ˴����ʵ�ִ���.
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


// �ж�ǰ�������㹹�ɵ��������Ƿ�������ĸ���
bool Common::TriangleIncludeSubpoint(Point pi, Point pj, Point bestP, Point subpoint)
{
	// TODO: �ڴ˴����ʵ�ִ���.
	float area = TriangleArea(pi, pj, bestP);
	float area_new = TriangleArea(pi, pj, subpoint) + TriangleArea(pi, bestP, subpoint) + TriangleArea(pj, bestP, subpoint);
	if (area_new > area)
		return false;
	else
		return true;
}


// �õ�������ABC�����
float Common::TriangleArea(Point A, Point B, Point C)
{
	// TODO: �ڴ˴����ʵ�ִ���.
	float s = ((B._y - A._y) * (C._z - A._z) + (B._z - A._z) * (C._x - A._x) + (B._x - A._x) * (C._y - A._y)) -
		((C._y - A._y) * (B._z - A._z) + (C._z - A._z) * (B._x - A._x) + (C._x - A._x) * (B._y - A._y));
	return 0.5 * s;
}


// ���»����
void Common::UpdateActiveList(std::vector<CLine> ActiveE, CLine CurrentE, Point bestP, std::vector<CLine> InnerE, std::vector<Point> FreeP, std::vector<Point> ActiveP, std::map<Point, bool> flag)
{
	// �ж���ѵ���ӵ�λ��
	int type = BestPositionType(ActiveE, CurrentE, bestP, FreeP);
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
		UpdateMode3(ActiveE, CurrentE, bestP, InnerE, FreeP, ActiveP, flag);
		break;
	default:
		break;
	}
}


// �ж���ѵ���ӵ�λ������
/*
0 -- ��ѵ������ɵ�
1 -- ��ѵ�λ�ڻ������Ϊ��ǰ���ǰ���ڱߵĶ˵�
2 -- ��ѵ�λ�ڻ������Ϊ��ǰ��ߺ����ڱߵĶ˵�
3 -- ��ѵ�λ�ڻ�������뵱ǰ���û�����ڹ�ϵ
*/
int Common::BestPositionType(std::vector<CLine> ActiveE, CLine CurrentE, Point bestP, std::vector<Point> FreeP)
{
	// TODO: �ڴ˴����ʵ�ִ���.
	std::vector<Point>::iterator it = find(FreeP.begin(), FreeP.end(), bestP);
	if (it != FreeP.end()) //��ѵ������ɵ�
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
		front = its--;
		behind = its++;
	}
	if ((front->getPointStart()._x == bestP._x)
		&& (front->getPointStart()._y == bestP._y)
		&& (front->getPointStart()._z == bestP._z))
		return 1;
	if ((behind->getPointEnd()._x == bestP._x)
		&& (behind->getPointEnd()._y == bestP._y)
		&& (behind->getPointEnd()._z == bestP._z))
		return 2;
	return 3;
}


/*************************���»����*************************************/
// ��ѵ������ɵ�
void Common::UpdateMode(std::vector<CLine> ActiveE, CLine CurrentE, Point bestP, std::vector<CLine> InnerE, std::vector<Point> FreeP, std::vector<Point> ActiveP)
{
	// TODO: �ڴ˴����ʵ�ִ���.
	InnerE.push_back(CurrentE);
	std::vector<Point>::iterator its = find(FreeP.begin(), FreeP.end(), bestP);
	FreeP.erase(its);
	ActiveP.push_back(bestP);

	std::vector<CLine>::iterator it = find(ActiveE.begin(), ActiveE.end(), CurrentE);
	CLine line1(CurrentE.getPointStart(), bestP);
	CLine line2(bestP, CurrentE.getPointEnd());
	ActiveE.insert(it, line1);
	ActiveE.insert(it, line2);
	ActiveE.erase(it);	
}


// ��ѵ�λ�ڻ������Ϊ��ǰ���ǰ���ڱߵĶ˵�
void Common::UpdateMode1(std::vector<CLine> ActiveE, CLine CurrentE, Point bestP, std::vector<CLine> InnerE, std::vector<Point> ActiveP, std::map<Point, bool> flag)
{
	// TODO: �ڴ˴����ʵ�ִ���.
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
	
	flag[CurrentE.getPointStart()] = true;
}


// ��ѵ�λ�ڻ������Ϊ��ǰ��ߺ����ڱߵĶ˵�
void Common::UpdateMode2(std::vector<CLine> ActiveE, CLine CurrentE, Point bestP, std::vector<CLine> InnerE, std::vector<Point> ActiveP, std::map<Point, bool> flag)
{
	// TODO: �ڴ˴����ʵ�ִ���.
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

	flag[CurrentE.getPointEnd()] = true;
}


//��ѵ�λ�ڻ�������뵱ǰ���û�����ڹ�ϵ
void Common::UpdateMode3(std::vector<CLine> ActiveE, CLine CurrentE, Point bestP, std::vector<CLine> InnerE, std::vector<Point> FreeP, std::vector<Point> ActiveP, std::map<Point, bool> flag)
{
	// TODO: �ڴ˴����ʵ�ִ���.

}
