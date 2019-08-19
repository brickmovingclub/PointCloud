#pragma once

//	�ռ�����
class CVector
{
public:
	CVector();
	CVector(const CVector &vector);
	CVector(const float &x, const float &y, const float &z);
	~CVector();

	//	����
	float GetX()const;
	float GetY()const;
	float GetZ()const;

	void SetVector(const float &x, const float &y, const float &z);

	//	����
	CVector operator *(const CVector &vector);				//	�����Ĳ��
	void  operator =(const CVector &vector);
private:
	float _x, _y, _z;
public:
	float vectorMag(CVector &a);
	float vectorInnerProduct(CVector &a, CVector &b);
	CVector GetNormal(Point &p1, Point &p2, Point &p3);
<<<<<<< HEAD
	CVector GetNormal(pcl::PointXYZ &p1, pcl::PointXYZ &p2, pcl::PointXYZ &p3);

=======
	// ��������˺��������ģ��������˵ı�ֵ
	float MultiplicationCross(const CVector &vector);
>>>>>>> origin/dev_hhy
};

