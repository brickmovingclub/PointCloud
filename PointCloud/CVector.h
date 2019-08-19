#pragma once

//	空间向量
class CVector
{
public:
	CVector();
	CVector(const CVector &vector);
	CVector(const float &x, const float &y, const float &z);
	~CVector();

	//	方法
	float GetX()const;
	float GetY()const;
	float GetZ()const;

	void SetVector(const float &x, const float &y, const float &z);

	//	重载
	CVector operator *(const CVector &vector);				//	向量的叉积
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
	// 两向量叉乘后的向量的模与向量点乘的比值
	float MultiplicationCross(const CVector &vector);
>>>>>>> origin/dev_hhy
};

