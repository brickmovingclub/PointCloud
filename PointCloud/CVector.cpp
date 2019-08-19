#include "stdafx.h"
#include "CVector.h"


CVector::CVector() :_x(0), _y(0), _z(0)
{
}

CVector::CVector(const CVector &vector)
{
	this->_x = vector.GetX();
	this->_y = vector.GetY();
	this->_z = vector.GetZ();
}

CVector::CVector(const float &x, const float &y, const float &z) :_x(x), _y(y), _z(z)
{

}


CVector::~CVector()
{
}

float CVector::GetX()const
{
	return _x;
}
float CVector::GetY()const
{
	return _y;
}
float CVector::GetZ()const
{
	return _z;
}
void CVector::SetVector(const float &x, const float &y, const float &z)
{
	this->_x = x;
	this->_y = y;
	this->_z = z;
}
CVector CVector::operator *(const CVector &vector)
{
	//(a1,a2,a3)x(b1,b2,b3)=(a2b3-a3b2,a3b1-a1b3,a1b2-a2b1)
	float normalx, normaly, normalz;
	normalx = (this->GetY()*vector.GetZ() - this->GetZ()*vector.GetY());
	normaly = (this->GetZ()*vector.GetX() - this->GetX() *vector.GetZ());
	normalz = (this->GetX()*vector.GetY() - this->GetY() *vector.GetX());
	float  len = (float)sqrt(normalx *normalx + normaly * normaly + normalz * normalz);
	return CVector(normalx / len, normaly / len, normalz / len);
}
void  CVector::operator =(const CVector &vector)
{
	this->_x = vector._x;
	this->_y = vector._y;
	this->_z = vector._z;
}


//计算向量的模
float CVector::vectorMag(CVector &a)
{
	// TODO: 在此处添加实现代码.
	return sqrt(a._x * a._x + a._y * a._y + a._z * a._z);
}


//计算向量点乘
float CVector::vectorInnerProduct(CVector &a, CVector &b)
{
	// TODO: 在此处添加实现代码.
	return a._x * b._x + a._y * b._y + a._z * b._z;
}


CVector CVector::GetNormal(Point &p1, Point &p2, Point &p3)
{
	// TODO: 在此处添加实现代码.
	CVector v1(p2._x - p1._x, p2._y - p1._y, p2._z - p1._z);
	CVector v2(p3._x - p2._x, p3._y - p2._y, p3._z - p2._z);
	CVector v3(p1._x - p3._x, p1._y - p3._y, p1._z - p3._z);

	float na = (v2._y - v1._y)*(v3._z - v1._z) - (v2._z - v1._z)*(v3._y - v1._y);
	float nb = (v2._z - v1._z)*(v3._x - v1._x) - (v2._x - v1._x)*(v3._z - v1._z);
	float nc = (v2._x - v1._x)*(v3._y - v1._y) - (v2._y - v1._y)*(v3._x - v1._x);

	return CVector(na, nb, nc);
}

CVector CVector::GetNormal(pcl::PointXYZ &p1, pcl::PointXYZ &p2, pcl::PointXYZ &p3)
{
	CVector v1(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);
	CVector v2(p3.x - p2.x, p3.y - p2.y, p3.z - p2.z);
	CVector v3(p1.x - p3.x, p1.y - p3.y, p1.z - p3.z);

	float na = (v2._y - v1._y)*(v3._z - v1._z) - (v2._z - v1._z)*(v3._y - v1._y);
	float nb = (v2._z - v1._z)*(v3._x - v1._x) - (v2._x - v1._x)*(v3._z - v1._z);
	float nc = (v2._x - v1._x)*(v3._y - v1._y) - (v2._y - v1._y)*(v3._x - v1._x);

	return CVector(na, nb, nc);
}
