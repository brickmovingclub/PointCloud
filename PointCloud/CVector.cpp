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
