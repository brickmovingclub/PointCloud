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
};

