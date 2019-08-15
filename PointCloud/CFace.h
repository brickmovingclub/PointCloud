#pragma once

//	中頭
class CFace
{
public:
	CFace(const struct Point &point1, const struct Point &point2, const struct Point &point3);
	~CFace();

	Point GetPoint1()const;
	Point GetPoint2()const;
	Point GetPoint3()const;

	void  CFace::operator =(const CFace &p);
private:
	struct Point _point1;
	struct Point _point2;
	struct Point _point3;
};

// 利鯉爆中 std::list<CFace> ST