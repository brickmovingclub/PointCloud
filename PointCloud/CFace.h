#pragma once

//	��Ƭ
class CFace
{
public:
	CFace(const struct Point &point1, const struct Point &point2, const struct Point &point3);
	~CFace();
private:
	struct Point _point1;
	struct Point _point2;
	struct Point _point3;
};

// �������� std::list<CFace> ST