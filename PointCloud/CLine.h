#pragma once


//	œﬂ¿‡
class CLine
{
public:
	CLine(struct Point pointStart,struct Point pointEnd);
	CLine();
	~CLine();
private:
	struct Point _pointStart;
	struct Point _pointEnd;
public:
	float LineLength(CLine line);
	float LineLength_Point(Point p1, Point p2);
};

