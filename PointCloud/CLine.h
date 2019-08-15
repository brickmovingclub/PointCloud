#pragma once


//	œﬂ¿‡
class CLine
{
public:
	CLine(struct Point pointStart,struct Point pointEnd);
	CLine();
	~CLine();
	bool CLine::operator ==(const CLine& a)const;
private:
	struct Point _pointStart;
	struct Point _pointEnd;

public:
	float LineLength(CLine line);
	float LineLength_Point(Point p1, Point p2);
	Point getPointStart();
	Point getPointEnd();
};

