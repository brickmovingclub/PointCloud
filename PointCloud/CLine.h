#pragma once


//	œﬂ¿‡
class CLine
{
public:
	CLine(struct Point pointStart,struct Point pointEnd);
	~CLine();
private:
	struct Point _pointStart;
	struct Point _pointEnd;
};

