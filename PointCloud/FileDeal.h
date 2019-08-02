#pragma once
class FileDeal
{
public:
	FileDeal();
	~FileDeal();
	void fileChange();
private:
	void readFile();
	void writeFile();
public:
	int test();
};
