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
private:
	vector<string> file;
	int point_count = 0;

};
