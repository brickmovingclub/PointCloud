#include "stdafx.h"
#include "PCLViewer.h"


PCLViewer::PCLViewer()
{

	red = 128;
	green = 128;
	blue = 128;
}


PCLViewer::~PCLViewer()
{
}

void PCLViewer::ReadPcdFile(pcl::PointCloud<pcl::PointXYZRGBA>::Ptr cloud, const fs::path &fileName)
{
	char *buffer = nullptr;				//	��ʱ�洢�ļ�����
	int i = 0,pointSize = 0;
	float x, y, z;				//	��������
	string useless;				//	��Ч���ݣ��ļ���ͷ���Լ�һЩ������Ϣ��

	//	��ȡ�ļ�����
	ReadBuffer(fileName, &buffer);
	stringstream ss(buffer);
	ss.precision(std::numeric_limits<long double>::digits10);		//	���ö�ȡdouble�������ݵĽ��ȣ�stringstreamĬ�Ͼ���Ϊ6

	// ȥ��������������
	do {
		ss >> useless;
		if (useless == "DATA")
		{
			getline(ss, useless);
			break;
		}
		else if (useless == "POINTS")
			ss >> pointSize;
		getline(ss, useless);
	} while (1);

	cloud->resize(pointSize);
	// Fill the cloud with some points
	while (!ss.eof() && i < pointSize)
	{
		ss >> x >> y >> z;
		cloud->points[i].x = x;
		cloud->points[i].y = y;
		cloud->points[i].z = z;

		cloud->points[i].r = red;
		cloud->points[i].g = green;
		cloud->points[i].b = blue;

		i++;

	}
	free(buffer);

}


void PCLViewer::ReadBuffer(const fs::path &fileDir, char **buffer)
{
	FILE * pFile;
	long lSize;
	size_t result;

	string fileName = fileDir.string();
	fopen_s(&pFile, fileName.c_str(), "rb");
	if (pFile == NULL)
	{
		fputs("File error", stderr);
		exit(1);
	}

	/* ��ȡ�ļ���С */
	fseek(pFile, 0, SEEK_END);
	lSize = ftell(pFile);
	rewind(pFile);

	/* �����ڴ�洢�����ļ� */
	*buffer = (char*)malloc(sizeof(char)*lSize);
	if (*buffer == NULL)
	{
		return;
	}

	result = fread(*buffer, 1, lSize, pFile);
	if (result != lSize)
	{
		return;
	}
	fclose(pFile);
}

void PCLViewer::AscToPcd(const fs::path &fileName)
{
	if (fileName.extension() != ".asc")
		return;
	vector<string> file;
	int point_count = 0;

	// TODO: �ڴ˴����ʵ�ִ���.
	string s;
	ifstream f;
	f.open(fileName.string(), ios::in);
	while (!f.eof())
	{
		getline(f, s);
		file.push_back(s);
		point_count++;
	}
	if (s == "\n" || s == "")
		point_count--;
	f.close();



	// TODO: �ڴ˴����ʵ�ִ���.
	fs::path  pathFileNamePCD = fileName.filename().string() + ".pcd";
	ofstream file_write(pathFileNamePCD.string(), ios::out);
	file_write << "# .PCD v.5 - Point Cloud Data file format" << "\n";
	file_write << "VERSION .5" << "\n";
	file_write << "FIELDS x y z" << "\n";
	file_write << "SIZE 4 4 4" << "\n";
	file_write << "TYPE F F F" << "\n";
	file_write << "COUNT 1 1 1" << "\n";
	file_write << "WIDTH " << point_count << "\n";
	file_write << "HEIGHT 1" << "\n";
	file_write << "POINTS " << point_count << "\n";
	file_write << "DATA ascii" << "\n";
	for (auto it = file.begin(); it != file.end(); it++)
		file_write << *it << '\n';
	//file_write << it << endl;

	file_write.close();
}