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