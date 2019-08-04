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
	char *buffer = nullptr;				//	临时存储文件内容
	int i = 0,pointSize = 0;
	float x, y, z;				//	点云坐标
	string useless;				//	无效数据（文件的头部以及一些其它信息）

	//	读取文件内容
	ReadBuffer(fileName, &buffer);
	stringstream ss(buffer);
	ss.precision(std::numeric_limits<long double>::digits10);		//	设置读取double类型数据的进度，stringstream默认精度为6

	// 去掉首行无用数据
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

	/* 获取文件大小 */
	fseek(pFile, 0, SEEK_END);
	lSize = ftell(pFile);
	rewind(pFile);

	/* 分配内存存储整个文件 */
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