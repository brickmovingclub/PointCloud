#include "stdafx.h"

#include "FileDeal.h"




using namespace std;


FileDeal::FileDeal()
{

}


FileDeal::~FileDeal()
{
}


//�ļ�����ת��
void FileDeal::fileChange()
{
	// TODO: �ڴ˴����ʵ�ִ���.
	readFile();
	writeFile();
}

//��asc�ļ�
void FileDeal::readFile()
{
	// TODO: �ڴ˴����ʵ�ִ���.
	string s;
	ifstream f;
	f.open("D:\\code\\c++\\PointCloud\\bunny.asc", ios::in);
	while (!f.eof())
	{
		getline(f, s);
		file.push_back(s);
		point_count++;
	}
	if (s == "\n" || s == "")
		point_count--;
	f.close();
}


//дpcd�ļ�
void FileDeal::writeFile()
{
	// TODO: �ڴ˴����ʵ�ִ���.
	ofstream file_write("D:\\code\\c++\\PointCloud\\bunny.pcd", ios::out);
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


int FileDeal::test()
{
	// TODO: �ڴ˴����ʵ�ִ���.
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>); // �������ƣ�ָ�룩

	/*C:\\Users\\DR\\Desktop\\pcdFile\\bunny.pcd*/
	/*D:\\code\\c++\\PointCloud\\PointCloudDeal\\bunny.pcd*/
	if (pcl::io::loadPCDFile<pcl::PointXYZ>("D:\\code\\c++\\PointCloud\\bunny.pcd", *cloud) == -1) //* ����PCD��ʽ���ļ�������ļ������ڣ�����-1
	{
		PCL_ERROR("Couldn't read file bunny.pcd \n"); //�ļ�������ʱ�����ش�����ֹ����
		return (-1);
	}
	pcl::visualization::CloudViewer viewer("Simple Cloud Viewer");//ֱ�Ӵ���һ����ʾ����
	viewer.showCloud(cloud);//�����������ʾ����
	while (!viewer.wasStopped())
	{

	}
	return (0);
}
