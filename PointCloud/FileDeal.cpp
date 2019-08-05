#include "stdafx.h"

#include "FileDeal.h"




using namespace std;


FileDeal::FileDeal()
{

}


FileDeal::~FileDeal()
{
}


//文件类型转换
void FileDeal::fileChange()
{
	// TODO: 在此处添加实现代码.
	readFile();
	writeFile();
}

//读asc文件
void FileDeal::readFile()
{
	// TODO: 在此处添加实现代码.
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


//写pcd文件
void FileDeal::writeFile()
{
	// TODO: 在此处添加实现代码.
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
	// TODO: 在此处添加实现代码.
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>); // 创建点云（指针）

	/*C:\\Users\\DR\\Desktop\\pcdFile\\bunny.pcd*/
	/*D:\\code\\c++\\PointCloud\\PointCloudDeal\\bunny.pcd*/
	if (pcl::io::loadPCDFile<pcl::PointXYZ>("D:\\code\\c++\\PointCloud\\bunny.pcd", *cloud) == -1) //* 读入PCD格式的文件，如果文件不存在，返回-1
	{
		PCL_ERROR("Couldn't read file bunny.pcd \n"); //文件不存在时，返回错误，终止程序。
		return (-1);
	}
	pcl::visualization::CloudViewer viewer("Simple Cloud Viewer");//直接创造一个显示窗口
	viewer.showCloud(cloud);//再这个窗口显示点云
	while (!viewer.wasStopped())
	{

	}
	return (0);
}
