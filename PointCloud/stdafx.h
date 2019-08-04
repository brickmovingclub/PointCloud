#include <QtWidgets>



#include <pcl/io/pcd_io.h> //PCL的PCD格式文件的输入输出头文件
#include <pcl/point_types.h> //PCL对各种格式的点的支持头文件
#include <pcl/visualization/cloud_viewer.h>//点云查看窗口头文件


// Visualization Toolkit (VTK)
#include <vtkRenderWindow.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>

namespace fs = std::experimental::filesystem;
using namespace std;