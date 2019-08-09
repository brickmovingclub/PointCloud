#include <QtWidgets>



#include <pcl/io/pcd_io.h> //PCL��PCD��ʽ�ļ����������ͷ�ļ�
#include <pcl/point_types.h> //PCL�Ը��ָ�ʽ�ĵ��֧��ͷ�ļ�
#include <pcl/visualization/cloud_viewer.h>//���Ʋ鿴����ͷ�ļ�

#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/features/normal_3d.h>
#include <pcl/surface/gp3.h>
#include <pcl/visualization/pcl_visualizer.h>
//#include <boost/thread/thread.hpp>
#include <pcl/surface/poisson.h>
#include <pcl/surface/marching_cubes_hoppe.h>
#include <pcl/surface/marching_cubes_rbf.h>
#include <Eigen/Core>
#include <pcl/common/transforms.h>
#include <pcl/common/common.h>

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