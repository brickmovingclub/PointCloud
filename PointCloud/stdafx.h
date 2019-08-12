#include <QtWidgets>

#include <pcl/io/ply_io.h>
#include <pcl/io/pcd_io.h> //PCL的PCD格式文件的输入输出头文件
#include <pcl/point_types.h> //PCL对各种格式的点的支持头文件
#include <pcl/visualization/cloud_viewer.h>//点云查看窗口头文件

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

//画包围盒
#include <pcl/ModelCoefficients.h>
#include <pcl/filters/project_inliers.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/passthrough.h>
#include <pcl/features/normal_3d.h>
#include <pcl/filters/radius_outlier_removal.h>
#include <pcl/segmentation/extract_clusters.h>
#include <Eigen/Core>
#include <pcl/common/transforms.h>
#include <pcl/common/common.h>
#include <pcl/common/time.h>
#include <pcl/common/angles.h>
#include <pcl/registration/transformation_estimation_svd.h>


// Visualization Toolkit (VTK)
#include <vtkRenderWindow.h>
#include <vtkMath.h>
#include <vtkPlane.h>
#include <vtkSmartPointer.h>


#include <limits>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>

namespace fs = std::experimental::filesystem;
using namespace std;


//	空间点结构体
struct Point
{
	float _x;
	float _y;
	float _z;
	Point()
	{
		this->_x = 0.0f;
		this->_y = 0.0f;
		this->_z = 0.0f;
	}
	Point(float x,float y,float z)
	{
		_x = x;
		_y = y;
		_z = z;
	}
	Point(const Point &point)
	{
		this->_x = point._x;
		this->_y = point._y;
		this->_z = point._z;
	}
	~Point()
	{
	}

};

