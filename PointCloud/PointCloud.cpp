#include "stdafx.h"

#include "FileDeal.h"
#include "KNearWidget.h"

#include "PointCloud.h"

PointCloud::PointCloud(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);



	this->setWindowTitle("PCL viewer");

	// Setup the cloud pointer
	_cloud.reset(new pcl::PointCloud<pcl::PointXYZ>);
	
	// Set up the QVTK window
	_viewer.reset(new pcl::visualization::PCLVisualizer("viewer", false));
	ui.qvtkWidget->SetRenderWindow(_viewer->getRenderWindow());
	_viewer->setupInteractor(ui.qvtkWidget->GetInteractor(), ui.qvtkWidget->GetRenderWindow());
	ui.qvtkWidget->update();
	this->setCentralWidget(ui.qvtkWidget);
	this->setWindowTitle("PCL viewer");
}


void PointCloud::OnReadFile()
{
	//_cloud->clear();
	//_cloud.reset(new pcl::PointCloud<pcl::PointXYZ>);
	fs::path tempFileName = "";
	QString fileName = QFileDialog::getOpenFileName(this, tr("Open"), ".");
	if (fileName.isEmpty())
		return;
	fs::path fileFullName = fileName.toStdWString().c_str();
	if (fileFullName.extension() == ".pcd")
	{
		_pclViewer.ReadPcdFile(_cloud, fileFullName);
		tempFileName = fileFullName;
	}
	else if (fileFullName.extension() == ".asc")
	{
		tempFileName = fileFullName.filename().string() + ".pcd";

		_pclViewer.AscToPcd(fileFullName);
		_pclViewer.ReadPcdFile(_cloud, tempFileName);
	}


	pcl::PCLPointCloud2 cloud_blob;
	pcl::io::loadPCDFile(tempFileName.string(), cloud_blob);
	pcl::fromPCLPointCloud2(cloud_blob, *_cloud);
	_viewer->addPointCloud<pcl::PointXYZ>(_cloud, "cloud");

	ui.qvtkWidget->update();
//	_viewer->addPointCloud(_cloud, "cloud");
//	_viewer->resetCamera();
	//ui.qvtkWidget->update();
}

void PointCloud::SaveAsPlY()
{
	fs::path tempFileName = "";
	QString fileName = QFileDialog::getOpenFileName(this, tr("Open"), "ply(.ply)");
	if (fileName.isEmpty())
		return;
	fs::path fileFullName = fileName.toStdWString().c_str();
	pcl::io::savePLYFile(fileFullName.string(), *_cloud);
}


void PointCloud::greedyTriangulation_reconstruct()
{

	//cloud_with_noramls=cloud+normals
	pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals = getPointNormal();

	//create search tree
	pcl::search::KdTree<pcl::PointNormal>::Ptr tree2(new pcl::search::KdTree<pcl::PointNormal>);
	tree2->setInputCloud(cloud_with_normals);

	//initialize objects
	pcl::GreedyProjectionTriangulation<pcl::PointNormal> gp3;
	pcl::PolygonMesh triangles;

	//set the maximum distance between connected points (maximum edge length);
	gp3.setSearchRadius(0.025);

	//set typical values for the parameters
	gp3.setMu(2.5);
	gp3.setMaximumNearestNeighbors(100);
	gp3.setMaximumSurfaceAngle(M_PI / 4);//45degrees
	gp3.setMinimumAngle(M_PI / 18);
	gp3.setMaximumAngle(2 * M_PI / 3);//120degrees
	gp3.setNormalConsistency(false);

	//get result
	gp3.setInputCloud(cloud_with_normals);
	gp3.setSearchMethod(tree2);
	gp3.reconstruct(triangles);

	//std::count<<triangles
	//additional vertex information
	//std::vector<int> parts = gp3.getPartIDs();
	//std::vector<int>states = gp3.getPointStates();
	_viewer->removePointCloud("cloud");
	_viewer->removePolygonMesh("my");

	_viewer->addPolygonMesh(triangles, "my");
	ui.qvtkWidget->update();

	std::cout << "三角面片划分完成" << std::endl;

}
void PointCloud::poisson_reconstruct()
{


	//cloud_with_noramls=cloud+normals
	pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals = getPointNormal();

	//create search tree
	pcl::search::KdTree<pcl::PointNormal>::Ptr tree2(new pcl::search::KdTree<pcl::PointNormal>);
	tree2->setInputCloud(cloud_with_normals);

	pcl::Poisson<pcl::PointNormal> pn;
	pcl::PolygonMesh triangles;

	pn.setConfidence(false);
	pn.setIsoDivide(7);


	pn.setInputCloud(cloud_with_normals);
	pn.setSearchMethod(tree2);
	pn.reconstruct(triangles);

	_viewer->removePointCloud("cloud");
	_viewer->removePolygonMesh("my");
	_viewer->addPolygonMesh(triangles, "my");
	ui.qvtkWidget->update();
}


pcl::PointCloud<pcl::PointNormal>::Ptr PointCloud::getPointNormal()
{
	//Normal estimation
	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal>n;
	pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
	tree->setInputCloud(_cloud);
	n.setInputCloud(_cloud);
	n.setSearchMethod(tree);
	n.setKSearch(20);
	n.compute(*normals);

	//normals should not contain the point normals+surface curvatures
	//concatenate the xyz and normal fields
	pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals(new pcl::PointCloud<pcl::PointNormal>);
	pcl::concatenateFields(*_cloud, *normals, *cloud_with_normals);
	return cloud_with_normals;
}

void  PointCloud::open_pcd_file()
{
	QString filter;
	filter = "PCD file(*.pcd)";

	QDir dir;
	QString fileName = QFileDialog::getOpenFileName(this, QString(tr("open PCD file"))
		, dir.absolutePath(), filter);
	if (fileName.isEmpty() == true)
	{

		std::cout << "empty pcd files";
		return;
	}


	//支持带中文路径的读取
	QByteArray ba = fileName.toLocal8Bit();
	const char * fileName_str = ba.data();

	//_cloud.reset(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PCLPointCloud2 cloud_blob;
	pcl::io::loadPCDFile(fileName_str, cloud_blob);
	pcl::fromPCLPointCloud2(cloud_blob, *_cloud);
	_viewer->addPointCloud<pcl::PointXYZ>(_cloud, "cloud");

	ui.qvtkWidget->update();
}

void	PointCloud::OnClear()
{
	_cloud->clear();
	//_cloud.reset(new pcl::PointCloud<pcl::PointXYZ>);
	_viewer->addPointCloud<pcl::PointXYZ>(_cloud, "cloud");
	ui.qvtkWidget->update();

}

void PointCloud::OnActionSearchKNear()
{
	KNearWidget *widget = new KNearWidget();
	widget->show();
	//connect(widget, SIGNAL(SignalSearchKNear(float x, float y, float z, int &k)), this, SLOT(SearchKNear(int &k)));
}

void PointCloud::SearchKNear(float x, float y, float z, int &k)
{
	// Neighbors within voxel search
	KNearWidget *widget = new KNearWidget();

	widget->show();
	
}


//找到八叉树中的叶子节点并显示
//画包围盒
void PointCloud::ShowLeafNode()
{
	//获取点云质心
	Eigen::Vector4f pcaCentroid;
	pcl::compute3DCentroid(*_cloud, pcaCentroid); 
	//计算协方差，获取协方差矩阵
	Eigen::Matrix3f covariance;
	pcl::computeCovarianceMatrixNormalized(*_cloud, pcaCentroid, covariance); 
	//计算特征向量、特征值
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigen_solver(covariance, Eigen::ComputeEigenvectors);
	Eigen::Matrix3f eigenVectorsPCA = eigen_solver.eigenvectors();
	Eigen::Vector3f eigenValuesPCA = eigen_solver.eigenvalues();
	//校正主方向间垂直，特征向量即为主方向
	eigenVectorsPCA.col(2) = eigenVectorsPCA.col(0).cross(eigenVectorsPCA.col(1)); 
	eigenVectorsPCA.col(0) = eigenVectorsPCA.col(1).cross(eigenVectorsPCA.col(2));
	eigenVectorsPCA.col(1) = eigenVectorsPCA.col(2).cross(eigenVectorsPCA.col(0));

	Eigen::Matrix4f tm = Eigen::Matrix4f::Identity(); //变换矩阵
	Eigen::Matrix4f tm_inv = Eigen::Matrix4f::Identity(); //变换逆矩阵
	tm.block<3, 3>(0, 0) = eigenVectorsPCA.transpose();   
	tm.block<3, 1>(0, 3) = -1.0f * (eigenVectorsPCA.transpose()) *(pcaCentroid.head<3>());
	tm_inv = tm.inverse();

	pcl::PointCloud<pcl::PointXYZ>::Ptr transformedCloud(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::transformPointCloud(*_cloud, *transformedCloud, tm);

	pcl::PointXYZ min_p, max_p;   //点云的最大值与最小值点
	Eigen::Vector3f c1, c;
	pcl::getMinMax3D(*transformedCloud, min_p, max_p);
	c1 = 0.5f*(min_p.getVector3fMap() + max_p.getVector3fMap());

	Eigen::Affine3f tm_inv_aff(tm_inv);
	pcl::transformPoint(c1, c, tm_inv_aff);

	Eigen::Vector3f whd, whd1;
	whd1 = max_p.getVector3fMap() - min_p.getVector3fMap();
	whd = whd1;
	float sc1 = (whd1(0) + whd1(1) + whd1(2)) / 3;  //点云平均尺度，用于设置主方向箭头大小

	const Eigen::Quaternionf bboxQ1(Eigen::Quaternionf::Identity());
	const Eigen::Vector3f    bboxT1(c1);

	const Eigen::Quaternionf bboxQ(tm_inv.block<3, 3>(0, 0));
	const Eigen::Vector3f    bboxT(c);


	//变换到原点的点云主方向
	pcl::PointXYZ op;
	op.x = 0.0;
	op.y = 0.0;
	op.z = 0.0;
	Eigen::Vector3f px, py, pz;
	Eigen::Affine3f tm_aff(tm);
	pcl::transformVector(eigenVectorsPCA.col(0), px, tm_aff);
	pcl::transformVector(eigenVectorsPCA.col(1), py, tm_aff);
	pcl::transformVector(eigenVectorsPCA.col(2), pz, tm_aff);
	pcl::PointXYZ pcaX;
	pcaX.x = sc1 * px(0);
	pcaX.y = sc1 * px(1);
	pcaX.z = sc1 * px(2);
	pcl::PointXYZ pcaY;
	pcaY.x = sc1 * py(0);
	pcaY.y = sc1 * py(1);
	pcaY.z = sc1 * py(2);
	pcl::PointXYZ pcaZ;
	pcaZ.x = sc1 * pz(0);
	pcaZ.y = sc1 * pz(1);
	pcaZ.z = sc1 * pz(2);

	//visualization
	pcl::visualization::PCLVisualizer viewer;

	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> color_handler(_cloud, 255, 255, 0); //输入的初始点云相关
	viewer.addPointCloud(_cloud, color_handler, "cloud");
	viewer.addCube(bboxT, bboxQ, whd(0), whd(1), whd(2), "bbox");
	viewer.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_REPRESENTATION, pcl::visualization::PCL_VISUALIZER_REPRESENTATION_WIREFRAME, "bbox");
	viewer.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0.0, 0.0, "bbox");

	viewer.addCoordinateSystem(0.5f*sc1);
	viewer.setBackgroundColor(0.0, 0.0, 0.0);
	while (!viewer.wasStopped())
	{
		viewer.spinOnce();
	}
}
