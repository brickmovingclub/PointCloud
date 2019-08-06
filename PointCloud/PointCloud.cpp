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
	

	
}

void PointCloud::ShowLeafNode()
{
	// 何洪玉
}
