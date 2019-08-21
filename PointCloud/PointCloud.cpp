#include "stdafx.h"

#include "CVector.h"
#include "CLine.h"
#include "CFace.h"

#include "Common.h"



#include "KNearWidget.h"

#include "PointCloud.h"

PointCloud::PointCloud(QWidget *parent)
	: QMainWindow(parent),_cloud(NULL),_viewer(NULL)
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


//读文件--asc/pcd
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
	
	//标记所有的点均为自由点（非排除点与固定点）
	for (auto it = _cloud->begin(); it != _cloud->end(); it++)
	{
		Point p;
		p._x = it->x;
		p._y = it->y;
		p._z = it->z;
		flag.insert(pair<Point, bool>(p, false));
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


//另存为PLY文件
void PointCloud::SaveAsPlY()
{
	fs::path tempFileName = "";
	QString fileName = QFileDialog::getOpenFileName(this, tr("Open"), "ply(.ply)");
	if (fileName.isEmpty())
		return;
	fs::path fileFullName = fileName.toStdWString().c_str();
	pcl::io::savePLYFile(fileFullName.string(), *_cloud);
}


//贪婪三角形重建
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


//泊松三角形重建
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


//打开pcd文件
void PointCloud::open_pcd_file()
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


void PointCloud::OnClear()
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
void PointCloud::ShowLeafNode()
{
	float resolu = 0.01f;
	pcl::octree::OctreePointCloudSearch<pcl::PointXYZ> tree(resolu);
	tree.setInputCloud(_cloud);
	tree.addPointsFromInputCloud();
	std::cout << "叶子节点个数：" << tree.getLeafCount() << std::endl;
	int depth = tree.getTreeDepth();
	std::vector<Eigen::Vector3f> min, max;
	for (auto it = tree.begin(depth); it != tree.end(); it++)
	{
		if (it.isLeafNode())
		{
			Eigen::Vector3f min_pt, max_pt;
			tree.getVoxelBounds(it, min_pt, max_pt);
			min.push_back(min_pt);
			max.push_back(max_pt);
		}
	}

	float r = 0.0f, g = 0.0f, b = 1.0f;
	pcl::visualization::PCLVisualizer viewer;
	int id = 0;
	for (auto it = min.begin(), its = max.begin(); it != min.end(); it++, its++)
	{
		std::cout << "极小值：" << it->x() << "\t" << it->y() << "\t" << it->z() << std::endl;
		std::cout << "极大值：" << its->x() << "\t" << its->y() << "\t" << its->z() << std::endl;
		std::cout << std::endl;
		viewer.addCube(it->x(), its->x(), it->y(), its->y(), it->z(), its->z(), r, g, b, std::to_string(id));
		id++;
	}
	viewer.setBackgroundColor(0.0, 0.0, 0.0);
	while (!viewer.wasStopped())
	{
		viewer.spinOnce();
	}
}


//为整个模型画包围盒
void PointCloud::DrawBoundingBox()
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


void PointCloud::Triangulation()
{
	
	std::list<CFace> ST;					//	三角网格
	std::vector<CLine> activeList; //活动边
	CLine CurrentE; // 当前活动边
	std::vector<CLine> InnerE; //固定边
	std::vector<Point> FreeP; //自由点
	std::vector<Point> ActiveP; //活动点

	int i = 0; int j = 0;
	
	int pointsize = _cloud->points.size();
	pcl::PointXYZ min;//用于存放三个轴的最小值
	pcl::PointXYZ max;//用于存放三个轴的最大值
	pcl::getMinMax3D(*_cloud, min, max);
	double temp = (double)(((max.x - min.x)*(max.y - min.y)*(max.z - min.z) / _cloud->points.size()) );
	double radius = sqrt(temp);	//	搜索半径r

	//	1、求种子三角形
	pcl::PointXYZ pk;		// 默认点云中的第一个点为随机点,第一个点不合适，则自增
	pcl::PointXYZ pi, pj;						//	即将添加的两个点组成三角面片
	do
	{
		pk = _cloud->points[j++];
		std::vector<std::pair< double, pcl::PointXYZ>>	nearPoint;	//	领域点集(double：距离， pcl::PointXYZ：领域点)
		//pcl::PointXYZ origin(0, 0, 0);

		Common::NearRadiusSearch(_cloud, pk, radius , nearPoint);

		std::sort(nearPoint.begin(), nearPoint.end(), [&](const std::pair< double, pcl::PointXYZ> &Pair1, const std::pair< double, pcl::PointXYZ> &Pair2) {return (Pair1.first < Pair2.first ? true : false); });//	按领域点到当前点的距离从小到大排序

		std::vector<std::pair< double, pcl::PointXYZ>>	nearPointBack;	//	不包括pi,pj,pk的领域点备份
		

		bool bfind = false;				//	是否找到种子三角形
		while ((i + 1) < nearPoint.size())
		{
			pi = nearPoint.at(i).second; 
			pj = nearPoint.at(i + 1).second; 

			nearPointBack.clear();
			for (auto iter : nearPoint)
			{
				if ((iter.second.x == pi.x && iter.second.y == pi.y && iter.second.z == pi.z) ||
					(iter.second.x == pj.x && iter.second.y == pj.y && iter.second.z == pj.z) ||
					(iter.second.x == pk.x && iter.second.y == pk.y && iter.second.z == pk.z))
					continue;
				else
					nearPointBack.push_back(iter);
			}

			if (!Common::Condition_a_b(pi, pj, pk, nearPointBack))	//	检测选取得三点是否在同一直线上或经过三点得圆内不包含领域中的其它点
			{
				float dx, dy, dz;
				CVector vector;
				//Common::GetNormal(pi, pj, pk, vector);
				Common::CalNormalVector(pi.x, pi.y, pi.z, pj.x, pj.y, pj.z, pk.x, pk.y, pk.z, dx, dy, dz);
				vector.SetVector(dx, dy, dz);
				float temp = vector.GetX();
				//vector.SetVector(dx, dy, dz);
				if (Common::OnTheSameSide(vector, pk, nearPointBack))
				{
					//	程序走到此步表示当前选取的领域点 pi,pj 可以作为种子三角形的另外两个点
					CLine lineij(pi, pj), linejk(pj, pk), lineki(pk, pi);
					activeList.push_back(lineij);	//	插入种子三角形的活动bian
					activeList.push_back(linejk);
					activeList.push_back(lineki);
					bfind = true;
					break;
				}
				
			}
			i += 2;
		}
		if (bfind)
			break;

	} while (j < _cloud->points.size());
	

	Common::PCLDrawLine(_cloud, _viewer, activeList);

	ui.qvtkWidget->update();

	// 更新自由点
	for (int i = 0; i < _cloud->points.size(); ++i)
	{
		FreeP.push_back(Point(_cloud->points[i].x, _cloud->points[i].y, _cloud->points[i].z));
	}
	
	//	更新活动点
	ActiveP.push_back(Point(pi.x, pi.y, pi.z));
	ActiveP.push_back(Point(pj.x, pj.y, pj.z));
	ActiveP.push_back(Point(pk.x, pk.y, pk.z));

	std::vector<CLine>::iterator itercurretntE;		//	指向当前活动边
	itercurretntE = activeList.begin();

	// 将种子三角形加入到三角网格
	CFace face(pi, pj, pk);
	ST.push_back(face);

	do
	{
		// 初始化当前边
		CurrentE = *itercurretntE;
		Point pointi, pointj, pointk;
		pointi = CurrentE.getPointStart();
		pointj = CurrentE.getPointEnd();
		// 查找精简后选点集
		std::vector<std::pair< double, pcl::PointXYZ>> result;
 		Common::findCandidatePoints(_cloud, pointi, pointj, face.GetOtherPoint(pointi, pointj), flag, activeList, CurrentE, result);
		if (result.size() > 0)
		{
			// 筛选最佳节点
			Point bestP;
			bestP =  Common::FindBestPoint(pointi, pointj, face.GetOtherPoint(pointi, pointj), result, CurrentE, ST, InnerE, flag);
			if (bestP._x == 0 && bestP._y == 0 && bestP._z == 0)
				itercurretntE++;
			else
			{
				//	更新活动边表
				Common::UpdateActiveList(activeList, CurrentE, bestP, InnerE, FreeP, ActiveP, flag, ST);
				itercurretntE = activeList.begin();
			}			
		}
		else
			itercurretntE++;
	}while(!activeList.empty());
	
	
}


