#pragma once
#include "ui_Widget_SearchKNear.h"

/*****************寻找点云数据空间划分后某一点的k阶领域点********************/
enum SearchType
{
	VOVEXSEARCH = 1,				
	KNEARESTNODESEARCH,				//	k个最近领域点搜索
	RADIUSSEARCH					//	半径搜索
};
class KNearWidget :
	public QWidget
{
	Q_OBJECT

public:
	KNearWidget(QWidget  *parent = Q_NULLPTR);
	~KNearWidget();
	void SetCloudAndViewer(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, pcl::visualization::PCLVisualizer::Ptr viewer);
private:
	Ui::Widget_SearchKNear ui;
	pcl::PointCloud<pcl::PointXYZ>::Ptr _cloud;
	pcl::visualization::PCLVisualizer::Ptr _viewer;


			//	k个最近领域点搜索
signals:
	//void SignalSearchKNear(int searchType,float x,float y,float z,int &k);
private slots:
	void OnPushButton_voxelSearch();
	void OnPushButton_kNearNodes();
	void OnPushButton_NearRadiusSearch();
};


