#pragma once
#include "ui_Widget_SearchKNear.h"

/*****************Ѱ�ҵ������ݿռ仮�ֺ�ĳһ���k�������********************/
enum SearchType
{
	VOVEXSEARCH = 1,				
	KNEARESTNODESEARCH,				//	k��������������
	RADIUSSEARCH					//	�뾶����
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


			//	k��������������
signals:
	//void SignalSearchKNear(int searchType,float x,float y,float z,int &k);
private slots:
	void OnPushButton_voxelSearch();
	void OnPushButton_kNearNodes();
	void OnPushButton_NearRadiusSearch();
};


