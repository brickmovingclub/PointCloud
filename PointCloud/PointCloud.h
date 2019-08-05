#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_PointCloud.h"

#include "PCLViewer.h"

class PCLViewer;
class PointCloud : public QMainWindow
{
	Q_OBJECT

public:
	PointCloud(QWidget *parent = Q_NULLPTR);

private:
	Ui::PointCloudClass ui;

	pcl::visualization::PCLVisualizer::Ptr _viewer;
	pcl::PointCloud<pcl::PointXYZRGBA>::Ptr _cloud;
	PCLViewer _pclViewer;
private slots:
	void OnReadFile();								//	��pcd�ļ��ۺ���
	void PCL();							//	pclֱ����ʾ
	void SaveAsPlY();					//		����Ϊply��ʽ

};
