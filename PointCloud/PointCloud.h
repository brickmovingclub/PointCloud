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
	pcl::PointCloud<pcl::PointNormal>::Ptr getPointNormal();

private:
	Ui::PointCloudClass ui;

	pcl::visualization::PCLVisualizer::Ptr _viewer;
	pcl::PointCloud<pcl::PointXYZ>::Ptr _cloud;
	PCLViewer _pclViewer;
private slots:

	//	
	void OnReadFile();										//	��pcd�ļ��ۺ���
	void  open_pcd_file();
	void	OnClear();
	void SaveAsPlY();										//		����Ϊply��ʽ
	void OnActionSearchKNear();

	/**********KT��ʵ�ֵ������ݵĿռ仮��*************/
		
	void	greedyTriangulation_reconstruct();				//	̰�������㷨
	void	poisson_reconstruct();							//	�����ؽ�

	/*********�˲���ʵ�ֵ������ݵĿռ仮��*********/
	void SearchKNear(float x, float y, float z, int &k);									//	�����˲�����ָ������Ұ뾶Ϊk�����
	void ShowLeafNode();										//	��ʾ�ռ仮�ֵ�Ҷ�ӽڵ�

};
