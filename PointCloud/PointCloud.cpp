#include "stdafx.h"

#include "FileDeal.h"

#include "PointCloud.h"

PointCloud::PointCloud(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);



	this->setWindowTitle("PCL viewer");

	// Setup the cloud pointer
	_cloud.reset(new pcl::PointCloud<pcl::PointXYZRGBA>);
	
	// Set up the QVTK window
	_viewer.reset(new pcl::visualization::PCLVisualizer("viewer", false));
	ui.qvtkWidget->SetRenderWindow(_viewer->getRenderWindow());
	_viewer->setupInteractor(ui.qvtkWidget->GetInteractor(), ui.qvtkWidget->GetRenderWindow());
	ui.qvtkWidget->update();
	this->setWindowTitle("PCL viewer");
}


void PointCloud::OnReadFile()
{
	QString fileName = QFileDialog::getOpenFileName(this, tr("Open"), ".");
	if (fileName.isEmpty())
		return;
	fs::path fileFullName = fileName.toStdWString().c_str();
	if(fileFullName.extension() == ".pcd")
		_pclViewer.ReadPcdFile(_cloud, fileFullName);
	else if (fileFullName.extension() == ".asc")
	{
		_pclViewer.AscToPcd(fileFullName);
		_pclViewer.ReadPcdFile(_cloud, fileFullName.filename().string() + ".pcd");
	}

	_viewer->addPointCloud(_cloud, "cloud");
	_viewer->resetCamera();
	ui.qvtkWidget->update();
}

void PointCloud::PCL()
{
	FileDeal f;
	f.fileChange();
	f.test();
}

void PointCloud::SaveAsPlY()
{

}