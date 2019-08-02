#include "stdafx.h"

#include "FileDeal.h"

#include "PointCloud.h"

PointCloud::PointCloud(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	FileDeal f;
	f.fileChange();
	f.test();
}
