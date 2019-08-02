#include "stdafx.h"
#include "PointCloud.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	PointCloud w;
	w.show();
	return a.exec();
}
