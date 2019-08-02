/*editbin /SUBSYSTEM:CONSOLE $(OUTDIR)\$(ProjectName).exe*/

#include "vtkAutoInit.h" 
VTK_MODULE_INIT(vtkRenderingOpenGL2); // VTK was built with vtkRenderingOpenGL2
VTK_MODULE_INIT(vtkInteractionStyle);

#include "stdafx.h"
#include "FileDeal.h"

#include "PointCloud.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	PointCloud w;
	w.show();
	return a.exec();
}
