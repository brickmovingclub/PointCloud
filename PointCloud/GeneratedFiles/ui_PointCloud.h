/********************************************************************************
** Form generated from reading UI file 'PointCloud.ui'
**
** Created by: Qt User Interface Compiler version 5.12.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_POINTCLOUD_H
#define UI_POINTCLOUD_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QWidget>
#include "QVTKWidget.h"

QT_BEGIN_NAMESPACE

class Ui_PointCloudClass
{
public:
    QAction *actionImport_file;
    QAction *actionPclShow;
    QWidget *centralWidget;
    QVTKWidget *qvtkWidget;
    QMenuBar *menuBar;
    QMenu *menuFile;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *PointCloudClass)
    {
        if (PointCloudClass->objectName().isEmpty())
            PointCloudClass->setObjectName(QString::fromUtf8("PointCloudClass"));
        PointCloudClass->resize(600, 400);
        actionImport_file = new QAction(PointCloudClass);
        actionImport_file->setObjectName(QString::fromUtf8("actionImport_file"));
        actionPclShow = new QAction(PointCloudClass);
        actionPclShow->setObjectName(QString::fromUtf8("actionPclShow"));
        centralWidget = new QWidget(PointCloudClass);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        qvtkWidget = new QVTKWidget(centralWidget);
        qvtkWidget->setObjectName(QString::fromUtf8("qvtkWidget"));
        qvtkWidget->setGeometry(QRect(-10, 9, 591, 361));
        PointCloudClass->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(PointCloudClass);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 600, 23));
        menuFile = new QMenu(menuBar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        PointCloudClass->setMenuBar(menuBar);
        mainToolBar = new QToolBar(PointCloudClass);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        PointCloudClass->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(PointCloudClass);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        PointCloudClass->setStatusBar(statusBar);

        menuBar->addAction(menuFile->menuAction());
        menuFile->addAction(actionImport_file);
        menuFile->addAction(actionPclShow);

        retranslateUi(PointCloudClass);
        QObject::connect(actionImport_file, SIGNAL(triggered(bool)), PointCloudClass, SLOT(OnReadFile()));
        QObject::connect(actionPclShow, SIGNAL(triggered(bool)), PointCloudClass, SLOT(PCL()));

        QMetaObject::connectSlotsByName(PointCloudClass);
    } // setupUi

    void retranslateUi(QMainWindow *PointCloudClass)
    {
        PointCloudClass->setWindowTitle(QApplication::translate("PointCloudClass", "PointCloud", nullptr));
        actionImport_file->setText(QApplication::translate("PointCloudClass", "Import file", nullptr));
        actionPclShow->setText(QApplication::translate("PointCloudClass", "PclShow", nullptr));
        menuFile->setTitle(QApplication::translate("PointCloudClass", "File", nullptr));
    } // retranslateUi

};

namespace Ui {
    class PointCloudClass: public Ui_PointCloudClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_POINTCLOUD_H
