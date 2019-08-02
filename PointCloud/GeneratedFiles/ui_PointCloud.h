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
#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_PointCloudClass
{
public:
    QMenuBar *menuBar;
    QToolBar *mainToolBar;
    QWidget *centralWidget;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *PointCloudClass)
    {
        if (PointCloudClass->objectName().isEmpty())
            PointCloudClass->setObjectName(QString::fromUtf8("PointCloudClass"));
        PointCloudClass->resize(600, 400);
        menuBar = new QMenuBar(PointCloudClass);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        PointCloudClass->setMenuBar(menuBar);
        mainToolBar = new QToolBar(PointCloudClass);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        PointCloudClass->addToolBar(mainToolBar);
        centralWidget = new QWidget(PointCloudClass);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        PointCloudClass->setCentralWidget(centralWidget);
        statusBar = new QStatusBar(PointCloudClass);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        PointCloudClass->setStatusBar(statusBar);

        retranslateUi(PointCloudClass);

        QMetaObject::connectSlotsByName(PointCloudClass);
    } // setupUi

    void retranslateUi(QMainWindow *PointCloudClass)
    {
        PointCloudClass->setWindowTitle(QApplication::translate("PointCloudClass", "PointCloud", nullptr));
    } // retranslateUi

};

namespace Ui {
    class PointCloudClass: public Ui_PointCloudClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_POINTCLOUD_H
