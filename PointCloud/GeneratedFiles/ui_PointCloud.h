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
    QAction *actionSave_as;
    QAction *action_drawBox;
    QAction *actiontanlan;
    QAction *actionbosong;
    QAction *actionyidong;
    QAction *actionoepnPcdFile;
    QAction *actionqing_kong;
    QAction *actionSearchKNear;
    QAction *actionSpaceDiv;
    QAction *actionShowLeafNode;
    QAction *actionOctreeSpaceDev;
    QWidget *centralWidget;
    QVTKWidget *qvtkWidget;
    QMenuBar *menuBar;
    QMenu *menuFile;
    QMenu *menuEdit;
    QMenu *menu;
    QMenu *menuOcTree;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *PointCloudClass)
    {
        if (PointCloudClass->objectName().isEmpty())
            PointCloudClass->setObjectName(QString::fromUtf8("PointCloudClass"));
        PointCloudClass->resize(981, 799);
        actionImport_file = new QAction(PointCloudClass);
        actionImport_file->setObjectName(QString::fromUtf8("actionImport_file"));
        actionPclShow = new QAction(PointCloudClass);
        actionPclShow->setObjectName(QString::fromUtf8("actionPclShow"));
        actionSave_as = new QAction(PointCloudClass);
        actionSave_as->setObjectName(QString::fromUtf8("actionSave_as"));
        action_drawBox = new QAction(PointCloudClass);
        action_drawBox->setObjectName(QString::fromUtf8("action_drawBox"));
        actiontanlan = new QAction(PointCloudClass);
        actiontanlan->setObjectName(QString::fromUtf8("actiontanlan"));
        actionbosong = new QAction(PointCloudClass);
        actionbosong->setObjectName(QString::fromUtf8("actionbosong"));
        actionyidong = new QAction(PointCloudClass);
        actionyidong->setObjectName(QString::fromUtf8("actionyidong"));
        actionoepnPcdFile = new QAction(PointCloudClass);
        actionoepnPcdFile->setObjectName(QString::fromUtf8("actionoepnPcdFile"));
        actionqing_kong = new QAction(PointCloudClass);
        actionqing_kong->setObjectName(QString::fromUtf8("actionqing_kong"));
        actionSearchKNear = new QAction(PointCloudClass);
        actionSearchKNear->setObjectName(QString::fromUtf8("actionSearchKNear"));
        actionSpaceDiv = new QAction(PointCloudClass);
        actionSpaceDiv->setObjectName(QString::fromUtf8("actionSpaceDiv"));
        actionShowLeafNode = new QAction(PointCloudClass);
        actionShowLeafNode->setObjectName(QString::fromUtf8("actionShowLeafNode"));
        actionOctreeSpaceDev = new QAction(PointCloudClass);
        actionOctreeSpaceDev->setObjectName(QString::fromUtf8("actionOctreeSpaceDev"));
        centralWidget = new QWidget(PointCloudClass);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        qvtkWidget = new QVTKWidget(centralWidget);
        qvtkWidget->setObjectName(QString::fromUtf8("qvtkWidget"));
        qvtkWidget->setGeometry(QRect(10, 0, 911, 731));
        PointCloudClass->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(PointCloudClass);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 981, 23));
        menuFile = new QMenu(menuBar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        menuEdit = new QMenu(menuBar);
        menuEdit->setObjectName(QString::fromUtf8("menuEdit"));
        menu = new QMenu(menuBar);
        menu->setObjectName(QString::fromUtf8("menu"));
        menuOcTree = new QMenu(menuBar);
        menuOcTree->setObjectName(QString::fromUtf8("menuOcTree"));
        PointCloudClass->setMenuBar(menuBar);
        mainToolBar = new QToolBar(PointCloudClass);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        PointCloudClass->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(PointCloudClass);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        PointCloudClass->setStatusBar(statusBar);

        menuBar->addAction(menuFile->menuAction());
        menuBar->addAction(menuEdit->menuAction());
        menuBar->addAction(menu->menuAction());
        menuBar->addAction(menuOcTree->menuAction());
        menuFile->addAction(actionImport_file);
        menuFile->addAction(actionPclShow);
        menuFile->addAction(actionSave_as);
        menuFile->addAction(actionoepnPcdFile);
        menuFile->addAction(actionqing_kong);
        menuEdit->addAction(action_drawBox);
        menu->addAction(actiontanlan);
        menu->addAction(actionbosong);
        menuOcTree->addAction(actionSearchKNear);
        menuOcTree->addAction(actionShowLeafNode);
        menuOcTree->addAction(actionOctreeSpaceDev);

        retranslateUi(PointCloudClass);
        QObject::connect(actionImport_file, SIGNAL(triggered(bool)), PointCloudClass, SLOT(OnReadFile()));
        QObject::connect(actionSave_as, SIGNAL(triggered(bool)), PointCloudClass, SLOT(SaveAsPlY()));
        QObject::connect(actionbosong, SIGNAL(triggered(bool)), PointCloudClass, SLOT(poisson_reconstruct()));
        QObject::connect(actiontanlan, SIGNAL(triggered(bool)), PointCloudClass, SLOT(greedyTriangulation_reconstruct()));
        QObject::connect(actionoepnPcdFile, SIGNAL(triggered(bool)), PointCloudClass, SLOT(open_pcd_file()));
        QObject::connect(actionSearchKNear, SIGNAL(triggered(bool)), PointCloudClass, SLOT(OnSearchKNear()));
        QObject::connect(actionShowLeafNode, SIGNAL(triggered(bool)), PointCloudClass, SLOT(ShowLeafNode()));
        QObject::connect(action_drawBox, SIGNAL(triggered(bool)), PointCloudClass, SLOT(DrawBoundingBox()));
        QObject::connect(actionOctreeSpaceDev, SIGNAL(triggered(bool)), PointCloudClass, SLOT(Triangulation()));

        QMetaObject::connectSlotsByName(PointCloudClass);
    } // setupUi

    void retranslateUi(QMainWindow *PointCloudClass)
    {
        PointCloudClass->setWindowTitle(QApplication::translate("PointCloudClass", "PointCloud", nullptr));
        actionImport_file->setText(QApplication::translate("PointCloudClass", "Import file", nullptr));
        actionPclShow->setText(QApplication::translate("PointCloudClass", "PclShow", nullptr));
        actionSave_as->setText(QApplication::translate("PointCloudClass", "Save as", nullptr));
        action_drawBox->setText(QApplication::translate("PointCloudClass", "Draw Box", nullptr));
        actiontanlan->setText(QApplication::translate("PointCloudClass", "\350\264\252\345\251\252\344\270\211\350\247\222\347\256\227\346\263\225", nullptr));
        actionbosong->setText(QApplication::translate("PointCloudClass", "\346\263\212\346\235\276\347\256\227\346\263\225", nullptr));
        actionyidong->setText(QApplication::translate("PointCloudClass", "\347\247\273\345\212\250\347\253\213\344\275\223\347\256\227\346\263\225", nullptr));
        actionoepnPcdFile->setText(QApplication::translate("PointCloudClass", "\346\211\223\345\274\200pcd\346\226\207\344\273\266", nullptr));
        actionqing_kong->setText(QApplication::translate("PointCloudClass", "\346\270\205\347\251\272", nullptr));
        actionSearchKNear->setText(QApplication::translate("PointCloudClass", "\346\220\234\347\264\242k\351\230\266\351\242\206\345\237\237\347\202\271", nullptr));
        actionSpaceDiv->setText(QApplication::translate("PointCloudClass", "\347\251\272\351\227\264\345\210\222\345\210\206", nullptr));
        actionShowLeafNode->setText(QApplication::translate("PointCloudClass", "\346\230\276\347\244\272\345\217\266\345\255\220\350\212\202\347\202\271", nullptr));
        actionOctreeSpaceDev->setText(QApplication::translate("PointCloudClass", "\345\205\253\345\217\211\346\240\221\347\251\272\351\227\264\345\210\222\345\210\206", nullptr));
        menuFile->setTitle(QApplication::translate("PointCloudClass", "File", nullptr));
        menuEdit->setTitle(QApplication::translate("PointCloudClass", "Edit", nullptr));
        menu->setTitle(QApplication::translate("PointCloudClass", "KtTree", nullptr));
        menuOcTree->setTitle(QApplication::translate("PointCloudClass", "OcTree", nullptr));
    } // retranslateUi

};

namespace Ui {
    class PointCloudClass: public Ui_PointCloudClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_POINTCLOUD_H
