/********************************************************************************
** Form generated from reading UI file 'Widget_SearchKNear.ui'
**
** Created by: Qt User Interface Compiler version 5.12.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_WIDGET_SEARCHKNEAR_H
#define UI_WIDGET_SEARCHKNEAR_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_Widget_SearchKNear
{
public:
    QWidget *layoutWidget;
    QHBoxLayout *horizontalLayout;
    QLabel *label_x;
    QLineEdit *lineEdit_X;
    QLabel *label_Y;
    QLineEdit *lineEdit_Y;
    QLabel *label_Z;
    QLineEdit *lineEdit_Z;
    QWidget *widget;
    QHBoxLayout *horizontalLayout_2;
    QPushButton *pushButton_voxelSearch;
    QPushButton *pushButton_kNearNodes;
    QPushButton *pushButton_NearRadiusSearch;
    QWidget *widget1;
    QHBoxLayout *horizontalLayout_3;
    QLabel *label_K;
    QLineEdit *lineEdit;

    void setupUi(QWidget *Widget_SearchKNear)
    {
        if (Widget_SearchKNear->objectName().isEmpty())
            Widget_SearchKNear->setObjectName(QString::fromUtf8("Widget_SearchKNear"));
        Widget_SearchKNear->resize(250, 143);
        layoutWidget = new QWidget(Widget_SearchKNear);
        layoutWidget->setObjectName(QString::fromUtf8("layoutWidget"));
        layoutWidget->setGeometry(QRect(20, 10, 185, 22));
        horizontalLayout = new QHBoxLayout(layoutWidget);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        horizontalLayout->setContentsMargins(0, 0, 0, 0);
        label_x = new QLabel(layoutWidget);
        label_x->setObjectName(QString::fromUtf8("label_x"));

        horizontalLayout->addWidget(label_x);

        lineEdit_X = new QLineEdit(layoutWidget);
        lineEdit_X->setObjectName(QString::fromUtf8("lineEdit_X"));

        horizontalLayout->addWidget(lineEdit_X);

        label_Y = new QLabel(layoutWidget);
        label_Y->setObjectName(QString::fromUtf8("label_Y"));

        horizontalLayout->addWidget(label_Y);

        lineEdit_Y = new QLineEdit(layoutWidget);
        lineEdit_Y->setObjectName(QString::fromUtf8("lineEdit_Y"));

        horizontalLayout->addWidget(lineEdit_Y);

        label_Z = new QLabel(layoutWidget);
        label_Z->setObjectName(QString::fromUtf8("label_Z"));

        horizontalLayout->addWidget(label_Z);

        lineEdit_Z = new QLineEdit(layoutWidget);
        lineEdit_Z->setObjectName(QString::fromUtf8("lineEdit_Z"));

        horizontalLayout->addWidget(lineEdit_Z);

        widget = new QWidget(Widget_SearchKNear);
        widget->setObjectName(QString::fromUtf8("widget"));
        widget->setGeometry(QRect(0, 100, 250, 25));
        horizontalLayout_2 = new QHBoxLayout(widget);
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        horizontalLayout_2->setContentsMargins(0, 0, 0, 0);
        pushButton_voxelSearch = new QPushButton(widget);
        pushButton_voxelSearch->setObjectName(QString::fromUtf8("pushButton_voxelSearch"));

        horizontalLayout_2->addWidget(pushButton_voxelSearch);

        pushButton_kNearNodes = new QPushButton(widget);
        pushButton_kNearNodes->setObjectName(QString::fromUtf8("pushButton_kNearNodes"));

        horizontalLayout_2->addWidget(pushButton_kNearNodes);

        pushButton_NearRadiusSearch = new QPushButton(widget);
        pushButton_NearRadiusSearch->setObjectName(QString::fromUtf8("pushButton_NearRadiusSearch"));

        horizontalLayout_2->addWidget(pushButton_NearRadiusSearch);

        widget1 = new QWidget(Widget_SearchKNear);
        widget1->setObjectName(QString::fromUtf8("widget1"));
        widget1->setGeometry(QRect(21, 61, 195, 22));
        horizontalLayout_3 = new QHBoxLayout(widget1);
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        horizontalLayout_3->setContentsMargins(0, 0, 0, 0);
        label_K = new QLabel(widget1);
        label_K->setObjectName(QString::fromUtf8("label_K"));

        horizontalLayout_3->addWidget(label_K);

        lineEdit = new QLineEdit(widget1);
        lineEdit->setObjectName(QString::fromUtf8("lineEdit"));

        horizontalLayout_3->addWidget(lineEdit);


        retranslateUi(Widget_SearchKNear);
        QObject::connect(pushButton_voxelSearch, SIGNAL(clicked(bool)), Widget_SearchKNear, SLOT(OnPushButton_voxelSearch()));
        QObject::connect(pushButton_kNearNodes, SIGNAL(clicked(bool)), Widget_SearchKNear, SLOT(OnPushButton_kNearNodes()));
        QObject::connect(pushButton_NearRadiusSearch, SIGNAL(clicked(bool)), Widget_SearchKNear, SLOT(OnPushButton_NearRadiusSearch()));

        QMetaObject::connectSlotsByName(Widget_SearchKNear);
    } // setupUi

    void retranslateUi(QWidget *Widget_SearchKNear)
    {
        Widget_SearchKNear->setWindowTitle(QApplication::translate("Widget_SearchKNear", "Form", nullptr));
        label_x->setText(QApplication::translate("Widget_SearchKNear", "X:", nullptr));
        label_Y->setText(QApplication::translate("Widget_SearchKNear", "Y:", nullptr));
        label_Z->setText(QApplication::translate("Widget_SearchKNear", "Z:", nullptr));
        pushButton_voxelSearch->setText(QApplication::translate("Widget_SearchKNear", "\344\275\223\347\264\240\346\220\234\347\264\242", nullptr));
        pushButton_kNearNodes->setText(QApplication::translate("Widget_SearchKNear", "K\351\242\206\345\237\237\350\212\202\347\202\271\346\220\234\347\264\242", nullptr));
        pushButton_NearRadiusSearch->setText(QApplication::translate("Widget_SearchKNear", "\345\215\212\345\276\204\346\220\234\347\264\242", nullptr));
        label_K->setText(QApplication::translate("Widget_SearchKNear", "K\351\230\266\351\242\206\345\237\237\357\274\232", nullptr));
    } // retranslateUi

};

namespace Ui {
    class Widget_SearchKNear: public Ui_Widget_SearchKNear {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_WIDGET_SEARCHKNEAR_H
