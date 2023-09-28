/********************************************************************************
** Form generated from reading UI file 'addbrush.ui'
**
** Created by: Qt User Interface Compiler version 5.15.3
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_ADDBRUSH_H
#define UI_ADDBRUSH_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>

QT_BEGIN_NAMESPACE

class Ui_AddBrush
{
public:
    QDialogButtonBox *buttonBox;
    QLabel *label_4;
    QDoubleSpinBox *brushOccupancy;
    QLabel *label_5;
    QLabel *label_6;
    QLineEdit *brushRadialDist;
    QLabel *label;
    QComboBox *brushBeadID;
    QLabel *label_7;
    QLineEdit *densityOutBox;
    QLineEdit *numbeadsOutBox;
    QLabel *label_8;
    QLineEdit *outerLayerTarget;
    QCheckBox *forceAttachBox;

    void setupUi(QDialog *AddBrush)
    {
        if (AddBrush->objectName().isEmpty())
            AddBrush->setObjectName(QString::fromUtf8("AddBrush"));
        AddBrush->resize(400, 369);
        buttonBox = new QDialogButtonBox(AddBrush);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setGeometry(QRect(40, 320, 341, 32));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);
        label_4 = new QLabel(AddBrush);
        label_4->setObjectName(QString::fromUtf8("label_4"));
        label_4->setGeometry(QRect(50, 30, 101, 17));
        brushOccupancy = new QDoubleSpinBox(AddBrush);
        brushOccupancy->setObjectName(QString::fromUtf8("brushOccupancy"));
        brushOccupancy->setGeometry(QRect(280, 150, 65, 26));
        brushOccupancy->setMaximum(1.000000000000000);
        brushOccupancy->setSingleStep(0.100000000000000);
        brushOccupancy->setValue(1.000000000000000);
        label_5 = new QLabel(AddBrush);
        label_5->setObjectName(QString::fromUtf8("label_5"));
        label_5->setGeometry(QRect(20, 60, 241, 41));
        label_6 = new QLabel(AddBrush);
        label_6->setObjectName(QString::fromUtf8("label_6"));
        label_6->setGeometry(QRect(60, 150, 151, 17));
        brushRadialDist = new QLineEdit(AddBrush);
        brushRadialDist->setObjectName(QString::fromUtf8("brushRadialDist"));
        brushRadialDist->setGeometry(QRect(300, 70, 71, 25));
        label = new QLabel(AddBrush);
        label->setObjectName(QString::fromUtf8("label"));
        label->setGeometry(QRect(20, 110, 271, 17));
        brushBeadID = new QComboBox(AddBrush);
        brushBeadID->setObjectName(QString::fromUtf8("brushBeadID"));
        brushBeadID->setGeometry(QRect(300, 30, 86, 25));
        label_7 = new QLabel(AddBrush);
        label_7->setObjectName(QString::fromUtf8("label_7"));
        label_7->setGeometry(QRect(20, 230, 221, 17));
        densityOutBox = new QLineEdit(AddBrush);
        densityOutBox->setObjectName(QString::fromUtf8("densityOutBox"));
        densityOutBox->setEnabled(true);
        densityOutBox->setGeometry(QRect(260, 220, 113, 25));
        densityOutBox->setAcceptDrops(true);
        densityOutBox->setReadOnly(true);
        numbeadsOutBox = new QLineEdit(AddBrush);
        numbeadsOutBox->setObjectName(QString::fromUtf8("numbeadsOutBox"));
        numbeadsOutBox->setEnabled(true);
        numbeadsOutBox->setGeometry(QRect(260, 270, 113, 25));
        numbeadsOutBox->setAcceptDrops(true);
        numbeadsOutBox->setReadOnly(true);
        label_8 = new QLabel(AddBrush);
        label_8->setObjectName(QString::fromUtf8("label_8"));
        label_8->setGeometry(QRect(20, 270, 221, 17));
        outerLayerTarget = new QLineEdit(AddBrush);
        outerLayerTarget->setObjectName(QString::fromUtf8("outerLayerTarget"));
        outerLayerTarget->setEnabled(true);
        outerLayerTarget->setGeometry(QRect(260, 110, 113, 25));
        outerLayerTarget->setAcceptDrops(true);
        outerLayerTarget->setReadOnly(true);
        forceAttachBox = new QCheckBox(AddBrush);
        forceAttachBox->setObjectName(QString::fromUtf8("forceAttachBox"));
        forceAttachBox->setGeometry(QRect(30, 190, 171, 23));

        retranslateUi(AddBrush);
        QObject::connect(buttonBox, SIGNAL(accepted()), AddBrush, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), AddBrush, SLOT(reject()));

        QMetaObject::connectSlotsByName(AddBrush);
    } // setupUi

    void retranslateUi(QDialog *AddBrush)
    {
        AddBrush->setWindowTitle(QCoreApplication::translate("AddBrush", "Dialog", nullptr));
        label_4->setText(QCoreApplication::translate("AddBrush", "BeadTypeID", nullptr));
        label_5->setText(QCoreApplication::translate("AddBrush", "Layer center radial distance from 0", nullptr));
        label_6->setText(QCoreApplication::translate("AddBrush", "Brush occupancy", nullptr));
        brushRadialDist->setText(QCoreApplication::translate("AddBrush", "5", nullptr));
        label->setText(QCoreApplication::translate("AddBrush", "Suggested for outer layer:", nullptr));
        label_7->setText(QCoreApplication::translate("AddBrush", "Expected density (beads/nm^2)", nullptr));
        densityOutBox->setText(QCoreApplication::translate("AddBrush", "1", nullptr));
        numbeadsOutBox->setText(QCoreApplication::translate("AddBrush", "1", nullptr));
        label_8->setText(QCoreApplication::translate("AddBrush", "Expected total beads", nullptr));
        outerLayerTarget->setText(QCoreApplication::translate("AddBrush", "1", nullptr));
        forceAttachBox->setText(QCoreApplication::translate("AddBrush", "Allow floating beads", nullptr));
    } // retranslateUi

};

namespace Ui {
    class AddBrush: public Ui_AddBrush {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_ADDBRUSH_H
