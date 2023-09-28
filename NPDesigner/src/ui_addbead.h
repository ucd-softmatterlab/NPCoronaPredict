/********************************************************************************
** Form generated from reading UI file 'addbead.ui'
**
** Created by: Qt User Interface Compiler version 5.15.3
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_ADDBEAD_H
#define UI_ADDBEAD_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>

QT_BEGIN_NAMESPACE

class Ui_AddBead
{
public:
    QDialogButtonBox *buttonBox;
    QLineEdit *newBeadX;
    QLabel *label;
    QLabel *label_2;
    QLineEdit *newBeadY;
    QLabel *label_3;
    QLabel *label_4;
    QLineEdit *newBeadZ;
    QComboBox *beadTypeID;

    void setupUi(QDialog *AddBead)
    {
        if (AddBead->objectName().isEmpty())
            AddBead->setObjectName(QString::fromUtf8("AddBead"));
        AddBead->resize(400, 300);
        buttonBox = new QDialogButtonBox(AddBead);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setGeometry(QRect(30, 240, 341, 32));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);
        newBeadX = new QLineEdit(AddBead);
        newBeadX->setObjectName(QString::fromUtf8("newBeadX"));
        newBeadX->setGeometry(QRect(20, 120, 113, 25));
        label = new QLabel(AddBead);
        label->setObjectName(QString::fromUtf8("label"));
        label->setGeometry(QRect(50, 90, 21, 17));
        label_2 = new QLabel(AddBead);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setGeometry(QRect(180, 90, 21, 17));
        newBeadY = new QLineEdit(AddBead);
        newBeadY->setObjectName(QString::fromUtf8("newBeadY"));
        newBeadY->setGeometry(QRect(150, 120, 113, 25));
        label_3 = new QLabel(AddBead);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        label_3->setGeometry(QRect(310, 90, 21, 17));
        label_4 = new QLabel(AddBead);
        label_4->setObjectName(QString::fromUtf8("label_4"));
        label_4->setGeometry(QRect(80, 60, 101, 17));
        newBeadZ = new QLineEdit(AddBead);
        newBeadZ->setObjectName(QString::fromUtf8("newBeadZ"));
        newBeadZ->setGeometry(QRect(270, 120, 113, 25));
        beadTypeID = new QComboBox(AddBead);
        beadTypeID->setObjectName(QString::fromUtf8("beadTypeID"));
        beadTypeID->setGeometry(QRect(210, 50, 86, 25));

        retranslateUi(AddBead);
        QObject::connect(buttonBox, SIGNAL(accepted()), AddBead, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), AddBead, SLOT(reject()));

        QMetaObject::connectSlotsByName(AddBead);
    } // setupUi

    void retranslateUi(QDialog *AddBead)
    {
        AddBead->setWindowTitle(QCoreApplication::translate("AddBead", "Dialog", nullptr));
        newBeadX->setText(QCoreApplication::translate("AddBead", "0", nullptr));
        label->setText(QCoreApplication::translate("AddBead", "x", nullptr));
        label_2->setText(QCoreApplication::translate("AddBead", "y", nullptr));
        newBeadY->setText(QCoreApplication::translate("AddBead", "0", nullptr));
        label_3->setText(QCoreApplication::translate("AddBead", "z", nullptr));
        label_4->setText(QCoreApplication::translate("AddBead", "BeadTypeID", nullptr));
        newBeadZ->setText(QCoreApplication::translate("AddBead", "0", nullptr));
    } // retranslateUi

};

namespace Ui {
    class AddBead: public Ui_AddBead {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_ADDBEAD_H
