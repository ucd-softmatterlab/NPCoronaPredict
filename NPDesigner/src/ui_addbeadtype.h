/********************************************************************************
** Form generated from reading UI file 'addbeadtype.ui'
**
** Created by: Qt User Interface Compiler version 5.15.3
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_ADDBEADTYPE_H
#define UI_ADDBEADTYPE_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPlainTextEdit>
#include <QtWidgets/QPushButton>

QT_BEGIN_NAMESPACE

class Ui_AddBeadType
{
public:
    QDialogButtonBox *buttonBox;
    QPushButton *surfFindDirButton;
    QLabel *label_2;
    QLabel *label_5;
    QLabel *label_6;
    QLabel *label_3;
    QPlainTextEdit *hamakerFile;
    QPushButton *hamFindButton;
    QPlainTextEdit *surfDir;
    QLabel *label_8;
    QLabel *label_4;
    QLabel *label_7;
    QLineEdit *ljCutoff;
    QLineEdit *surfFactor;
    QLineEdit *coreFactor;
    QLineEdit *zeta;
    QLineEdit *radius;

    void setupUi(QDialog *AddBeadType)
    {
        if (AddBeadType->objectName().isEmpty())
            AddBeadType->setObjectName(QString::fromUtf8("AddBeadType"));
        AddBeadType->resize(400, 507);
        buttonBox = new QDialogButtonBox(AddBeadType);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setGeometry(QRect(30, 450, 341, 32));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);
        surfFindDirButton = new QPushButton(AddBeadType);
        surfFindDirButton->setObjectName(QString::fromUtf8("surfFindDirButton"));
        surfFindDirButton->setGeometry(QRect(290, 50, 89, 25));
        label_2 = new QLabel(AddBeadType);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setGeometry(QRect(20, 180, 131, 17));
        label_5 = new QLabel(AddBeadType);
        label_5->setObjectName(QString::fromUtf8("label_5"));
        label_5->setGeometry(QRect(20, 90, 171, 17));
        label_6 = new QLabel(AddBeadType);
        label_6->setObjectName(QString::fromUtf8("label_6"));
        label_6->setGeometry(QRect(30, 260, 131, 17));
        label_3 = new QLabel(AddBeadType);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        label_3->setGeometry(QRect(20, 220, 131, 17));
        hamakerFile = new QPlainTextEdit(AddBeadType);
        hamakerFile->setObjectName(QString::fromUtf8("hamakerFile"));
        hamakerFile->setGeometry(QRect(20, 120, 251, 31));
        hamFindButton = new QPushButton(AddBeadType);
        hamFindButton->setObjectName(QString::fromUtf8("hamFindButton"));
        hamFindButton->setGeometry(QRect(290, 120, 89, 25));
        surfDir = new QPlainTextEdit(AddBeadType);
        surfDir->setObjectName(QString::fromUtf8("surfDir"));
        surfDir->setGeometry(QRect(20, 50, 251, 31));
        label_8 = new QLabel(AddBeadType);
        label_8->setObjectName(QString::fromUtf8("label_8"));
        label_8->setGeometry(QRect(20, 370, 131, 17));
        label_4 = new QLabel(AddBeadType);
        label_4->setObjectName(QString::fromUtf8("label_4"));
        label_4->setGeometry(QRect(20, 20, 171, 17));
        label_7 = new QLabel(AddBeadType);
        label_7->setObjectName(QString::fromUtf8("label_7"));
        label_7->setGeometry(QRect(20, 320, 131, 17));
        ljCutoff = new QLineEdit(AddBeadType);
        ljCutoff->setObjectName(QString::fromUtf8("ljCutoff"));
        ljCutoff->setGeometry(QRect(240, 260, 113, 25));
        surfFactor = new QLineEdit(AddBeadType);
        surfFactor->setObjectName(QString::fromUtf8("surfFactor"));
        surfFactor->setGeometry(QRect(240, 320, 113, 25));
        coreFactor = new QLineEdit(AddBeadType);
        coreFactor->setObjectName(QString::fromUtf8("coreFactor"));
        coreFactor->setGeometry(QRect(250, 370, 113, 25));
        zeta = new QLineEdit(AddBeadType);
        zeta->setObjectName(QString::fromUtf8("zeta"));
        zeta->setGeometry(QRect(240, 220, 113, 25));
        radius = new QLineEdit(AddBeadType);
        radius->setObjectName(QString::fromUtf8("radius"));
        radius->setGeometry(QRect(240, 180, 113, 25));

        retranslateUi(AddBeadType);
        QObject::connect(buttonBox, SIGNAL(accepted()), AddBeadType, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), AddBeadType, SLOT(reject()));

        QMetaObject::connectSlotsByName(AddBeadType);
    } // setupUi

    void retranslateUi(QDialog *AddBeadType)
    {
        AddBeadType->setWindowTitle(QCoreApplication::translate("AddBeadType", "Dialog", nullptr));
        surfFindDirButton->setText(QCoreApplication::translate("AddBeadType", "Find", nullptr));
        label_2->setText(QCoreApplication::translate("AddBeadType", "Radius [nm]", nullptr));
        label_5->setText(QCoreApplication::translate("AddBeadType", "Hamaker File", nullptr));
        label_6->setText(QCoreApplication::translate("AddBeadType", "LJ cutoff [nm]", nullptr));
        label_3->setText(QCoreApplication::translate("AddBeadType", "Zeta potential [mV]", nullptr));
        hamFindButton->setText(QCoreApplication::translate("AddBeadType", "Find", nullptr));
        label_8->setText(QCoreApplication::translate("AddBeadType", "Hamaker factor", nullptr));
        label_4->setText(QCoreApplication::translate("AddBeadType", "Surface PMF Directory", nullptr));
        label_7->setText(QCoreApplication::translate("AddBeadType", "Surface factor", nullptr));
        ljCutoff->setText(QCoreApplication::translate("AddBeadType", "1", nullptr));
        surfFactor->setText(QCoreApplication::translate("AddBeadType", "1", nullptr));
        coreFactor->setText(QCoreApplication::translate("AddBeadType", "1", nullptr));
        zeta->setText(QCoreApplication::translate("AddBeadType", "0", nullptr));
        radius->setText(QCoreApplication::translate("AddBeadType", "5", nullptr));
    } // retranslateUi

};

namespace Ui {
    class AddBeadType: public Ui_AddBeadType {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_ADDBEADTYPE_H
