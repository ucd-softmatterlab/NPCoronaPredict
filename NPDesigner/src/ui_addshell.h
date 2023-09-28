/********************************************************************************
** Form generated from reading UI file 'addshell.ui'
**
** Created by: Qt User Interface Compiler version 5.15.3
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_ADDSHELL_H
#define UI_ADDSHELL_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPlainTextEdit>
#include <QtWidgets/QPushButton>

QT_BEGIN_NAMESPACE

class Ui_AddShell
{
public:
    QDialogButtonBox *buttonBox;
    QPushButton *findSurf;
    QPlainTextEdit *surfDirIn;
    QLabel *label_3;
    QLabel *label_4;
    QLabel *label_5;
    QLabel *label_2;
    QLabel *label_6;
    QPushButton *findHamaker;
    QPlainTextEdit *hamakerIn;
    QLabel *label_7;
    QLineEdit *innerRadius;
    QLineEdit *outerRadius;
    QLineEdit *zeta;
    QLineEdit *ljCutoff;

    void setupUi(QDialog *AddShell)
    {
        if (AddShell->objectName().isEmpty())
            AddShell->setObjectName(QString::fromUtf8("AddShell"));
        AddShell->resize(400, 416);
        buttonBox = new QDialogButtonBox(AddShell);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setGeometry(QRect(40, 340, 341, 32));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);
        findSurf = new QPushButton(AddShell);
        findSurf->setObjectName(QString::fromUtf8("findSurf"));
        findSurf->setGeometry(QRect(300, 40, 89, 25));
        surfDirIn = new QPlainTextEdit(AddShell);
        surfDirIn->setObjectName(QString::fromUtf8("surfDirIn"));
        surfDirIn->setGeometry(QRect(30, 40, 251, 31));
        label_3 = new QLabel(AddShell);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        label_3->setGeometry(QRect(30, 250, 131, 17));
        label_4 = new QLabel(AddShell);
        label_4->setObjectName(QString::fromUtf8("label_4"));
        label_4->setGeometry(QRect(30, 10, 171, 17));
        label_5 = new QLabel(AddShell);
        label_5->setObjectName(QString::fromUtf8("label_5"));
        label_5->setGeometry(QRect(30, 80, 171, 17));
        label_2 = new QLabel(AddShell);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setGeometry(QRect(30, 200, 131, 17));
        label_6 = new QLabel(AddShell);
        label_6->setObjectName(QString::fromUtf8("label_6"));
        label_6->setGeometry(QRect(40, 290, 131, 17));
        findHamaker = new QPushButton(AddShell);
        findHamaker->setObjectName(QString::fromUtf8("findHamaker"));
        findHamaker->setGeometry(QRect(300, 110, 89, 25));
        hamakerIn = new QPlainTextEdit(AddShell);
        hamakerIn->setObjectName(QString::fromUtf8("hamakerIn"));
        hamakerIn->setGeometry(QRect(30, 110, 251, 31));
        label_7 = new QLabel(AddShell);
        label_7->setObjectName(QString::fromUtf8("label_7"));
        label_7->setGeometry(QRect(30, 160, 131, 17));
        innerRadius = new QLineEdit(AddShell);
        innerRadius->setObjectName(QString::fromUtf8("innerRadius"));
        innerRadius->setGeometry(QRect(260, 150, 113, 25));
        outerRadius = new QLineEdit(AddShell);
        outerRadius->setObjectName(QString::fromUtf8("outerRadius"));
        outerRadius->setGeometry(QRect(260, 190, 113, 25));
        zeta = new QLineEdit(AddShell);
        zeta->setObjectName(QString::fromUtf8("zeta"));
        zeta->setGeometry(QRect(250, 250, 113, 25));
        ljCutoff = new QLineEdit(AddShell);
        ljCutoff->setObjectName(QString::fromUtf8("ljCutoff"));
        ljCutoff->setGeometry(QRect(240, 290, 113, 25));

        retranslateUi(AddShell);
        QObject::connect(buttonBox, SIGNAL(accepted()), AddShell, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), AddShell, SLOT(reject()));

        QMetaObject::connectSlotsByName(AddShell);
    } // setupUi

    void retranslateUi(QDialog *AddShell)
    {
        AddShell->setWindowTitle(QCoreApplication::translate("AddShell", "Dialog", nullptr));
        findSurf->setText(QCoreApplication::translate("AddShell", "Find", nullptr));
        label_3->setText(QCoreApplication::translate("AddShell", "Zeta potential [mV]", nullptr));
        label_4->setText(QCoreApplication::translate("AddShell", "Surface PMF Directory", nullptr));
        label_5->setText(QCoreApplication::translate("AddShell", "Hamaker File", nullptr));
        label_2->setText(QCoreApplication::translate("AddShell", "Outer Radius [nm]", nullptr));
        label_6->setText(QCoreApplication::translate("AddShell", "LJ cutoff [nm]", nullptr));
        findHamaker->setText(QCoreApplication::translate("AddShell", "Find", nullptr));
        label_7->setText(QCoreApplication::translate("AddShell", "Inner Radius [nm]", nullptr));
        innerRadius->setText(QCoreApplication::translate("AddShell", "2", nullptr));
        outerRadius->setText(QCoreApplication::translate("AddShell", "5", nullptr));
        zeta->setText(QCoreApplication::translate("AddShell", "0", nullptr));
        ljCutoff->setText(QCoreApplication::translate("AddShell", "1", nullptr));
    } // retranslateUi

};

namespace Ui {
    class AddShell: public Ui_AddShell {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_ADDSHELL_H
