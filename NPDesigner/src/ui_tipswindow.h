/********************************************************************************
** Form generated from reading UI file 'tipswindow.ui'
**
** Created by: Qt User Interface Compiler version 5.15.3
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_TIPSWINDOW_H
#define UI_TIPSWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QLabel>
#include <QtWidgets/QTextEdit>

QT_BEGIN_NAMESPACE

class Ui_tipsWindow
{
public:
    QDialogButtonBox *buttonBox;
    QLabel *label;
    QTextEdit *textEdit;

    void setupUi(QDialog *tipsWindow)
    {
        if (tipsWindow->objectName().isEmpty())
            tipsWindow->setObjectName(QString::fromUtf8("tipsWindow"));
        tipsWindow->resize(467, 640);
        buttonBox = new QDialogButtonBox(tipsWindow);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setGeometry(QRect(100, 590, 341, 32));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);
        label = new QLabel(tipsWindow);
        label->setObjectName(QString::fromUtf8("label"));
        label->setGeometry(QRect(20, 30, 421, 51));
        textEdit = new QTextEdit(tipsWindow);
        textEdit->setObjectName(QString::fromUtf8("textEdit"));
        textEdit->setGeometry(QRect(30, 100, 411, 441));

        retranslateUi(tipsWindow);
        QObject::connect(buttonBox, SIGNAL(accepted()), tipsWindow, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), tipsWindow, SLOT(reject()));

        QMetaObject::connectSlotsByName(tipsWindow);
    } // setupUi

    void retranslateUi(QDialog *tipsWindow)
    {
        tipsWindow->setWindowTitle(QCoreApplication::translate("tipsWindow", "Dialog", nullptr));
        label->setText(QCoreApplication::translate("tipsWindow", "NPDesigner - quick tips", nullptr));
        textEdit->setHtml(QCoreApplication::translate("tipsWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'Ubuntu'; font-size:11pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Set the UA directory first to make sure all the generated pathnames are set to make UA look in the right place.</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Simple NPs: Add a bead type first, then add the bead at 0,0,0 using the NP bead button.</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; "
                        "margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Adding shells: This interface takes care of adding the positive bead for the outer layer and the negative bead for the inner layer.</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Adding brushes: This needs a pre-existing bead type.</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Saving: If the number of orientations is set to 1,"
                        " only a single file is created. Else N random orientations are saved with numbers 1...N , plus the original with number 0.</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Changes made directly to the tables aren't confirmed until you press the &quot;Update from tables&quot; button. </p></body></html>", nullptr));
    } // retranslateUi

};

namespace Ui {
    class tipsWindow: public Ui_tipsWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_TIPSWINDOW_H
