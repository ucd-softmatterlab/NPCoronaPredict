/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.15.3
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QGraphicsView>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *actionNew;
    QAction *actionQuit;
    QAction *actionTips;
    QAction *actionLoad;
    QWidget *centralwidget;
    QHBoxLayout *horizontalLayout_2;
    QGridLayout *gridLayout;
    QPushButton *newBeadTypeButton;
    QGraphicsView *graphicsView;
    QTableWidget *beadTypeTable;
    QPushButton *recenterNP;
    QHBoxLayout *horizontalLayout;
    QLabel *label;
    QLineEdit *uaPath;
    QPushButton *findUADir;
    QPushButton *newBeadButton;
    QSpinBox *numNPOrients;
    QPushButton *updateTables;
    QPushButton *newBrushButton;
    QPushButton *newShellButton;
    QLabel *label_2;
    QTableWidget *beadTable;
    QVBoxLayout *verticalLayout;
    QPushButton *saveNPButton;
    QHBoxLayout *horizontalLayout_6;
    QCheckBox *autoBounds;
    QLabel *label_3;
    QLineEdit *lowerBoundLine;
    QLabel *label_4;
    QLineEdit *upperBoundLine;
    QSpacerItem *horizontalSpacer_2;
    QSpacerItem *horizontalSpacer;
    QMenuBar *menubar;
    QMenu *menuFile;
    QMenu *menuHelp;
    QStatusBar *statusbar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(1262, 842);
        actionNew = new QAction(MainWindow);
        actionNew->setObjectName(QString::fromUtf8("actionNew"));
        actionQuit = new QAction(MainWindow);
        actionQuit->setObjectName(QString::fromUtf8("actionQuit"));
        actionTips = new QAction(MainWindow);
        actionTips->setObjectName(QString::fromUtf8("actionTips"));
        actionLoad = new QAction(MainWindow);
        actionLoad->setObjectName(QString::fromUtf8("actionLoad"));
        centralwidget = new QWidget(MainWindow);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        horizontalLayout_2 = new QHBoxLayout(centralwidget);
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        newBeadTypeButton = new QPushButton(centralwidget);
        newBeadTypeButton->setObjectName(QString::fromUtf8("newBeadTypeButton"));

        gridLayout->addWidget(newBeadTypeButton, 2, 1, 1, 1);

        graphicsView = new QGraphicsView(centralwidget);
        graphicsView->setObjectName(QString::fromUtf8("graphicsView"));

        gridLayout->addWidget(graphicsView, 0, 0, 1, 4);

        beadTypeTable = new QTableWidget(centralwidget);
        if (beadTypeTable->columnCount() < 9)
            beadTypeTable->setColumnCount(9);
        QTableWidgetItem *__qtablewidgetitem = new QTableWidgetItem();
        beadTypeTable->setHorizontalHeaderItem(0, __qtablewidgetitem);
        QTableWidgetItem *__qtablewidgetitem1 = new QTableWidgetItem();
        beadTypeTable->setHorizontalHeaderItem(1, __qtablewidgetitem1);
        QTableWidgetItem *__qtablewidgetitem2 = new QTableWidgetItem();
        beadTypeTable->setHorizontalHeaderItem(2, __qtablewidgetitem2);
        QTableWidgetItem *__qtablewidgetitem3 = new QTableWidgetItem();
        beadTypeTable->setHorizontalHeaderItem(3, __qtablewidgetitem3);
        QTableWidgetItem *__qtablewidgetitem4 = new QTableWidgetItem();
        beadTypeTable->setHorizontalHeaderItem(4, __qtablewidgetitem4);
        QTableWidgetItem *__qtablewidgetitem5 = new QTableWidgetItem();
        beadTypeTable->setHorizontalHeaderItem(5, __qtablewidgetitem5);
        QTableWidgetItem *__qtablewidgetitem6 = new QTableWidgetItem();
        beadTypeTable->setHorizontalHeaderItem(6, __qtablewidgetitem6);
        QTableWidgetItem *__qtablewidgetitem7 = new QTableWidgetItem();
        beadTypeTable->setHorizontalHeaderItem(7, __qtablewidgetitem7);
        QTableWidgetItem *__qtablewidgetitem8 = new QTableWidgetItem();
        beadTypeTable->setHorizontalHeaderItem(8, __qtablewidgetitem8);
        beadTypeTable->setObjectName(QString::fromUtf8("beadTypeTable"));
        beadTypeTable->setAlternatingRowColors(true);

        gridLayout->addWidget(beadTypeTable, 4, 0, 1, 6);

        recenterNP = new QPushButton(centralwidget);
        recenterNP->setObjectName(QString::fromUtf8("recenterNP"));

        gridLayout->addWidget(recenterNP, 3, 0, 1, 1);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        label = new QLabel(centralwidget);
        label->setObjectName(QString::fromUtf8("label"));

        horizontalLayout->addWidget(label);

        uaPath = new QLineEdit(centralwidget);
        uaPath->setObjectName(QString::fromUtf8("uaPath"));

        horizontalLayout->addWidget(uaPath);

        findUADir = new QPushButton(centralwidget);
        findUADir->setObjectName(QString::fromUtf8("findUADir"));

        horizontalLayout->addWidget(findUADir);


        gridLayout->addLayout(horizontalLayout, 2, 5, 1, 1);

        newBeadButton = new QPushButton(centralwidget);
        newBeadButton->setObjectName(QString::fromUtf8("newBeadButton"));
        newBeadButton->setEnabled(false);

        gridLayout->addWidget(newBeadButton, 2, 0, 1, 1);

        numNPOrients = new QSpinBox(centralwidget);
        numNPOrients->setObjectName(QString::fromUtf8("numNPOrients"));
        numNPOrients->setMinimum(1);
        numNPOrients->setMaximum(999);

        gridLayout->addWidget(numNPOrients, 5, 2, 1, 1);

        updateTables = new QPushButton(centralwidget);
        updateTables->setObjectName(QString::fromUtf8("updateTables"));

        gridLayout->addWidget(updateTables, 3, 1, 1, 1);

        newBrushButton = new QPushButton(centralwidget);
        newBrushButton->setObjectName(QString::fromUtf8("newBrushButton"));
        newBrushButton->setEnabled(false);

        gridLayout->addWidget(newBrushButton, 2, 3, 1, 1);

        newShellButton = new QPushButton(centralwidget);
        newShellButton->setObjectName(QString::fromUtf8("newShellButton"));

        gridLayout->addWidget(newShellButton, 2, 2, 1, 1);

        label_2 = new QLabel(centralwidget);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        gridLayout->addWidget(label_2, 5, 1, 1, 1);

        beadTable = new QTableWidget(centralwidget);
        if (beadTable->columnCount() < 4)
            beadTable->setColumnCount(4);
        QTableWidgetItem *__qtablewidgetitem9 = new QTableWidgetItem();
        beadTable->setHorizontalHeaderItem(0, __qtablewidgetitem9);
        QTableWidgetItem *__qtablewidgetitem10 = new QTableWidgetItem();
        beadTable->setHorizontalHeaderItem(1, __qtablewidgetitem10);
        QTableWidgetItem *__qtablewidgetitem11 = new QTableWidgetItem();
        beadTable->setHorizontalHeaderItem(2, __qtablewidgetitem11);
        QTableWidgetItem *__qtablewidgetitem12 = new QTableWidgetItem();
        beadTable->setHorizontalHeaderItem(3, __qtablewidgetitem12);
        beadTable->setObjectName(QString::fromUtf8("beadTable"));
        beadTable->setAlternatingRowColors(true);
        beadTable->setColumnCount(4);

        gridLayout->addWidget(beadTable, 0, 5, 1, 1);

        verticalLayout = new QVBoxLayout();
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        saveNPButton = new QPushButton(centralwidget);
        saveNPButton->setObjectName(QString::fromUtf8("saveNPButton"));

        verticalLayout->addWidget(saveNPButton);


        gridLayout->addLayout(verticalLayout, 5, 0, 1, 1);

        horizontalLayout_6 = new QHBoxLayout();
        horizontalLayout_6->setObjectName(QString::fromUtf8("horizontalLayout_6"));
        autoBounds = new QCheckBox(centralwidget);
        autoBounds->setObjectName(QString::fromUtf8("autoBounds"));
        autoBounds->setChecked(true);

        horizontalLayout_6->addWidget(autoBounds);

        label_3 = new QLabel(centralwidget);
        label_3->setObjectName(QString::fromUtf8("label_3"));

        horizontalLayout_6->addWidget(label_3);

        lowerBoundLine = new QLineEdit(centralwidget);
        lowerBoundLine->setObjectName(QString::fromUtf8("lowerBoundLine"));
        lowerBoundLine->setEnabled(false);

        horizontalLayout_6->addWidget(lowerBoundLine);

        label_4 = new QLabel(centralwidget);
        label_4->setObjectName(QString::fromUtf8("label_4"));

        horizontalLayout_6->addWidget(label_4);

        upperBoundLine = new QLineEdit(centralwidget);
        upperBoundLine->setObjectName(QString::fromUtf8("upperBoundLine"));
        upperBoundLine->setEnabled(false);

        horizontalLayout_6->addWidget(upperBoundLine);

        horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_6->addItem(horizontalSpacer_2);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_6->addItem(horizontalSpacer);


        gridLayout->addLayout(horizontalLayout_6, 5, 5, 1, 1);


        horizontalLayout_2->addLayout(gridLayout);

        MainWindow->setCentralWidget(centralwidget);
        menubar = new QMenuBar(MainWindow);
        menubar->setObjectName(QString::fromUtf8("menubar"));
        menubar->setGeometry(QRect(0, 0, 1262, 22));
        menuFile = new QMenu(menubar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        menuHelp = new QMenu(menubar);
        menuHelp->setObjectName(QString::fromUtf8("menuHelp"));
        MainWindow->setMenuBar(menubar);
        statusbar = new QStatusBar(MainWindow);
        statusbar->setObjectName(QString::fromUtf8("statusbar"));
        MainWindow->setStatusBar(statusbar);

        menubar->addAction(menuFile->menuAction());
        menubar->addAction(menuHelp->menuAction());
        menuFile->addAction(actionNew);
        menuFile->addAction(actionLoad);
        menuFile->addAction(actionQuit);
        menuHelp->addAction(actionTips);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QCoreApplication::translate("MainWindow", "NPDesigner", nullptr));
        actionNew->setText(QCoreApplication::translate("MainWindow", "New", nullptr));
        actionQuit->setText(QCoreApplication::translate("MainWindow", "Quit", nullptr));
        actionTips->setText(QCoreApplication::translate("MainWindow", "Tips", nullptr));
        actionLoad->setText(QCoreApplication::translate("MainWindow", "Load", nullptr));
        newBeadTypeButton->setText(QCoreApplication::translate("MainWindow", "Add Bead Type", nullptr));
        QTableWidgetItem *___qtablewidgetitem = beadTypeTable->horizontalHeaderItem(0);
        ___qtablewidgetitem->setText(QCoreApplication::translate("MainWindow", "BeadTypeID", nullptr));
        QTableWidgetItem *___qtablewidgetitem1 = beadTypeTable->horizontalHeaderItem(1);
        ___qtablewidgetitem1->setText(QCoreApplication::translate("MainWindow", "SurfaceDirectory", nullptr));
        QTableWidgetItem *___qtablewidgetitem2 = beadTypeTable->horizontalHeaderItem(2);
        ___qtablewidgetitem2->setText(QCoreApplication::translate("MainWindow", "HamakerFile", nullptr));
        QTableWidgetItem *___qtablewidgetitem3 = beadTypeTable->horizontalHeaderItem(3);
        ___qtablewidgetitem3->setText(QCoreApplication::translate("MainWindow", "Radius [nm]", nullptr));
        QTableWidgetItem *___qtablewidgetitem4 = beadTypeTable->horizontalHeaderItem(4);
        ___qtablewidgetitem4->setText(QCoreApplication::translate("MainWindow", "Zeta [mV]", nullptr));
        QTableWidgetItem *___qtablewidgetitem5 = beadTypeTable->horizontalHeaderItem(5);
        ___qtablewidgetitem5->setText(QCoreApplication::translate("MainWindow", "Surface-factor", nullptr));
        QTableWidgetItem *___qtablewidgetitem6 = beadTypeTable->horizontalHeaderItem(6);
        ___qtablewidgetitem6->setText(QCoreApplication::translate("MainWindow", "Core-factor", nullptr));
        QTableWidgetItem *___qtablewidgetitem7 = beadTypeTable->horizontalHeaderItem(7);
        ___qtablewidgetitem7->setText(QCoreApplication::translate("MainWindow", "LJ-cutoff [nm]", nullptr));
        QTableWidgetItem *___qtablewidgetitem8 = beadTypeTable->horizontalHeaderItem(8);
        ___qtablewidgetitem8->setText(QCoreApplication::translate("MainWindow", "Correction-override", nullptr));
        recenterNP->setText(QCoreApplication::translate("MainWindow", "Re-center", nullptr));
        label->setText(QCoreApplication::translate("MainWindow", "UA directory", nullptr));
        findUADir->setText(QCoreApplication::translate("MainWindow", "Find UA Dir", nullptr));
        newBeadButton->setText(QCoreApplication::translate("MainWindow", "Add NP Bead", nullptr));
        updateTables->setText(QCoreApplication::translate("MainWindow", "Update from tables", nullptr));
        newBrushButton->setText(QCoreApplication::translate("MainWindow", "Add Brush", nullptr));
        newShellButton->setText(QCoreApplication::translate("MainWindow", "Add Shell", nullptr));
        label_2->setText(QCoreApplication::translate("MainWindow", "Orientations", nullptr));
        QTableWidgetItem *___qtablewidgetitem9 = beadTable->horizontalHeaderItem(0);
        ___qtablewidgetitem9->setText(QCoreApplication::translate("MainWindow", "BeadType", nullptr));
        QTableWidgetItem *___qtablewidgetitem10 = beadTable->horizontalHeaderItem(1);
        ___qtablewidgetitem10->setText(QCoreApplication::translate("MainWindow", "x [nm]", nullptr));
        QTableWidgetItem *___qtablewidgetitem11 = beadTable->horizontalHeaderItem(2);
        ___qtablewidgetitem11->setText(QCoreApplication::translate("MainWindow", "y [nm]", nullptr));
        QTableWidgetItem *___qtablewidgetitem12 = beadTable->horizontalHeaderItem(3);
        ___qtablewidgetitem12->setText(QCoreApplication::translate("MainWindow", "z [nm]", nullptr));
        saveNPButton->setText(QCoreApplication::translate("MainWindow", "Save NP(s)", nullptr));
        autoBounds->setText(QCoreApplication::translate("MainWindow", "Auto-set binding radii?", nullptr));
        label_3->setText(QCoreApplication::translate("MainWindow", "R0", nullptr));
        lowerBoundLine->setText(QCoreApplication::translate("MainWindow", "1", nullptr));
        label_4->setText(QCoreApplication::translate("MainWindow", "R1", nullptr));
        upperBoundLine->setText(QCoreApplication::translate("MainWindow", "3", nullptr));
        menuFile->setTitle(QCoreApplication::translate("MainWindow", "File", nullptr));
        menuHelp->setTitle(QCoreApplication::translate("MainWindow", "Help", nullptr));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
