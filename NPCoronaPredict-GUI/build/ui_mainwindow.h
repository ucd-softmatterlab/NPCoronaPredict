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
#include <QtWidgets/QApplication>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDial>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QGraphicsView>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPlainTextEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSlider>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralwidget;
    QGridLayout *gridLayout_6;
    QGridLayout *gridLayout;
    QTabWidget *tabWidget;
    QWidget *tab;
    QGridLayout *gridLayout_3;
    QVBoxLayout *verticalLayout_2;
    QHBoxLayout *horizontalLayout;
    QLineEdit *uaFolderBox;
    QPushButton *findUAButton;
    QSpacerItem *verticalSpacer_5;
    QHBoxLayout *horizontalLayout_14;
    QHBoxLayout *horizontalLayout_15;
    QCheckBox *npcpModeBox;
    QComboBox *npcpModeOptions;
    QHBoxLayout *horizontalLayout_2;
    QLineEdit *pdbTargetLine;
    QPushButton *pdbTargetButton;
    QSpacerItem *verticalSpacer_4;
    QHBoxLayout *horizontalLayout_4;
    QCheckBox *autoNPBox;
    QGridLayout *gridLayout_2;
    QSpinBox *zetaSpinBox;
    QLabel *label;
    QLabel *label_2;
    QSpinBox *radiusSpinBox;
    QVBoxLayout *verticalLayout;
    QHBoxLayout *horizontalLayout_3;
    QLabel *label_3;
    QComboBox *materialDropdown;
    QHBoxLayout *horizontalLayout_20;
    QComboBox *npTargetShapeOverride;
    QPushButton *loadMaterialButton;
    QHBoxLayout *horizontalLayout_5;
    QLineEdit *npTargetBox;
    QPushButton *npTargetButton;
    QHBoxLayout *horizontalLayout_16;
    QLabel *label_14;
    QDoubleSpinBox *jitterSpinBox;
    QCheckBox *boltzModeCheckBox;
    QSpacerItem *horizontalSpacer_4;
    QSpacerItem *verticalSpacer_3;
    QHBoxLayout *horizontalLayout_6;
    QLineEdit *resultFolderBox;
    QPushButton *resultFolderButton;
    QSpacerItem *verticalSpacer;
    QHBoxLayout *horizontalLayout_7;
    QSpacerItem *horizontalSpacer_2;
    QPushButton *runUAButton;
    QPushButton *cancelRunButton;
    QSpacerItem *horizontalSpacer;
    QPlainTextEdit *uaOutputBox;
    QWidget *tab_2;
    QHBoxLayout *horizontalLayout_12;
    QVBoxLayout *verticalLayout_5;
    QHBoxLayout *horizontalLayout_9;
    QGraphicsView *heatmapView;
    QGraphicsView *pdbView;
    QHBoxLayout *horizontalLayout_13;
    QGraphicsView *heatmapScaleBar;
    QGraphicsView *pdbColourScaleBar;
    QHBoxLayout *horizontalLayout_10;
    QVBoxLayout *verticalLayout_3;
    QLineEdit *loadUAMBox;
    QPushButton *loadUAMButton;
    QVBoxLayout *verticalLayout_4;
    QLineEdit *loadPDBBox;
    QPushButton *loadPDBButton;
    QHBoxLayout *horizontalLayout_11;
    QGridLayout *gridLayout_5;
    QLabel *label_4;
    QLabel *label_9;
    QLabel *label_5;
    QSpinBox *phiInputBox;
    QSpinBox *thetaInputBox;
    QSpinBox *npViewRadius;
    QDial *omegaDial;
    QLabel *label_10;
    QGridLayout *gridLayout_4;
    QLineEdit *energyOutBox;
    QPushButton *findMinEnergyButton;
    QLabel *label_7;
    QLineEdit *simpleOutBox;
    QPushButton *colourEnergy;
    QPushButton *loadBeadmapButton;
    QPushButton *findBoltzMinButton;
    QLineEdit *comDistOutBox;
    QPushButton *colourBoltz;
    QLineEdit *boltzOutBox;
    QLabel *label_6;
    QLabel *label_8;
    QLabel *label_13;
    QHBoxLayout *horizontalLayout_17;
    QCheckBox *showChargeBox;
    QCheckBox *showCOMBox;
    QHBoxLayout *horizontalLayout_18;
    QLabel *label_12;
    QSlider *opacitySlider;
    QHBoxLayout *horizontalLayout_19;
    QLabel *label_15;
    QComboBox *npShapeBox;
    QWidget *tab_3;
    QGridLayout *gridLayout_7;
    QTableWidget *mediumEditTable;
    QHBoxLayout *horizontalLayout_8;
    QLabel *label_11;
    QPushButton *checkStructureButton;
    QSpacerItem *horizontalSpacer_3;
    QPushButton *mediumNewButton;
    QPushButton *mediumLoadButton;
    QPushButton *mediumSaveButton;
    QMenuBar *menubar;
    QStatusBar *statusbar;
    QToolBar *toolBar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(919, 875);
        centralwidget = new QWidget(MainWindow);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        gridLayout_6 = new QGridLayout(centralwidget);
        gridLayout_6->setObjectName(QString::fromUtf8("gridLayout_6"));
        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        tabWidget = new QTabWidget(centralwidget);
        tabWidget->setObjectName(QString::fromUtf8("tabWidget"));
        tab = new QWidget();
        tab->setObjectName(QString::fromUtf8("tab"));
        gridLayout_3 = new QGridLayout(tab);
        gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));
        verticalLayout_2 = new QVBoxLayout();
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        uaFolderBox = new QLineEdit(tab);
        uaFolderBox->setObjectName(QString::fromUtf8("uaFolderBox"));

        horizontalLayout->addWidget(uaFolderBox);

        findUAButton = new QPushButton(tab);
        findUAButton->setObjectName(QString::fromUtf8("findUAButton"));

        horizontalLayout->addWidget(findUAButton);


        verticalLayout_2->addLayout(horizontalLayout);

        verticalSpacer_5 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Fixed);

        verticalLayout_2->addItem(verticalSpacer_5);

        horizontalLayout_14 = new QHBoxLayout();
        horizontalLayout_14->setObjectName(QString::fromUtf8("horizontalLayout_14"));

        verticalLayout_2->addLayout(horizontalLayout_14);

        horizontalLayout_15 = new QHBoxLayout();
        horizontalLayout_15->setObjectName(QString::fromUtf8("horizontalLayout_15"));
        npcpModeBox = new QCheckBox(tab);
        npcpModeBox->setObjectName(QString::fromUtf8("npcpModeBox"));

        horizontalLayout_15->addWidget(npcpModeBox);

        npcpModeOptions = new QComboBox(tab);
        npcpModeOptions->addItem(QString());
        npcpModeOptions->addItem(QString());
        npcpModeOptions->addItem(QString());
        npcpModeOptions->addItem(QString());
        npcpModeOptions->setObjectName(QString::fromUtf8("npcpModeOptions"));
        npcpModeOptions->setEnabled(false);

        horizontalLayout_15->addWidget(npcpModeOptions);


        verticalLayout_2->addLayout(horizontalLayout_15);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        pdbTargetLine = new QLineEdit(tab);
        pdbTargetLine->setObjectName(QString::fromUtf8("pdbTargetLine"));

        horizontalLayout_2->addWidget(pdbTargetLine);

        pdbTargetButton = new QPushButton(tab);
        pdbTargetButton->setObjectName(QString::fromUtf8("pdbTargetButton"));
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(pdbTargetButton->sizePolicy().hasHeightForWidth());
        pdbTargetButton->setSizePolicy(sizePolicy);

        horizontalLayout_2->addWidget(pdbTargetButton);


        verticalLayout_2->addLayout(horizontalLayout_2);

        verticalSpacer_4 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Fixed);

        verticalLayout_2->addItem(verticalSpacer_4);

        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setObjectName(QString::fromUtf8("horizontalLayout_4"));
        autoNPBox = new QCheckBox(tab);
        autoNPBox->setObjectName(QString::fromUtf8("autoNPBox"));
        autoNPBox->setMaximumSize(QSize(500, 16777215));
        autoNPBox->setChecked(true);

        horizontalLayout_4->addWidget(autoNPBox);

        gridLayout_2 = new QGridLayout();
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        zetaSpinBox = new QSpinBox(tab);
        zetaSpinBox->setObjectName(QString::fromUtf8("zetaSpinBox"));
        zetaSpinBox->setMinimum(-50);
        zetaSpinBox->setMaximum(50);

        gridLayout_2->addWidget(zetaSpinBox, 1, 1, 1, 1);

        label = new QLabel(tab);
        label->setObjectName(QString::fromUtf8("label"));

        gridLayout_2->addWidget(label, 0, 0, 1, 1);

        label_2 = new QLabel(tab);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        gridLayout_2->addWidget(label_2, 1, 0, 1, 1);

        radiusSpinBox = new QSpinBox(tab);
        radiusSpinBox->setObjectName(QString::fromUtf8("radiusSpinBox"));
        radiusSpinBox->setMinimum(1);
        radiusSpinBox->setMaximum(500);
        radiusSpinBox->setValue(10);

        gridLayout_2->addWidget(radiusSpinBox, 0, 1, 1, 1);


        horizontalLayout_4->addLayout(gridLayout_2);

        verticalLayout = new QVBoxLayout();
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        label_3 = new QLabel(tab);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        label_3->setMaximumSize(QSize(80, 16777215));

        horizontalLayout_3->addWidget(label_3);

        materialDropdown = new QComboBox(tab);
        materialDropdown->setObjectName(QString::fromUtf8("materialDropdown"));

        horizontalLayout_3->addWidget(materialDropdown);

        horizontalLayout_20 = new QHBoxLayout();
        horizontalLayout_20->setObjectName(QString::fromUtf8("horizontalLayout_20"));
        npTargetShapeOverride = new QComboBox(tab);
        npTargetShapeOverride->addItem(QString());
        npTargetShapeOverride->addItem(QString());
        npTargetShapeOverride->addItem(QString());
        npTargetShapeOverride->addItem(QString());
        npTargetShapeOverride->addItem(QString());
        npTargetShapeOverride->setObjectName(QString::fromUtf8("npTargetShapeOverride"));
        npTargetShapeOverride->setMaximumSize(QSize(150, 16777215));

        horizontalLayout_20->addWidget(npTargetShapeOverride);


        horizontalLayout_3->addLayout(horizontalLayout_20);

        horizontalLayout_3->setStretch(1, 2);
        horizontalLayout_3->setStretch(2, 1);

        verticalLayout->addLayout(horizontalLayout_3);

        loadMaterialButton = new QPushButton(tab);
        loadMaterialButton->setObjectName(QString::fromUtf8("loadMaterialButton"));

        verticalLayout->addWidget(loadMaterialButton);


        horizontalLayout_4->addLayout(verticalLayout);

        horizontalLayout_4->setStretch(1, 1);
        horizontalLayout_4->setStretch(2, 2);

        verticalLayout_2->addLayout(horizontalLayout_4);

        horizontalLayout_5 = new QHBoxLayout();
        horizontalLayout_5->setObjectName(QString::fromUtf8("horizontalLayout_5"));
        npTargetBox = new QLineEdit(tab);
        npTargetBox->setObjectName(QString::fromUtf8("npTargetBox"));
        npTargetBox->setEnabled(false);

        horizontalLayout_5->addWidget(npTargetBox);

        npTargetButton = new QPushButton(tab);
        npTargetButton->setObjectName(QString::fromUtf8("npTargetButton"));
        npTargetButton->setEnabled(false);

        horizontalLayout_5->addWidget(npTargetButton);


        verticalLayout_2->addLayout(horizontalLayout_5);

        horizontalLayout_16 = new QHBoxLayout();
        horizontalLayout_16->setObjectName(QString::fromUtf8("horizontalLayout_16"));
        label_14 = new QLabel(tab);
        label_14->setObjectName(QString::fromUtf8("label_14"));

        horizontalLayout_16->addWidget(label_14);

        jitterSpinBox = new QDoubleSpinBox(tab);
        jitterSpinBox->setObjectName(QString::fromUtf8("jitterSpinBox"));
        jitterSpinBox->setSingleStep(0.010000000000000);

        horizontalLayout_16->addWidget(jitterSpinBox);

        boltzModeCheckBox = new QCheckBox(tab);
        boltzModeCheckBox->setObjectName(QString::fromUtf8("boltzModeCheckBox"));

        horizontalLayout_16->addWidget(boltzModeCheckBox);

        horizontalSpacer_4 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_16->addItem(horizontalSpacer_4);


        verticalLayout_2->addLayout(horizontalLayout_16);

        verticalSpacer_3 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Fixed);

        verticalLayout_2->addItem(verticalSpacer_3);

        horizontalLayout_6 = new QHBoxLayout();
        horizontalLayout_6->setObjectName(QString::fromUtf8("horizontalLayout_6"));
        resultFolderBox = new QLineEdit(tab);
        resultFolderBox->setObjectName(QString::fromUtf8("resultFolderBox"));

        horizontalLayout_6->addWidget(resultFolderBox);

        resultFolderButton = new QPushButton(tab);
        resultFolderButton->setObjectName(QString::fromUtf8("resultFolderButton"));

        horizontalLayout_6->addWidget(resultFolderButton);


        verticalLayout_2->addLayout(horizontalLayout_6);

        verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Fixed);

        verticalLayout_2->addItem(verticalSpacer);

        horizontalLayout_7 = new QHBoxLayout();
        horizontalLayout_7->setObjectName(QString::fromUtf8("horizontalLayout_7"));
        horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_7->addItem(horizontalSpacer_2);

        runUAButton = new QPushButton(tab);
        runUAButton->setObjectName(QString::fromUtf8("runUAButton"));
        QSizePolicy sizePolicy1(QSizePolicy::Expanding, QSizePolicy::Fixed);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(runUAButton->sizePolicy().hasHeightForWidth());
        runUAButton->setSizePolicy(sizePolicy1);

        horizontalLayout_7->addWidget(runUAButton);

        cancelRunButton = new QPushButton(tab);
        cancelRunButton->setObjectName(QString::fromUtf8("cancelRunButton"));
        cancelRunButton->setEnabled(false);

        horizontalLayout_7->addWidget(cancelRunButton);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_7->addItem(horizontalSpacer);


        verticalLayout_2->addLayout(horizontalLayout_7);

        uaOutputBox = new QPlainTextEdit(tab);
        uaOutputBox->setObjectName(QString::fromUtf8("uaOutputBox"));
        uaOutputBox->setEnabled(true);
        uaOutputBox->setReadOnly(true);

        verticalLayout_2->addWidget(uaOutputBox);


        gridLayout_3->addLayout(verticalLayout_2, 0, 0, 1, 1);

        tabWidget->addTab(tab, QString());
        tab_2 = new QWidget();
        tab_2->setObjectName(QString::fromUtf8("tab_2"));
        horizontalLayout_12 = new QHBoxLayout(tab_2);
        horizontalLayout_12->setObjectName(QString::fromUtf8("horizontalLayout_12"));
        verticalLayout_5 = new QVBoxLayout();
        verticalLayout_5->setObjectName(QString::fromUtf8("verticalLayout_5"));
        horizontalLayout_9 = new QHBoxLayout();
        horizontalLayout_9->setObjectName(QString::fromUtf8("horizontalLayout_9"));
        heatmapView = new QGraphicsView(tab_2);
        heatmapView->setObjectName(QString::fromUtf8("heatmapView"));

        horizontalLayout_9->addWidget(heatmapView);

        pdbView = new QGraphicsView(tab_2);
        pdbView->setObjectName(QString::fromUtf8("pdbView"));

        horizontalLayout_9->addWidget(pdbView);


        verticalLayout_5->addLayout(horizontalLayout_9);

        horizontalLayout_13 = new QHBoxLayout();
        horizontalLayout_13->setObjectName(QString::fromUtf8("horizontalLayout_13"));
        heatmapScaleBar = new QGraphicsView(tab_2);
        heatmapScaleBar->setObjectName(QString::fromUtf8("heatmapScaleBar"));
        sizePolicy1.setHeightForWidth(heatmapScaleBar->sizePolicy().hasHeightForWidth());
        heatmapScaleBar->setSizePolicy(sizePolicy1);
        heatmapScaleBar->setMaximumSize(QSize(16777215, 50));

        horizontalLayout_13->addWidget(heatmapScaleBar);

        pdbColourScaleBar = new QGraphicsView(tab_2);
        pdbColourScaleBar->setObjectName(QString::fromUtf8("pdbColourScaleBar"));
        sizePolicy1.setHeightForWidth(pdbColourScaleBar->sizePolicy().hasHeightForWidth());
        pdbColourScaleBar->setSizePolicy(sizePolicy1);
        pdbColourScaleBar->setMaximumSize(QSize(16777215, 50));

        horizontalLayout_13->addWidget(pdbColourScaleBar);


        verticalLayout_5->addLayout(horizontalLayout_13);

        horizontalLayout_10 = new QHBoxLayout();
        horizontalLayout_10->setObjectName(QString::fromUtf8("horizontalLayout_10"));
        verticalLayout_3 = new QVBoxLayout();
        verticalLayout_3->setObjectName(QString::fromUtf8("verticalLayout_3"));
        loadUAMBox = new QLineEdit(tab_2);
        loadUAMBox->setObjectName(QString::fromUtf8("loadUAMBox"));

        verticalLayout_3->addWidget(loadUAMBox);

        loadUAMButton = new QPushButton(tab_2);
        loadUAMButton->setObjectName(QString::fromUtf8("loadUAMButton"));

        verticalLayout_3->addWidget(loadUAMButton);


        horizontalLayout_10->addLayout(verticalLayout_3);

        verticalLayout_4 = new QVBoxLayout();
        verticalLayout_4->setObjectName(QString::fromUtf8("verticalLayout_4"));
        loadPDBBox = new QLineEdit(tab_2);
        loadPDBBox->setObjectName(QString::fromUtf8("loadPDBBox"));

        verticalLayout_4->addWidget(loadPDBBox);

        loadPDBButton = new QPushButton(tab_2);
        loadPDBButton->setObjectName(QString::fromUtf8("loadPDBButton"));

        verticalLayout_4->addWidget(loadPDBButton);


        horizontalLayout_10->addLayout(verticalLayout_4);


        verticalLayout_5->addLayout(horizontalLayout_10);

        horizontalLayout_11 = new QHBoxLayout();
        horizontalLayout_11->setObjectName(QString::fromUtf8("horizontalLayout_11"));
        gridLayout_5 = new QGridLayout();
        gridLayout_5->setObjectName(QString::fromUtf8("gridLayout_5"));
        label_4 = new QLabel(tab_2);
        label_4->setObjectName(QString::fromUtf8("label_4"));

        gridLayout_5->addWidget(label_4, 0, 0, 1, 1);

        label_9 = new QLabel(tab_2);
        label_9->setObjectName(QString::fromUtf8("label_9"));

        gridLayout_5->addWidget(label_9, 2, 0, 1, 1);

        label_5 = new QLabel(tab_2);
        label_5->setObjectName(QString::fromUtf8("label_5"));

        gridLayout_5->addWidget(label_5, 1, 0, 1, 1);

        phiInputBox = new QSpinBox(tab_2);
        phiInputBox->setObjectName(QString::fromUtf8("phiInputBox"));
        phiInputBox->setMaximum(360);

        gridLayout_5->addWidget(phiInputBox, 0, 1, 1, 1);

        thetaInputBox = new QSpinBox(tab_2);
        thetaInputBox->setObjectName(QString::fromUtf8("thetaInputBox"));
        thetaInputBox->setMaximum(180);

        gridLayout_5->addWidget(thetaInputBox, 1, 1, 1, 1);

        npViewRadius = new QSpinBox(tab_2);
        npViewRadius->setObjectName(QString::fromUtf8("npViewRadius"));
        npViewRadius->setMinimum(1);
        npViewRadius->setMaximum(999);
        npViewRadius->setValue(5);

        gridLayout_5->addWidget(npViewRadius, 2, 1, 1, 1);

        omegaDial = new QDial(tab_2);
        omegaDial->setObjectName(QString::fromUtf8("omegaDial"));
        omegaDial->setMaximum(360);
        omegaDial->setWrapping(true);
        omegaDial->setNotchesVisible(true);

        gridLayout_5->addWidget(omegaDial, 3, 1, 1, 1);

        label_10 = new QLabel(tab_2);
        label_10->setObjectName(QString::fromUtf8("label_10"));

        gridLayout_5->addWidget(label_10, 3, 0, 1, 1);


        horizontalLayout_11->addLayout(gridLayout_5);

        gridLayout_4 = new QGridLayout();
        gridLayout_4->setObjectName(QString::fromUtf8("gridLayout_4"));
        energyOutBox = new QLineEdit(tab_2);
        energyOutBox->setObjectName(QString::fromUtf8("energyOutBox"));
        energyOutBox->setAcceptDrops(false);
        energyOutBox->setReadOnly(true);

        gridLayout_4->addWidget(energyOutBox, 0, 3, 1, 1);

        findMinEnergyButton = new QPushButton(tab_2);
        findMinEnergyButton->setObjectName(QString::fromUtf8("findMinEnergyButton"));

        gridLayout_4->addWidget(findMinEnergyButton, 3, 1, 1, 1);

        label_7 = new QLabel(tab_2);
        label_7->setObjectName(QString::fromUtf8("label_7"));

        gridLayout_4->addWidget(label_7, 3, 2, 1, 1);

        simpleOutBox = new QLineEdit(tab_2);
        simpleOutBox->setObjectName(QString::fromUtf8("simpleOutBox"));
        simpleOutBox->setAcceptDrops(false);
        simpleOutBox->setReadOnly(true);

        gridLayout_4->addWidget(simpleOutBox, 3, 3, 1, 1);

        colourEnergy = new QPushButton(tab_2);
        colourEnergy->setObjectName(QString::fromUtf8("colourEnergy"));

        gridLayout_4->addWidget(colourEnergy, 1, 1, 1, 1);

        loadBeadmapButton = new QPushButton(tab_2);
        loadBeadmapButton->setObjectName(QString::fromUtf8("loadBeadmapButton"));

        gridLayout_4->addWidget(loadBeadmapButton, 0, 4, 1, 1);

        findBoltzMinButton = new QPushButton(tab_2);
        findBoltzMinButton->setObjectName(QString::fromUtf8("findBoltzMinButton"));

        gridLayout_4->addWidget(findBoltzMinButton, 4, 1, 1, 1);

        comDistOutBox = new QLineEdit(tab_2);
        comDistOutBox->setObjectName(QString::fromUtf8("comDistOutBox"));
        comDistOutBox->setReadOnly(true);

        gridLayout_4->addWidget(comDistOutBox, 1, 3, 1, 1);

        colourBoltz = new QPushButton(tab_2);
        colourBoltz->setObjectName(QString::fromUtf8("colourBoltz"));

        gridLayout_4->addWidget(colourBoltz, 0, 1, 1, 1);

        boltzOutBox = new QLineEdit(tab_2);
        boltzOutBox->setObjectName(QString::fromUtf8("boltzOutBox"));
        boltzOutBox->setAcceptDrops(false);
        boltzOutBox->setReadOnly(true);

        gridLayout_4->addWidget(boltzOutBox, 4, 3, 1, 1);

        label_6 = new QLabel(tab_2);
        label_6->setObjectName(QString::fromUtf8("label_6"));

        gridLayout_4->addWidget(label_6, 0, 2, 1, 1);

        label_8 = new QLabel(tab_2);
        label_8->setObjectName(QString::fromUtf8("label_8"));

        gridLayout_4->addWidget(label_8, 4, 2, 1, 1);

        label_13 = new QLabel(tab_2);
        label_13->setObjectName(QString::fromUtf8("label_13"));

        gridLayout_4->addWidget(label_13, 1, 2, 1, 1);

        horizontalLayout_17 = new QHBoxLayout();
        horizontalLayout_17->setObjectName(QString::fromUtf8("horizontalLayout_17"));
        showChargeBox = new QCheckBox(tab_2);
        showChargeBox->setObjectName(QString::fromUtf8("showChargeBox"));

        horizontalLayout_17->addWidget(showChargeBox);

        showCOMBox = new QCheckBox(tab_2);
        showCOMBox->setObjectName(QString::fromUtf8("showCOMBox"));
        showCOMBox->setChecked(true);

        horizontalLayout_17->addWidget(showCOMBox);


        gridLayout_4->addLayout(horizontalLayout_17, 1, 4, 1, 1);

        horizontalLayout_18 = new QHBoxLayout();
        horizontalLayout_18->setObjectName(QString::fromUtf8("horizontalLayout_18"));
        label_12 = new QLabel(tab_2);
        label_12->setObjectName(QString::fromUtf8("label_12"));

        horizontalLayout_18->addWidget(label_12);

        opacitySlider = new QSlider(tab_2);
        opacitySlider->setObjectName(QString::fromUtf8("opacitySlider"));
        opacitySlider->setMaximum(255);
        opacitySlider->setSingleStep(1);
        opacitySlider->setValue(255);
        opacitySlider->setOrientation(Qt::Horizontal);

        horizontalLayout_18->addWidget(opacitySlider);


        gridLayout_4->addLayout(horizontalLayout_18, 3, 4, 1, 1);

        horizontalLayout_19 = new QHBoxLayout();
        horizontalLayout_19->setObjectName(QString::fromUtf8("horizontalLayout_19"));
        label_15 = new QLabel(tab_2);
        label_15->setObjectName(QString::fromUtf8("label_15"));

        horizontalLayout_19->addWidget(label_15);

        npShapeBox = new QComboBox(tab_2);
        npShapeBox->addItem(QString());
        npShapeBox->addItem(QString());
        npShapeBox->addItem(QString());
        npShapeBox->setObjectName(QString::fromUtf8("npShapeBox"));

        horizontalLayout_19->addWidget(npShapeBox);


        gridLayout_4->addLayout(horizontalLayout_19, 4, 4, 1, 1);


        horizontalLayout_11->addLayout(gridLayout_4);


        verticalLayout_5->addLayout(horizontalLayout_11);


        horizontalLayout_12->addLayout(verticalLayout_5);

        tabWidget->addTab(tab_2, QString());
        tab_3 = new QWidget();
        tab_3->setObjectName(QString::fromUtf8("tab_3"));
        gridLayout_7 = new QGridLayout(tab_3);
        gridLayout_7->setObjectName(QString::fromUtf8("gridLayout_7"));
        mediumEditTable = new QTableWidget(tab_3);
        if (mediumEditTable->columnCount() < 2)
            mediumEditTable->setColumnCount(2);
        QTableWidgetItem *__qtablewidgetitem = new QTableWidgetItem();
        mediumEditTable->setHorizontalHeaderItem(0, __qtablewidgetitem);
        QTableWidgetItem *__qtablewidgetitem1 = new QTableWidgetItem();
        mediumEditTable->setHorizontalHeaderItem(1, __qtablewidgetitem1);
        mediumEditTable->setObjectName(QString::fromUtf8("mediumEditTable"));
        mediumEditTable->setContextMenuPolicy(Qt::CustomContextMenu);
        mediumEditTable->horizontalHeader()->setDefaultSectionSize(400);
        mediumEditTable->horizontalHeader()->setStretchLastSection(true);

        gridLayout_7->addWidget(mediumEditTable, 0, 0, 1, 1);

        horizontalLayout_8 = new QHBoxLayout();
        horizontalLayout_8->setObjectName(QString::fromUtf8("horizontalLayout_8"));
        label_11 = new QLabel(tab_3);
        label_11->setObjectName(QString::fromUtf8("label_11"));

        horizontalLayout_8->addWidget(label_11);

        checkStructureButton = new QPushButton(tab_3);
        checkStructureButton->setObjectName(QString::fromUtf8("checkStructureButton"));

        horizontalLayout_8->addWidget(checkStructureButton);

        horizontalSpacer_3 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_8->addItem(horizontalSpacer_3);

        mediumNewButton = new QPushButton(tab_3);
        mediumNewButton->setObjectName(QString::fromUtf8("mediumNewButton"));

        horizontalLayout_8->addWidget(mediumNewButton);

        mediumLoadButton = new QPushButton(tab_3);
        mediumLoadButton->setObjectName(QString::fromUtf8("mediumLoadButton"));

        horizontalLayout_8->addWidget(mediumLoadButton);

        mediumSaveButton = new QPushButton(tab_3);
        mediumSaveButton->setObjectName(QString::fromUtf8("mediumSaveButton"));

        horizontalLayout_8->addWidget(mediumSaveButton);


        gridLayout_7->addLayout(horizontalLayout_8, 1, 0, 1, 1);

        tabWidget->addTab(tab_3, QString());

        gridLayout->addWidget(tabWidget, 1, 0, 1, 1);


        gridLayout_6->addLayout(gridLayout, 0, 0, 1, 1);

        MainWindow->setCentralWidget(centralwidget);
        menubar = new QMenuBar(MainWindow);
        menubar->setObjectName(QString::fromUtf8("menubar"));
        menubar->setGeometry(QRect(0, 0, 919, 22));
        MainWindow->setMenuBar(menubar);
        statusbar = new QStatusBar(MainWindow);
        statusbar->setObjectName(QString::fromUtf8("statusbar"));
        MainWindow->setStatusBar(statusbar);
        toolBar = new QToolBar(MainWindow);
        toolBar->setObjectName(QString::fromUtf8("toolBar"));
        MainWindow->addToolBar(Qt::TopToolBarArea, toolBar);

        retranslateUi(MainWindow);

        tabWidget->setCurrentIndex(0);


        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QCoreApplication::translate("MainWindow", "NPCoronaPredict-GUI", nullptr));
        uaFolderBox->setText(QString());
        findUAButton->setText(QCoreApplication::translate("MainWindow", "Find UA Folder", nullptr));
        npcpModeBox->setText(QCoreApplication::translate("MainWindow", "Advanced mode (needs .csv medium file )", nullptr));
        npcpModeOptions->setItemText(0, QCoreApplication::translate("MainWindow", "Prepare only", nullptr));
        npcpModeOptions->setItemText(1, QCoreApplication::translate("MainWindow", "Prepare + UA", nullptr));
        npcpModeOptions->setItemText(2, QCoreApplication::translate("MainWindow", "Prepare + UA + BCP", nullptr));
        npcpModeOptions->setItemText(3, QCoreApplication::translate("MainWindow", "Prepare + UA + BCP + KMC", nullptr));

        pdbTargetButton->setText(QCoreApplication::translate("MainWindow", "PDB Target", nullptr));
        autoNPBox->setText(QCoreApplication::translate("MainWindow", "Auto-NP?", nullptr));
        label->setText(QCoreApplication::translate("MainWindow", "Radius [nm]", nullptr));
        label_2->setText(QCoreApplication::translate("MainWindow", "Zeta [mV]", nullptr));
        label_3->setText(QCoreApplication::translate("MainWindow", "Material", nullptr));
        npTargetShapeOverride->setItemText(0, QCoreApplication::translate("MainWindow", "Sphere", nullptr));
        npTargetShapeOverride->setItemText(1, QCoreApplication::translate("MainWindow", "Cylinder", nullptr));
        npTargetShapeOverride->setItemText(2, QCoreApplication::translate("MainWindow", "Cube", nullptr));
        npTargetShapeOverride->setItemText(3, QCoreApplication::translate("MainWindow", "SWCNT", nullptr));
        npTargetShapeOverride->setItemText(4, QCoreApplication::translate("MainWindow", "MWCNT", nullptr));

        loadMaterialButton->setText(QCoreApplication::translate("MainWindow", "Load Material Set", nullptr));
        npTargetButton->setText(QCoreApplication::translate("MainWindow", "NP Target", nullptr));
        label_14->setText(QCoreApplication::translate("MainWindow", "Jitter magnitude [nm]", nullptr));
        boltzModeCheckBox->setText(QCoreApplication::translate("MainWindow", "Local Boltzmann averaging?", nullptr));
        resultFolderButton->setText(QCoreApplication::translate("MainWindow", "Result folder", nullptr));
        runUAButton->setText(QCoreApplication::translate("MainWindow", "Run UA", nullptr));
        cancelRunButton->setText(QCoreApplication::translate("MainWindow", "Cancel run", nullptr));
        tabWidget->setTabText(tabWidget->indexOf(tab), QCoreApplication::translate("MainWindow", "Run", nullptr));
        loadUAMButton->setText(QCoreApplication::translate("MainWindow", "Load .uam", nullptr));
        loadPDBButton->setText(QCoreApplication::translate("MainWindow", "Load .pdb", nullptr));
        label_4->setText(QCoreApplication::translate("MainWindow", "phi", nullptr));
        label_9->setText(QCoreApplication::translate("MainWindow", "NP Radius", nullptr));
        label_5->setText(QCoreApplication::translate("MainWindow", "theta", nullptr));
        label_10->setText(QCoreApplication::translate("MainWindow", "Rotate view", nullptr));
        findMinEnergyButton->setText(QCoreApplication::translate("MainWindow", "Set to min. energy", nullptr));
        label_7->setText(QCoreApplication::translate("MainWindow", "EAds(simple)", nullptr));
        colourEnergy->setText(QCoreApplication::translate("MainWindow", "Colour by energy", nullptr));
        loadBeadmapButton->setText(QCoreApplication::translate("MainWindow", "Reload CG bead data", nullptr));
        findBoltzMinButton->setText(QCoreApplication::translate("MainWindow", "Set to Boltz minimum", nullptr));
        colourBoltz->setText(QCoreApplication::translate("MainWindow", "Colour by distance", nullptr));
        label_6->setText(QCoreApplication::translate("MainWindow", "EAds(phi,theta)", nullptr));
        label_8->setText(QCoreApplication::translate("MainWindow", "EAds(Boltz)", nullptr));
        label_13->setText(QCoreApplication::translate("MainWindow", "COM-dist", nullptr));
        showChargeBox->setText(QCoreApplication::translate("MainWindow", "Show charges", nullptr));
        showCOMBox->setText(QCoreApplication::translate("MainWindow", "Show COMs", nullptr));
        label_12->setText(QCoreApplication::translate("MainWindow", "Opacity", nullptr));
        label_15->setText(QCoreApplication::translate("MainWindow", "NPShape", nullptr));
        npShapeBox->setItemText(0, QCoreApplication::translate("MainWindow", "Sphere", nullptr));
        npShapeBox->setItemText(1, QCoreApplication::translate("MainWindow", "Cylinder", nullptr));
        npShapeBox->setItemText(2, QCoreApplication::translate("MainWindow", "Cube", nullptr));

        tabWidget->setTabText(tabWidget->indexOf(tab_2), QCoreApplication::translate("MainWindow", "View UA Results", nullptr));
        QTableWidgetItem *___qtablewidgetitem = mediumEditTable->horizontalHeaderItem(0);
        ___qtablewidgetitem->setText(QCoreApplication::translate("MainWindow", "MoleculeID", nullptr));
        QTableWidgetItem *___qtablewidgetitem1 = mediumEditTable->horizontalHeaderItem(1);
        ___qtablewidgetitem1->setText(QCoreApplication::translate("MainWindow", "Concentration [mol/L]", nullptr));
        label_11->setText(QCoreApplication::translate("MainWindow", "File locations: all_proteins/*.pdb", nullptr));
        checkStructureButton->setText(QCoreApplication::translate("MainWindow", "Check Structures", nullptr));
        mediumNewButton->setText(QCoreApplication::translate("MainWindow", "New", nullptr));
        mediumLoadButton->setText(QCoreApplication::translate("MainWindow", "Load", nullptr));
        mediumSaveButton->setText(QCoreApplication::translate("MainWindow", "Save", nullptr));
        tabWidget->setTabText(tabWidget->indexOf(tab_3), QCoreApplication::translate("MainWindow", "Molecule List Editor", nullptr));
        toolBar->setWindowTitle(QCoreApplication::translate("MainWindow", "toolBar", nullptr));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
