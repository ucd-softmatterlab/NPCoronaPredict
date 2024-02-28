#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QProcess>
#include <QMenu>
#include <QGraphicsScene>
#include "clickablescene.h"
#include <QNetworkAccessManager>
#include <QNetworkReply>

class MaterialType{
public:
    QString name;
    QString surfaceDir;
    QString hamakerFile;
    MaterialType( QString nameIn, QString surfaceDirIn, QString hamakerFileIn ){
        name = nameIn;
        surfaceDir = surfaceDirIn;
        hamakerFile = hamakerFileIn;
    }
};

class Atom{
public:
    std::string atomName;
    double x0;
    double y0;
    double z0;
    double xc;
    double yc;
    double zc;

    double xc2;
    double yc2;
    double zc2;

    double dAvn;
     double dAvd;
   double dAv;
   double colourParam;
   double orientationPhi;
   double orientationTheta;
   bool isShrinkWrap = false;
   double radius = 0.5;
    Atom( std::string nameIn, double x, double y, double z, bool swrap = false, double radiusIn = 0.5){
        atomName = nameIn;
        x0 = x;
        xc = x;
        y0 = y;
        yc = y;
        z0 = z;
        zc = z;
        dAv = 0.01;
        colourParam = 0.1;
        orientationPhi=0;
        orientationTheta=0;
        isShrinkWrap = swrap;
        radius = radiusIn;

    }
};


QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    QString uaGlobalPath;
    std::vector<MaterialType> materialTypes;
    QProcess processHandle;
    QPixmap heatmapPixmap;
    ClickableScene scene;
    QGraphicsScene pdbScene;
    QGraphicsScene heatmapBarScene;
    QPixmap heatmapScalePixmap;

    QGraphicsScene pdbBarScene;
    QPixmap pdbScalePixmap;

    QMenu tableMenu;

    QNetworkAccessManager networkManager;

    QString lastUAMPath = "";

    double energyData[72][36] = {{0.0}} ;
    double comData[72][36] = {{-1.0}} ;
    std::vector<Atom> atomList;
    std::vector<Atom> shrinkwrapList;
    int uamBoxHeight = 1;
    int uamBoxWidth = 1;
    bool showEnergyShrinkWrap = false;
    bool blockMediumColouring = false;

private slots:
    void on_loadUAMButton_clicked();
    void updateHeatmapPlot();
    void on_findUAButton_clicked();

    void on_resultFolderButton_clicked();

    void on_loadMaterialButton_clicked();

    void on_pdbTargetButton_clicked();

    void on_runUAButton_clicked();
    void updateUABox();
    void uaDoneAlert(int exitCode, QProcess::ExitStatus exitStatus);
    void updateEnergyBox();
    void on_phiInputBox_valueChanged(int arg1);
    void getSceneMouseClick(QPointF scenePosLoc);
    void on_thetaInputBox_valueChanged(int arg1);
    int phiToIndex(double phiVal, double deltaVal);

    void on_loadPDBButton_clicked();
    void updateMoleculeBox();

    void on_radiusSpinBox_valueChanged(int arg1);

    void on_npViewRadius_valueChanged(int arg1);

    void on_omegaDial_valueChanged(int value);
    void checkForMaterials();
    void loadMaterials(QString materialFile, bool doAlert);
    void calcBeadDistances(bool doRotate);
    void calcBeadBoltzDistances(bool doRotate);


    void on_colourBoltz_clicked();

    void on_findMinEnergyButton_clicked();

    void on_findBoltzMinButton_clicked();

    void on_colourEnergy_clicked();

    void on_autoNPBox_stateChanged(int arg1);

    void on_npTargetButton_clicked();

    void on_npcpModeBox_stateChanged(int arg1);

    void on_mediumEditTable_customContextMenuRequested(const QPoint &pos);

    void addMoleculeToMedium();
    void removeMoleculeFromMedium();

    void on_mediumNewButton_clicked();

    void on_mediumSaveButton_clicked();

    void on_mediumLoadButton_clicked();

    void on_cancelRunButton_clicked();

    void on_checkStructureButton_clicked();
    void colourStructures();
    void colourStructureRow(int row);

    void startDownload(QUrl targetURL);
    void downloadReplyFinished( QNetworkReply *reply );

    void on_mediumEditTable_cellChanged(int row, int column);

    void on_mediumEditTable_currentCellChanged(int currentRow, int currentColumn, int previousRow, int previousColumn);

    void on_showCOMBox_stateChanged(int arg1);

    void on_opacitySlider_sliderMoved(int position);

    void on_opacitySlider_valueChanged(int value);

private:
    Ui::MainWindow *ui;
};
#endif // MAINWINDOW_H
