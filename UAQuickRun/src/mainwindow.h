#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QProcess>
#include <QGraphicsScene>
#include "clickablescene.h"

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
    Atom( std::string nameIn, double x, double y, double z){
        atomName = nameIn;
        x0 = x;
        xc = x;
        y0 = y;
        yc = y;
        z0 = z;
        zc = z;
        dAv = 0.01;
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
    double energyData[72][36] = {{0.0}} ;
    std::vector<Atom> atomList;
    int uamBoxHeight = 1;
    int uamBoxWidth = 1;

private slots:
    void on_loadUAMButton_clicked();

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
    void loadMaterials(QString materialFile);
    void calcBeadDistances(bool doRotate);


    void on_colourBoltz_clicked();

    void on_findMinEnergyButton_clicked();

    void on_findBoltzMinButton_clicked();

private:
    Ui::MainWindow *ui;
};
#endif // MAINWINDOW_H
