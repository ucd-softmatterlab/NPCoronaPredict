#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QGraphicsScene>
#include "addbeadtype.h"
#include "addbead.h"
#include "addbrush.h"
#include "addshell.h"
#include "tipswindow.h"
#include <random>

class BeadType{
public:
    int beadID;
    std::string hamakerFile;
    std::string surfaceDir;
    float radius;
    float surfacePotential;
    float surfFactor;
    float coreFactor;
    float ljCutoffVal;
    int correctionOverride;


    BeadType(int beadIDIn, std::string hamakerFileIn, std::string surfaceDirIn, float   radiusIn, float   surfacePotentialIn, float surfFactorIn, float coreFactorIn, float ljCutoffIn, int correctionOverrideIn){
        beadID = beadIDIn;
        hamakerFile = hamakerFileIn;
        surfaceDir = surfaceDirIn;
        radius = radiusIn;
        surfacePotential = surfacePotentialIn;
        surfFactor = surfFactorIn;
        coreFactor = coreFactorIn;
        ljCutoffVal = ljCutoffIn;
        correctionOverride = correctionOverrideIn;
    }
};

class NPBead{
public:
    int beadTypeID;
    float x;
    float y;
    float z;
    NPBead( int beadTypeIDIn, float xIn, float yIn, double zIn ){
        beadTypeID = beadTypeIDIn;
        x = xIn;
        y = yIn;
        z = zIn;
    }
};

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


QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    std::vector<BeadType> beadTypes;
    std::vector<NPBead> npBeads;
    std::vector<MaterialType> materialTypes;
    QString uaGlobalPath;
    std::default_random_engine generator;
    QGraphicsScene scene;
    double npInnerRadius;
    double npOuterRadius;


private slots:
    void on_newBeadTypeButton_clicked();
    void on_newBeadButton_clicked();
    void on_newShellButton_clicked();
    void on_newBrushButton_clicked();
    void recieveNewNPBead( double x, double y, double z, int beadID, bool doUpdate);
    void recieveNewNPBeadType( std::string hamakerFileIn, std::string surfaceDirIn, float   radiusIn, float   surfacePotentialIn, float surfFactorIn, float coreFactorIn, float ljCutoffIn, int correctionOverrideIn);
    void recieveNewShell(std::string hamakerFileIn, std::string surfaceDirIn, float   innerRadius, float outerRadius, float   surfacePotentialIn,  double ljCutoff);
    void recieveNewBrush(double brushOccupancy,double brushRadialDist, int beadTypeID, bool forceAttach);
    void updateBeadTypeTable();
    void updateBeadTable();
    double updateM(double m, int numPoints);
    double updateTheta(double theta, double m, double j);
    void on_actionQuit_triggered();

    void on_findUADir_clicked();

    void on_saveNPButton_clicked();

    void on_actionTips_triggered();
    void updateGraphicsWindow();

    void on_updateTables_clicked();

    void on_beadTypeTable_cellChanged(int row, int column);

    void on_recenterNP_clicked();
    void updateBindingRadii();

    void on_autoBounds_stateChanged(int arg1);

    void on_lowerBoundLine_textEdited(const QString &arg1);

    void on_upperBoundLine_textEdited(const QString &arg1);

    void on_actionNew_triggered();

    void on_actionLoad_triggered();

    void on_pushButton_clicked();

    void on_loadMaterialSet_clicked();

private:
    Ui::MainWindow *ui;
    AddBead *addBead;
    AddBeadType *addBeadType;
    AddBrush *addBrush;
    AddShell *addShell;
    tipsWindow *tipsWindowI;
    void saveNP(QString filename, bool doRotate, bool isPDBFile);
};
#endif // MAINWINDOW_H
