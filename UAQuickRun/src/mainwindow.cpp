#include "mainwindow.h"
#include "ui_mainwindow.h"


#include <vector>
#include <QtDebug>
#include <QFile>
#include <QFileDialog>
#include <QDir>
#include <QComboBox>
#include <random>
#include <QMessageBox>
#include <QGraphicsTextItem>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/algorithm/string.hpp>
#include <QString>
#include <QProcess>


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);




    QString currentPath = QDir::currentPath() ;
    currentPath = QDir::cleanPath(currentPath+"/../");

    this->uaGlobalPath = currentPath;

    if(!QFile::exists(currentPath+"/UnitedAtom")){
        //alert - no UA found

        QMessageBox msgBox;
        msgBox.setText("Warning: UnitedAtom executable not found in the current location. Please check it exists or find the correct folder.");
           msgBox.exec();
           this->statusBar()->showMessage("UA not yet found");
    }
    else{
        this->statusBar()->showMessage("Ready");
    }


    this->findChild<QLineEdit *>("uaFolderBox")->setText(this->uaGlobalPath);
    this->findChild<QLineEdit *>("resultFolderBox")->setText(   QDir::cleanPath(  this->uaGlobalPath   +"/uaquickrunresults"));
    checkForMaterials();
    heatmapPixmap = QPixmap(360,180);
    heatmapPixmap.fill();
    heatmapScalePixmap = QPixmap(400,50);
    pdbScalePixmap = QPixmap(400,50);
    connect(&scene, &ClickableScene::sendMouseClickPos, this, &MainWindow::getSceneMouseClick);
    connect(&processHandle, SIGNAL ( readyReadStandardOutput() )  , this,  SLOT( updateUABox() )        );
    connect(&processHandle, SIGNAL(finished(int, QProcess::ExitStatus) ), this, SLOT(uaDoneAlert(int, QProcess::ExitStatus)));


    //this->findChild<QLabel *>("heatmapPlotLabel")->setPixmap(this->heatmapPixmap);
    //this->materialTypes.emplace_back( MaterialType( "Custom"   , "",""     ) ) ;
    //addBeadType->findChild<QComboBox *>("materialTypeBox")->addItem( "Custom" ,   QList<QVariant>() <<  QString::fromStdString("") <<  QString::fromStdString("")) ;
    //addShell->findChild<QComboBox *>("materialTypeBox")->addItem( "Custom" ,   QList<QVariant>() <<  QString::fromStdString("") <<  QString::fromStdString("")) ;
    //this->findChild<QComboBox *>("materialDropdown")->addItem( "Custom" ,   QList<QVariant>() <<  QString::fromStdString("") <<  QString::fromStdString("")) ;

}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_loadUAMButton_clicked()
{
    QGraphicsView* heatmapBox = this->findChild<QGraphicsView *>("heatmapView") ;
    int gbWidth = heatmapBox->width() ;
    int gbHeight = heatmapBox->height() ;

    QGraphicsView* heatmapScaleBox = this->findChild<QGraphicsView *>("heatmapScaleBar") ;
    int scaleBarWidth = heatmapScaleBox->width() ;
    int scaleBarHeight = heatmapScaleBox->height() ;


    uamBoxHeight = gbHeight;
    uamBoxWidth = gbWidth;
    double heatmapElemWidth = gbWidth/72;
    double heatmapElemHeight = gbHeight/36;

    double angleDelta = 5;
    double piVal = 3.14159265;
    double simpleAverageNum = 0;
    double simpleAverageDenom = 0.00001;
    double boltzAverageNum = 0;
    double boltzAverageDenom = 0.000001;
   double kbtVal = 1;
   double minPhiVal = 1;
   double minThetaVal = 1;

    QString targetUAMFile = QFileDialog::getOpenFileName(this, tr("UAM File"),  this->findChild<QLineEdit *>("resultFolderBox")->text(),  tr("UA-Output (*.uam)"));
    this->findChild<QLineEdit *>("loadUAMBox")->setText(targetUAMFile);
    //load in the .uam file
    if(targetUAMFile!=""){

        QFile file(targetUAMFile);
        if(file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QTextStream fileIn(&file);

        this->statusBar()->showMessage("Loading UAM");
        double minEnergy = 0;
        double maxEnergy = 1;
        //load in the data
        while(!fileIn.atEnd()){
            std::string lineIn = fileIn.readLine().toStdString();
            if(lineIn.substr(0,1) !="#" && lineIn.length() > 5) {

                std::vector<std::string> results;
                 boost::split(results, lineIn, [](char c){return c == ' ';}  , boost::token_compress_on   );
                 double thetaVal = std::stod(  results[1] );
                double phiVal = std::stod(  results[0] );
                double energyVal = std::stod(  results[2] );
       int  thetaIndex = (int)thetaVal/angleDelta;
       int phiIndex = (int)phiVal/angleDelta;
     //  qDebug() << "phi " << QString::number(phiVal) << " to index " << QString::number(phiIndex) << "\n";
       double thetaRad = (thetaVal + 2.5) * piVal/180.0;
       double sinTheta = sin(thetaRad);
       simpleAverageNum += sinTheta*energyVal ;
        boltzAverageNum += sinTheta*energyVal*exp( - energyVal/kbtVal) ;
        simpleAverageDenom += sinTheta;
        boltzAverageDenom += sinTheta*exp(-energyVal/kbtVal);
        energyData[phiIndex][thetaIndex] = energyVal;

        if(energyVal < minEnergy){
            minEnergy = energyVal;
            minPhiVal = phiVal + 2.5;
            minThetaVal = thetaVal + 2.5;
        }
        //minEnergy = std::min(minEnergy,energyVal);
        maxEnergy = std::max(maxEnergy, energyVal);
            }
        }
        heatmapPixmap.fill();
        QPainter *paint = new QPainter(&heatmapPixmap);
        if(abs( minEnergy - maxEnergy) < 0.0001 ){
            minEnergy = maxEnergy - 0.01;
        }
        for(int i = 0; i<72; i++){
            for(int j = 0; j<36;j++){
                double energyVal = energyData[i][j];

                double x1 = i*5 ;
                //double y1=j*5; //this plots theta = 0 at the top
                double y1=(35-j)*5;

                double eScale = (energyVal - maxEnergy)/(  minEnergy - maxEnergy );
                //qDebug() << QString::number(eScale) << " from " << QString::number(energyVal) << "\n";
                paint->setPen(QColor(255*eScale,0,0));
                 paint->setBrush(QColor(255*eScale,0,0));
                paint->drawRect(x1, y1, 5, 5  );


            }
        }
        delete paint;
        heatmapScalePixmap.fill();
        QPainter *scalePaint = new QPainter(&heatmapScalePixmap);
        QFont font = scalePaint->font();
        //font.setPixelSize(10 );
        scalePaint->setFont(font);
        int barWidth = 400/10;
        for(int i =0; i<10;i++){
            double energyVal = minEnergy - i*(minEnergy-maxEnergy)/10;
            double eScale = (energyVal - maxEnergy)/(  minEnergy - maxEnergy );
            //qDebug() << QString::number(eScale) << " from " << QString::number(energyVal) << "\n";
            double y1 = 0;
            double x1 = i*barWidth;
            scalePaint->setPen(QColor(255*eScale,0,0));
             scalePaint->setBrush(QColor(255*eScale,0,0));
            scalePaint->drawRect(x1, y1, barWidth, 20  );
             scalePaint->setPen(QColor(0,0,0));
             //scalePaint->drawText( x1, y1+20, 10, 20, Qt::AlignLeft,  QString::number(energyVal, 'f',1)) ;
             if(i%2 == 0){

                 if(energyVal < -99){
                      scalePaint->drawText( x1, y1+40,   QString::number((int)round(energyVal))) ;
                 }
                 else{
             scalePaint->drawText( x1, y1+40,   QString::number(energyVal, 'f',1)) ;
                 }
             }
        }
        delete scalePaint;


        //save the pixmap out to the label
       // "heatmapPlotLabel" ;


      this->findChild<QSpinBox *>("phiInputBox")->setValue(minPhiVal);
        this->findChild<QSpinBox *>("thetaInputBox")->setValue(minThetaVal);

      //update the boxes
        this->findChild<QLineEdit *>("boltzOutBox")->setText( QString::number(boltzAverageNum/boltzAverageDenom)) ;
        this->findChild<QLineEdit *>("simpleOutBox")->setText( QString::number(simpleAverageNum/simpleAverageDenom)) ;
        updateEnergyBox();


        //update graphics

 // scene.addPixmap(heatmapPixmap );
         scene.addPixmap(heatmapPixmap.scaled(gbWidth*0.95,gbHeight*0.95));
       // this->findChild<QLabel *>("heatmapPlotLabel")->setPixmap(this->heatmapPixmap);
       heatmapBox->setScene(&scene);

       heatmapBarScene.addPixmap(heatmapScalePixmap.scaled(scaleBarWidth*0.95,scaleBarHeight*0.95)   );
       heatmapScaleBox->setScene(&heatmapBarScene);

        //this->findChild<QGraphicsView *>("heatmapView")->addPixmap(this->heatmapPixmap);
      this->statusBar()->showMessage("UAM load complete");

    }
}


}

int MainWindow::phiToIndex(double phiVal, double deltaVal = 5.0){
    int index = (int)floor( (phiVal - deltaVal/2.0) / deltaVal) ;
    if(index < 0){
        index = index+72;
    }
    index = index % 72;
    return index;
}

void MainWindow::updateEnergyBox(){
    double targetPhiVal =( this->findChild<QSpinBox *>("phiInputBox")->value()   );
    double targetThetaVal = this->findChild<QSpinBox *>("thetaInputBox")->value();
    double deltaVal = 5.0;

    /*
    double targetPhiLower = deltaVal*floor(( targetPhiVal - deltaVal/2.0)/deltaVal) + deltaVal/2.0;
    double targetPhiUpper = targetPhiLower + deltaVal;

    int phiLowerIndex = phiToIndex( targetPhiLower);
    int phiUpperIndex = phiToIndex(targetPhiUpper) ;
        int thetaLowerIndex = (int)floor( (targetThetaVal-deltaVal/2.0)/ deltaVal) ;
   double foundEnergy = 0.0;
     */

     //find the two bracketing values for phi, taking into account wrapping of phi
    /*
    if(targetPhiVal <= deltaVal){
     double phiLower = deltaVal/2.0 + floor( (targetPhiVal-deltaVal/2.0)/deltaVal)*deltaVal;
      double phiUpper = phiLower + deltaVal;
      phiLower = phiLower + 360.0;
    }
    else if (targetPhiVal >= 360-deltaVal){

     double phiLower = deltaVal/2.0 + floor( (targetPhiVal-deltaVal/2.0)/deltaVal)*deltaVal;
     double phiUpper = phiLower + deltaVal - 360;
    }
    else{
    double phiLower = deltaVal/2.0 + floor( (targetPhiVal-deltaVal/2.0)/deltaVal)*deltaVal;
    double phiUpper = phiLower + deltaVal;
    }

   */






  //  qDebug() << QString::number(targetPhiVal) << " mapped to: "<< QString::number(phiLowerIndex) << " " << QString::number(phiUpperIndex) << "\n";
 //   qDebug() << QString::number(targetPhiVal) << " mapped to left-edges: "<< QString::number(phiLowerIndex*5) << " " << QString::number(phiUpperIndex*5) << "\n";
  //  qDebug() << QString::number(targetPhiVal) << " mapped to bin-centers: "<< QString::number(phiLowerIndex*5 + 2.5) << " " << QString::number(phiUpperIndex*5 + 2.5) << "\n";






    //double foundEnergy = energyData[phiLowerIndex][thetaLowerIndex];
    //for now skip interpolation and just set the box to the closest match


    int closestI = (int)floor( (targetPhiVal)/deltaVal ) ;
    int closestJ =(int)floor( (targetThetaVal)/deltaVal ) ;
    if(closestI > 71){
        closestI = 0;
    }
    if(closestJ > 35){
        closestJ = 35;
    }


    double foundEnergy = energyData[closestI][closestJ];
    this->findChild<QLineEdit *>("energyOutBox")->setText(QString::number(foundEnergy)) ;




}


void MainWindow::on_findUAButton_clicked()
{
    QString uaDir = QFileDialog::getExistingDirectory(this, tr("UA Install Directory"), this->uaGlobalPath, QFileDialog::ShowDirsOnly);
    if(uaDir != ""){
        this->findChild<QLineEdit *>("uaFolderBox")->setText(uaDir);
        this->uaGlobalPath = uaDir;
        checkForMaterials();
        this->findChild<QLineEdit *>("resultFolderBox")->setText(   QDir::cleanPath(  this->uaGlobalPath   +"/uaquickrunresults"));

        if(!QFile::exists(uaDir+"/UnitedAtom")){
            //alert - no UA found

            QMessageBox msgBox;
            msgBox.setText("Warning: UnitedAtom executable not found in this location.");
               msgBox.exec();
               this->statusBar()->showMessage("UA not yet found");
        }
        else{
            this->statusBar()->showMessage("Ready");
        }
    }
}


void MainWindow::on_resultFolderButton_clicked()
{
    QString resultDir = QFileDialog::getExistingDirectory(this, tr("Results Directory"), this->uaGlobalPath, QFileDialog::ShowDirsOnly);
    this->findChild<QLineEdit *>("resultFolderBox")->setText(resultDir);
}


void MainWindow::on_loadMaterialButton_clicked()
{

    QString targetMaterialFile = QFileDialog::getOpenFileName(this, tr("Material File"),  this->uaGlobalPath,  tr("CSV (*.csv)"));
    //Load in the file and process
    loadMaterials(targetMaterialFile);
}

void MainWindow::loadMaterials( QString materialFile){

    QFile file(materialFile);
    if(file.open(QIODevice::ReadOnly | QIODevice::Text)) {
    QTextStream fileIn(&file);
    int numFound = 0;

    while(!fileIn.atEnd()){
        std::string lineIn = fileIn.readLine().toStdString();


        if(lineIn.substr(0, 1) != "#") {
            std::vector<std::string> results;
             boost::split(results, lineIn, [](char c){return c == ',';});
             //silicaquartz,surface/SiO2-Quartz,hamaker/SiO2_Quartz.dat,1
             //qDebug() << QString::fromStdString(results[0]) << " " << QString::fromStdString(results[1]) << " " <<  QString::fromStdString(results[2] )<< "\n";
             if(results.size()==4){
             materialTypes.emplace_back( MaterialType( QString::fromStdString(results[0])   , QString::fromStdString(results[1]) ,QString::fromStdString(results[2])     ) ) ;
             //addBeadTypeMaterialBox->addItem( QString::fromStdString(results[0]) ,   QList<QVariant>() <<  QString::fromStdString(results[1]) <<  QString::fromStdString(results[2])) ;
             this->findChild<QComboBox *>("materialDropdown")->addItem( QString::fromStdString(results[0]) ,   QList<QVariant>() <<  QString::fromStdString(results[1]) <<  QString::fromStdString(results[2]))  ;

             numFound += 1;
             }

        }
    }

    QMessageBox msgBox;
    msgBox.setText("Loaded " + QString::number(numFound) + " materials");
       msgBox.exec();

}

}


void MainWindow::on_pdbTargetButton_clicked()
{
    QString targetPDB = QFileDialog::getOpenFileName(this, tr("PDB Target"),  this->uaGlobalPath,  tr("PDB (*.pdb)"));
    this->findChild<QLineEdit *>("pdbTargetLine")->setText(targetPDB);
}


void MainWindow::updateUABox(){
    QString newText = QString (this->processHandle.readAllStandardOutput()) ;

    this->findChild<QPlainTextEdit *>("uaOutputBox")->appendPlainText(newText+"\n");
}


void MainWindow::uaDoneAlert(int exitCode, QProcess::ExitStatus exitStatus){
  qDebug() << "UA run complete \n";
  QMessageBox::warning(this, tr("UA Quick Run"),       tr("UA Run complete.\n" ) );
  this->findChild<QPushButton *>("runUAButton")->setDisabled(false);
  this->findChild<QLineEdit *>("loadPDBBox")->setText(this->findChild<QLineEdit *>("pdbTargetLine")->text());

  int targetRadius = this->findChild<QSpinBox *>("radiusSpinBox")->value() ;
  int targetZetaMV = this->findChild<QSpinBox *>("zetaSpinBox")->value() ;
   this->findChild<QSpinBox *>("npViewRadius")->setValue(targetRadius);
  QString targetOutputFile = (this->findChild<QLineEdit *>("resultFolderBox")->text())+"/NP1R_"+QString::number(targetRadius)+"_ZP_"+QString::number(targetZetaMV)+"/" ;
  this->findChild<QLineEdit *>("loadUAMBox")->setText(   QDir::cleanPath(targetOutputFile  ));
this->statusBar()->showMessage("UA run complete");


//this->processHandle.close();
}

void MainWindow::on_runUAButton_clicked()
{

    QString targetPDB = this->findChild<QLineEdit *>("pdbTargetLine")->text();
    QString targetMaterial = this->findChild<QComboBox *>("materialDropdown")->currentText() ;
    int targetRadius = this->findChild<QSpinBox *>("radiusSpinBox")->value() ;
    int targetZetaMV = this->findChild<QSpinBox *>("zetaSpinBox")->value() ;
    QString outputFolder = this->findChild<QLineEdit *>("resultFolderBox")->text();
    double targetZeta = targetZetaMV/1000.0;
    bool canRun = true;
    if( targetPDB==""){
        qDebug() << "no pdb specified \n";
        canRun = false;
    }
    if(targetMaterial == ""){
        qDebug() << "no material specified \n";
        canRun = false;
    }

    if(canRun == true){
        this->statusBar()->showMessage("Beginning UA run");
    //qDebug() << "Attempting to run: " << targetPDB << " on NP of material " << targetMaterial << " radius: " << QString::number(targetRadius) <<   "zeta: " << QString::number(targetZetaMV) <<  "\n";
    QString argString = "--operation-type=pdb -p "+targetPDB+" -r "+QString::number(targetRadius)+" -z "+QString::number(targetZeta)+" -m "+targetMaterial+" -o "+outputFolder;
    //qDebug() << argString << "\n";
    QStringList commandArgs;
    commandArgs << "RunUA.py";
    //commandArgs << argString;
    commandArgs << "--operation-type=pdb";
    commandArgs << "-p";
    commandArgs << targetPDB;
    commandArgs <<"-r";
    commandArgs<< QString::number(targetRadius) ;
    commandArgs <<"-z";
    commandArgs <<QString::number(targetZeta) ;
    commandArgs <<"-m";
    commandArgs <<targetMaterial ;
     commandArgs <<"-o";
    commandArgs << outputFolder;
    commandArgs << "-P";
    commandArgs << "0";

    if(this->findChild<QCheckBox *>("autoNPBox")->isChecked() == false){
        QString npTarget = this->findChild<QLineEdit *>("npTargetBox")->text();
       commandArgs << "-N";
       commandArgs << npTarget;
    }


    this->findChild<QPlainTextEdit *>("uaOutputBox")->appendPlainText(commandArgs.join(" ")+"\n");
    this->processHandle.setWorkingDirectory(uaGlobalPath );
    processHandle.setProcessChannelMode(QProcess::MergedChannels) ;

    // connect(&processHandle, SIGNAL ( readyReadStandardOutput() )  , this,  SLOT( updateUABox() )        );
   //  connect(&processHandle, SIGNAL(finished(int, QProcess::ExitStatus) ), this, SLOT(uaDoneAlert(int, QProcess::ExitStatus)));
     this->findChild<QPushButton *>("runUAButton")->setDisabled(true);
    this->processHandle.start("python3", commandArgs) ;

    }
}


void MainWindow::on_phiInputBox_valueChanged(int arg1)
{
    updateEnergyBox() ;
    updateMoleculeBox();
}


void MainWindow::on_thetaInputBox_valueChanged(int arg1)
{
    updateEnergyBox() ;
    updateMoleculeBox();
}

void MainWindow::getSceneMouseClick(QPointF scenePosLoc ){
    QGraphicsView* heatmapBox = this->findChild<QGraphicsView *>("heatmapView") ;
    //int gbWidth = heatmapBox->width() ;
    //int gbHeight = heatmapBox->height() ;
    double scaleFactorH =   360.0/(uamBoxWidth*0.95);
    double scaleFactorV =    180.0/(uamBoxHeight*0.95);
    //qDebug() << "got click at " << QString::number(scenePosLoc.x()) << "," << QString::number(scenePosLoc.y()) <<  "\n";
    double equivPhi = scenePosLoc.x() * scaleFactorH  ;
    //double equivTheta = (scenePosLoc.y() ) * scaleFactorV;
    double equivTheta = 180.0 -  (scenePosLoc.y() ) * scaleFactorV;
     //qDebug() << "mapping to " << QString::number(equivPhi) << "," << QString::number(equivTheta) <<  "\n";

     equivPhi = std::min( 360.0, equivPhi);
     equivPhi = std::max(0.0, equivPhi);
     equivTheta = std::min(180.0, equivTheta);
     equivTheta = std::max(0.0, equivTheta);

     this->findChild<QSpinBox *>("phiInputBox")->setValue( (int)round(equivPhi) );
     this->findChild<QSpinBox *>("thetaInputBox")->setValue( (int)round(equivTheta) );
}

void MainWindow::on_loadPDBButton_clicked()
{
    QString pdbBasePath = this->uaGlobalPath;
    QString presetPDBPath = this->findChild<QLineEdit *>("pdbTargetLine")->text();
    if(presetPDBPath != ""){
    pdbBasePath = presetPDBPath;
    }

    QString currentTarget = this->findChild<QLineEdit *>("loadPDBBox")->text();
    if(currentTarget != ""){
        pdbBasePath = currentTarget;
    }


    QString targetPDBFile = QFileDialog::getOpenFileName(this, tr("PDB File"),  pdbBasePath , tr("PDB-file (*.pdb)"));
    this->findChild<QLineEdit *>("loadPDBBox")->setText(targetPDBFile);
    //load in the .uam file
    if(targetPDBFile!=""){
        this->statusBar()->showMessage("Loading PDB");
        QFile file(targetPDBFile);
        if(file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QTextStream fileIn(&file);

        atomList.clear();
        //load in the data
        double xcenter = 0;
        double ycenter = 0;
        double zcenter = 0;
        int numAtoms = 0;

        while(!fileIn.atEnd()){
            std::string lineIn = fileIn.readLine().toStdString();
            if(lineIn.substr(0,4) =="ATOM" && lineIn.substr(13, 2)=="CA") {

                //load in PDB lines
                std::string nameIn =lineIn.substr(17, 3);
                double x = (0.1 * std::stod(lineIn.substr(30, 8)));
                double y=(0.1 * std::stod(lineIn.substr(38, 8)));
                double z=(0.1 * std::stod(lineIn.substr(46, 8)));
                  atomList.emplace_back( Atom(nameIn,x,y,z)) ;
                  xcenter += x;
                  ycenter += y;
                  zcenter += z;
                  numAtoms +=1;
            }
        }
       xcenter = xcenter/(numAtoms+0.0001);
        ycenter = ycenter/(numAtoms+0.0001);
         zcenter = zcenter/(numAtoms+0.0001);
     //    qDebug() << QString::number( atomList[0].x0 ) << "\n";
        for( auto& atom: atomList){
          //  qDebug() << " x0 was: " << QString::number( atom.x0 ) << "\n";
            atom.x0 = atom.x0 - xcenter;
            atom.y0 = atom.y0 - ycenter;
            atom.z0 = atom.z0 - zcenter;

            atom.orientationTheta = acos( atom.z0 / sqrt(atom.x0*atom.x0 + atom.y0*atom.y0 + atom.z0*atom.z0) ) ;
            if(atom.y0 >= 0){
            atom.orientationPhi =   acos( atom.x0/(sqrt(atom.x0*atom.x0 + atom.y0*atom.y0))) ;
            }
            else{
             atom.orientationPhi =   2*3.1415 -   acos( atom.x0/(sqrt(atom.x0*atom.x0 + atom.y0*atom.y0))) ;
            }

         //   qDebug() << " x0 now: " << QString::number( atom.x0 ) << "\n";
        }
     //   qDebug() << QString::number( atomList[0].x0 ) << "\n";

}
}
    this->statusBar()->showMessage("Finished loading PDB");
updateMoleculeBox();
}

void MainWindow::calcBeadDistances(bool doRotate = false){
    //for each bead, calculate its boltzmann-weighted distance to the surface of the NP

    double npRadius  = this->findChild<QSpinBox *>("npViewRadius")->value();
    double beadDelta = 0.5 + 0.2 + npRadius ;

    this->statusBar()->showMessage("Beginning distance calculation");


     //first pass: initialise the working variables
       for( auto& atom: atomList){
           atom.xc2 = atom.x0;
           atom.yc2 = atom.y0;
           atom.zc2 = atom.z0;
           atom.dAv = 0;
           atom.dAvn = 0;
           atom.dAvd = 0.000001;
       }
//for loop i for loop j

    int i = 2;
    int j = 3;

    for( i = 0; i<72;++i){
        for(j = 0; j<36;++j){



    double currentPhi = (  i*5 + 2.5)* 3.1415/180.0;
    double currentTheta = (j*5 + 2.5)* 3.1415/180.0;
    double energyVal = energyData[i][j];
    double p = -1.0*currentPhi;
    double t = 3.1415 - currentTheta;
    double omega = this->findChild<QDial *>("omegaDial")->value() * 3.1415/180.0;
    double extent = 1;
    double minz = 0;

     double rotateMatrix[3][3];
     rotateMatrix[0][0] = std::cos(t) * std::cos(p) * std::cos(omega) - std::sin(omega) * std::sin(p);
     rotateMatrix[0][1] = -1.0 * std::cos(t) * std::sin(p) * std::cos(omega) - std::cos(p)*std::sin(omega);
     rotateMatrix[0][2] = std::sin(t)*std::cos(omega);
     rotateMatrix[1][0] = std::sin(p)*std::cos(omega) + std::cos(p) * std::cos(t)*std::sin(omega);
     rotateMatrix[1][1] = std::cos(p)*std::cos(omega) - std::cos(t)*std::sin(omega)*std::sin(p);
     rotateMatrix[1][2] = std::sin(omega) * std::sin(t);
     rotateMatrix[2][0] = -1.0 * std::sin(t) * std::cos(p);
     rotateMatrix[2][1] = std::sin(t) * std::sin(p);
     rotateMatrix[2][2] = std::cos(t);


double xAbsMax = 0;
double yAbsMax = 0;
double zAbsMax = 0;


//second pass: rotate the biomolecule and store its coordinates


    for( auto& atom: atomList){

        double xc = atom.x0*rotateMatrix[0][0] + atom.y0 * rotateMatrix[0][1] + atom.z0 * rotateMatrix[0][2];
        double yc = atom.x0*rotateMatrix[1][0] + atom.y0 *rotateMatrix[1][1] + atom.z0 * rotateMatrix[1][2];
        double zc = atom.x0*rotateMatrix[2][0] + atom.y0 * rotateMatrix[2][1] + atom.z0 * rotateMatrix[2][2];
        atom.xc2 = xc;
        atom.yc2 = yc;
        atom.zc2 = zc;
        minz = std::min( minz, zc);
        xAbsMax = std::max( xAbsMax, abs(xc));
        yAbsMax = std::max( yAbsMax, abs(yc));
        zAbsMax = std::max( zAbsMax, abs(zc));
        extent = std::max( extent, sqrt(xc*xc + yc*yc + zc*zc));
    }
   // double proteinRadius = extent;


//third pass: get the distance to the NP for this orientation

    for( auto& atom: atomList){
        double zOffset = atom.zc2 -minz + beadDelta;
     double beadDist = sqrt( atom.xc2*atom.xc2   + atom.yc2*atom.yc2 + zOffset*zOffset ) - npRadius;
     atom.dAvn += beadDist * exp(-energyVal)*std::sin(currentTheta);
     atom.dAvd += exp(-energyVal)*std::sin(currentTheta);

    }

        }
    }

    //final pass: set the average distance

    for( auto& atom: atomList){
     atom.dAv = atom.dAvn / (atom.dAvd + 0.01);
     atom.colourParam = atom.dAv;
    }
 this->statusBar()->showMessage("Distances calculated");
    updateMoleculeBox();
}

void MainWindow::updateMoleculeBox(){
 //"pdbScene" for the painting scene
    // "pdbView" for the QGraphicsView box

    QGraphicsView* graphicsBox = this->findChild<QGraphicsView *>("pdbView") ;
    int gbWidth = graphicsBox->width() ;
    int gbHeight = graphicsBox->height() ;
    int scaleBasis = gbWidth;
    if(gbHeight < gbWidth){
        scaleBasis = gbHeight;
    }

    std::vector< std::vector<double> > plotObjects;


    double currentPhi = this->findChild<QSpinBox *>("phiInputBox")->value() * 3.1415/180.0;
    double currentTheta = this->findChild<QSpinBox *>("thetaInputBox")->value() * 3.1415/180.0;

    double p = -1.0*currentPhi;
    double t = 3.1415 - currentTheta;
    double omega = this->findChild<QDial *>("omegaDial")->value() * 3.1415/180.0;
    double extent = 1;
    double minz = 0;

     double rotateMatrix[3][3];
     rotateMatrix[0][0] = std::cos(t) * std::cos(p) * std::cos(omega) - std::sin(omega) * std::sin(p);
     rotateMatrix[0][1] = -1.0 * std::cos(t) * std::sin(p) * std::cos(omega) - std::cos(p)*std::sin(omega);
     rotateMatrix[0][2] = std::sin(t)*std::cos(omega);
     rotateMatrix[1][0] = std::sin(p)*std::cos(omega) + std::cos(p) * std::cos(t)*std::sin(omega);
     rotateMatrix[1][1] = std::cos(p)*std::cos(omega) - std::cos(t)*std::sin(omega)*std::sin(p);
     rotateMatrix[1][2] = std::sin(omega) * std::sin(t);
     rotateMatrix[2][0] = -1.0 * std::sin(t) * std::cos(p);
     rotateMatrix[2][1] = std::sin(t) * std::sin(p);
     rotateMatrix[2][2] = std::cos(t);
pdbScene.clear( );

double xAbsMax = 0;
double yAbsMax = 0;
double zAbsMax = 0;

double npRadius  = this->findChild<QSpinBox *>("npViewRadius")->value();
//first pass: rotate the biomolecule and store its coordinates


    for( auto& atom: atomList){

        double xc = atom.x0*rotateMatrix[0][0] + atom.y0 * rotateMatrix[0][1] + atom.z0 * rotateMatrix[0][2];
        double yc = atom.x0*rotateMatrix[1][0] + atom.y0 *rotateMatrix[1][1] + atom.z0 * rotateMatrix[1][2];
        double zc = atom.x0*rotateMatrix[2][0] + atom.y0 * rotateMatrix[2][1] + atom.z0 * rotateMatrix[2][2];
        atom.xc = xc;
        atom.yc = yc;
        atom.zc = zc;
        minz = std::min( minz, zc);
        xAbsMax = std::max( xAbsMax, abs(xc));
        yAbsMax = std::max( yAbsMax, abs(yc));
        zAbsMax = std::max( zAbsMax, abs(zc));
        extent = std::max( extent, sqrt(xc*xc + yc*yc + zc*zc));



    }
    double proteinRadius = extent;

   double beadDelta = 0.5 + 0.2 + npRadius ; //offset to apply to beads to displace their center from the NP center at (0,0) by one bead radius + a separation + the NP radius;
double zMax = 0;

double colourParamMin = 0;
double colourParamMax = 0.01;

   //final pass: get the set of circles to plot, translating such that in the world-frame the molecule is above the NP and has its minimum point defined relative to the NP center at (0,0,0)
    for( auto& atom: atomList){
        double colourParam = atom.colourParam ;
        colourParamMin = std::min(colourParam,colourParamMin);
        colourParamMax = std::max(colourParam,colourParamMax);
        std::vector<double> plotCircle{ atom.xc,atom.yc,(atom.zc - minz) + beadDelta ,0.5, colourParam} ;
        zMax = std::max( zMax, (atom.zc - minz) + beadDelta ); //get the maximum z-coordinate used to redefine the uppermost point for transformation to graphics co-ords
        plotObjects.emplace_back(  plotCircle );


    }

 //sort by y-coord to get an approximation of depth
    std::sort(plotObjects.begin(), plotObjects.end(),
              [](const std::vector<double>& a, const std::vector<double>& b) {
      return a[1] < b[1];
    });


QPen outlinePen(Qt::black);
QBrush whiteFill(Qt::white);
QBrush greyFill(Qt::gray);

extent = std::max(1.0, zMax);
extent = std::max(extent, xAbsMax);
extent = std::max(extent, yAbsMax);

double scaleFactor = 0.5*scaleBasis/(extent + 0.01);
 double boundSize = scaleFactor*( 1 + 1.5*proteinRadius);
outlinePen.setWidth( boundSize/100);

for(const auto& pc: plotObjects){
 //qDebug() << " plotting bead at max-y-value " << pc[2] << "\n";
 double xleft = pc[0] - pc[3];
 double zup =  zMax -pc[2]  - pc[3];
 double colourVal = ( colourParamMax - pc[4]   )/(0.01 +  colourParamMax - colourParamMin ) ;
 QBrush colourFill(  QColor(255*colourVal,0,0)  );

 pdbScene.addEllipse(xleft*scaleFactor, zup*scaleFactor, 2*pc[3]*scaleFactor , 2*pc[3]*scaleFactor  ,outlinePen , colourFill);
}

double npLeft = -   npRadius;
double npUp = zMax -( 0 ) - npRadius;
 pdbScene.addEllipse( npLeft*scaleFactor, npUp*scaleFactor, 2*npRadius*scaleFactor , 2*npRadius*scaleFactor  ,outlinePen ,greyFill);
pdbScene.setSceneRect( (-proteinRadius)*scaleFactor, 0 , (2*proteinRadius)*scaleFactor, (2*proteinRadius)*scaleFactor);

 graphicsBox->setScene( &pdbScene );
 graphicsBox->show();


 graphicsBox->fitInView(  -boundSize, -boundSize, 2*boundSize, 2*boundSize  );



 //Add the scale bar

 QGraphicsView* pdbScaleBox = this->findChild<QGraphicsView *>("pdbColourScaleBar") ;
 int scaleBarWidth = pdbScaleBox->width() ;
 int scaleBarHeight = pdbScaleBox->height() ;



pdbScalePixmap.fill();
 QPainter *scalePaint = new QPainter(&pdbScalePixmap);
 QFont font = scalePaint->font();
 //font.setPixelSize(10 );
 scalePaint->setFont(font);
 int barWidth = 400/10;
 for(int i =0; i<10;i++){
     double distVal = colourParamMin + i*(colourParamMax - colourParamMin)/10.0;

     double colourVal = ( colourParamMax - distVal   )/(0.01 +  colourParamMax - colourParamMin ) ;
     QBrush colourFill(  QColor(255*colourVal,0,0)  );

     //qDebug() << QString::number(eScale) << " from " << QString::number(energyVal) << "\n";
     double y1 = 0;
     double x1 = i*barWidth;
     scalePaint->setPen(QColor(255*colourVal,0,0));
      scalePaint->setBrush(QColor(255*colourVal,0,0));
     scalePaint->drawRect(x1, y1, barWidth, 20  );
      scalePaint->setPen(QColor(0,0,0));
      //scalePaint->drawText( x1, y1+20, 10, 20, Qt::AlignLeft,  QString::number(energyVal, 'f',1)) ;
      if(i%2 == 0){

          if(distVal > 1){
               scalePaint->drawText( x1, y1+40,   QString::number((int)round(distVal))) ;
          }
          else{
      scalePaint->drawText( x1, y1+40,   QString::number(distVal, 'f',1)) ;
          }
      }
 }
 delete scalePaint;



 pdbBarScene.addPixmap(pdbScalePixmap.scaled(scaleBarWidth*0.95,scaleBarHeight*0.95)   );
 pdbScaleBox->setScene(&pdbBarScene);


}

void MainWindow::on_radiusSpinBox_valueChanged(int arg1)
{
    this->findChild<QSpinBox *>("npViewRadius")->setValue(arg1);
}


void MainWindow::on_npViewRadius_valueChanged(int arg1)
{
    updateMoleculeBox();
}


void MainWindow::on_omegaDial_valueChanged(int value)
{
    updateMoleculeBox();
}

void MainWindow::checkForMaterials(){
   QString baseMaterialSet = QDir::cleanPath(  uaGlobalPath+"/MaterialSet.csv" );
  if(QFile::exists(baseMaterialSet ) ){
    loadMaterials( baseMaterialSet);
  }

}

void MainWindow::on_colourBoltz_clicked()
{
    calcBeadDistances();
        for( auto& atom: atomList){
            atom.colourParam = atom.dAv;
        }


}


void MainWindow::on_findMinEnergyButton_clicked()
{
    double minEnergyVal = 50;
    double minPhi = 0;
    double minTheta = 0;
    for(int i = 0; i<72;++i){
        for(int j = 0; j<36;++j){



    double currentPhi = (  i*5 + 2.5) ;
    double currentTheta = (j*5 + 2.5) ;
    if(energyData[i][j] < minEnergyVal){
        minEnergyVal = energyData[i][j];
        minPhi  = currentPhi;
        minTheta = currentTheta;

        }
    }


}
    this->findChild<QSpinBox *>("phiInputBox")->setValue(minPhi);
      this->findChild<QSpinBox *>("thetaInputBox")->setValue(minTheta);

}

void MainWindow::on_findBoltzMinButton_clicked()
{

    calcBeadDistances(); //precalculate the optimum displacements and record these
    double npRadius  = this->findChild<QSpinBox *>("npViewRadius")->value();
    double beadDelta = 0.5 + 0.2 + npRadius ;
    double minTotalDisplacement = 500;
    double bestTheta = 0;
    double bestPhi = 0;

    this->statusBar()->showMessage("Beginning distance calculation");
     //first pass: initialise the working variables without overwriting the optimum
       for( auto& atom: atomList){
           atom.xc2 = atom.x0;
           atom.yc2 = atom.y0;
           atom.zc2 = atom.z0;
       }
//for loop i for loop j

    int i = 2;
    int j = 3;

    for( i = 0; i<72;++i){
        for(j = 0; j<36;++j){



    double currentPhi = (  i*5 + 2.5)* 3.1415/180.0;
    double currentTheta = (j*5 + 2.5)* 3.1415/180.0;
    double energyVal = energyData[i][j];
    double p = -1.0*currentPhi;
    double t = 3.1415 - currentTheta;
    double omega = this->findChild<QDial *>("omegaDial")->value() * 3.1415/180.0;
    double extent = 1;
    double minz = 0;

     double rotateMatrix[3][3];
     rotateMatrix[0][0] = std::cos(t) * std::cos(p) * std::cos(omega) - std::sin(omega) * std::sin(p);
     rotateMatrix[0][1] = -1.0 * std::cos(t) * std::sin(p) * std::cos(omega) - std::cos(p)*std::sin(omega);
     rotateMatrix[0][2] = std::sin(t)*std::cos(omega);
     rotateMatrix[1][0] = std::sin(p)*std::cos(omega) + std::cos(p) * std::cos(t)*std::sin(omega);
     rotateMatrix[1][1] = std::cos(p)*std::cos(omega) - std::cos(t)*std::sin(omega)*std::sin(p);
     rotateMatrix[1][2] = std::sin(omega) * std::sin(t);
     rotateMatrix[2][0] = -1.0 * std::sin(t) * std::cos(p);
     rotateMatrix[2][1] = std::sin(t) * std::sin(p);
     rotateMatrix[2][2] = std::cos(t);


double xAbsMax = 0;
double yAbsMax = 0;
double zAbsMax = 0;


//second pass: rotate the biomolecule and store its coordinates


    for( auto& atom: atomList){

        double xc = atom.x0*rotateMatrix[0][0] + atom.y0 * rotateMatrix[0][1] + atom.z0 * rotateMatrix[0][2];
        double yc = atom.x0*rotateMatrix[1][0] + atom.y0 *rotateMatrix[1][1] + atom.z0 * rotateMatrix[1][2];
        double zc = atom.x0*rotateMatrix[2][0] + atom.y0 * rotateMatrix[2][1] + atom.z0 * rotateMatrix[2][2];
        atom.xc2 = xc;
        atom.yc2 = yc;
        atom.zc2 = zc;
        minz = std::min( minz, zc);
        xAbsMax = std::max( xAbsMax, abs(xc));
        yAbsMax = std::max( yAbsMax, abs(yc));
        zAbsMax = std::max( zAbsMax, abs(zc));
        extent = std::max( extent, sqrt(xc*xc + yc*yc + zc*zc));
    }
   // double proteinRadius = extent;


//third pass: get the distance to the NP for this orientation
    double orientationTotalDisplacement = 0;
    for( auto& atom: atomList){
        double zOffset = atom.zc2 -minz + beadDelta;
     double beadDist = sqrt( atom.xc2*atom.xc2   + atom.yc2*atom.yc2 + zOffset*zOffset ) - npRadius;

     orientationTotalDisplacement += pow(beadDist - atom.dAv ,2);


    }


    if(orientationTotalDisplacement < minTotalDisplacement){
        minTotalDisplacement = orientationTotalDisplacement;
        bestPhi = (  i*5 + 2.5);
        bestTheta = (  j*5 + 2.5);
    }

    if(i == 0 && j ==0){
        minTotalDisplacement = orientationTotalDisplacement;
        bestPhi = 2.5;
        bestTheta = 2.5;
    }

        }
    }



 this->statusBar()->showMessage("Distances calculated and optimum found");
    this->findChild<QSpinBox *>("phiInputBox")->setValue(bestPhi);
      this->findChild<QSpinBox *>("thetaInputBox")->setValue(bestTheta);
    updateMoleculeBox();

}


void MainWindow::on_colourEnergy_clicked()
{

    for( auto& atom: atomList){
        int i =(int)( (atom.orientationPhi*180.0/3.1415 - 2.5)/5.0);
        int j = (int)( (atom.orientationTheta*180.0/3.1415- 2.5)/5.0);
        qDebug() << QString::number(atom.orientationPhi*180.0/3.1415) << " mapped to " << QString::number(i) << "\n";
        atom.colourParam = energyData[i][j];
    }
updateMoleculeBox();
}


void MainWindow::on_autoNPBox_stateChanged(int arg1)
{
   // qDebug() << "Checkbox new state: " << QString::number(arg1) << "\n";
    if (this->findChild<QCheckBox *>("autoNPBox")->isChecked() == true){
      //box is currently checked, disable the find NP window
        this->findChild<QLineEdit *>("npTargetBox")->setDisabled(true);
        this->findChild<QPushButton *>("npTargetButton")->setDisabled(true);
    }
    else{
this->findChild<QLineEdit *>("npTargetBox")->setDisabled(false);
         this->findChild<QPushButton *>("npTargetButton")->setDisabled(false);
    }

}


void MainWindow::on_npTargetButton_clicked()
{
    QString targetNPFile = QFileDialog::getOpenFileName(this, tr("NP File"),  this->uaGlobalPath,  tr("NP (*.np)"));
    this->findChild<QLineEdit *>("npTargetBox")->setText(targetNPFile);
}

