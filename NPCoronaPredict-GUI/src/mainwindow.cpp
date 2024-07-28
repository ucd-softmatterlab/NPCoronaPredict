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
#include <map>

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

    QAction *addMoleculeAction = new QAction("Add Molecule", this);
    QAction *removeMoleculeAction = new QAction("Remove Molecule", this);

    connect( addMoleculeAction, SIGNAL(triggered()) , this, SLOT( addMoleculeToMedium())   );
    connect( removeMoleculeAction, SIGNAL(triggered()) , this, SLOT( removeMoleculeFromMedium())   );
    tableMenu.addAction( addMoleculeAction );
    tableMenu.addAction(removeMoleculeAction);
    //this->findChild<QLabel *>("heatmapPlotLabel")->setPixmap(this->heatmapPixmap);
    //this->materialTypes.emplace_back( MaterialType( "Custom"   , "",""     ) ) ;
    //addBeadType->findChild<QComboBox *>("materialTypeBox")->addItem( "Custom" ,   QList<QVariant>() <<  QString::fromStdString("") <<  QString::fromStdString("")) ;
    //addShell->findChild<QComboBox *>("materialTypeBox")->addItem( "Custom" ,   QList<QVariant>() <<  QString::fromStdString("") <<  QString::fromStdString("")) ;
    //this->findChild<QComboBox *>("materialDropdown")->addItem( "Custom" ,   QList<QVariant>() <<  QString::fromStdString("") <<  QString::fromStdString("")) ;
    connect( &networkManager, SIGNAL(finished(QNetworkReply*)), this, SLOT(downloadReplyFinished(QNetworkReply*)) );
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
    double piVal = fPi;
    double simpleAverageNum = 0;
    double simpleAverageDenom = 0.00001;
    double boltzAverageNum = 0;
    double boltzAverageDenom = 0.000001;
   double kbtVal = 1;
   double minPhiVal = 1;
   double minThetaVal = 1;

   currentFileOmega = 0;

    QString targetUAMFileBasePath = this->findChild<QLineEdit *>("resultFolderBox")->text();
double npRadius  = this->findChild<QSpinBox *>("npViewRadius")->value();
    if(lastUAMPath != ""){
        targetUAMFileBasePath = lastUAMPath;
    }
    QString targetUAMFile = QFileDialog::getOpenFileName(this, tr("UAM File"), targetUAMFileBasePath,  tr("UA-Output (*.uam)"));
    this->findChild<QLineEdit *>("loadUAMBox")->setText(targetUAMFile);
    //load in the .uam file
    if(targetUAMFile!=""){
        lastUAMPath = targetUAMFile;
        QFile file(targetUAMFile);
        if(file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QTextStream fileIn(&file);
        QFileInfo fileInfo(targetUAMFile);
        QString assumedPDBName = fileInfo.fileName().split("_")[0];
        //qDebug() << fileInfo.fileName();
        //qDebug() << assumedPDBName;
        QString trialPDBPath = uaGlobalPath+"/all_proteins/"+assumedPDBName+".pdb";
        //attempt 1 - test directly in all_proteins
         //  qDebug() << "first trial: " << trialPDBPath;
        if( QFile::exists(trialPDBPath) ){

        this->findChild<QLineEdit *>("loadPDBBox")->setText(trialPDBPath);
        }
        else{
        //attempt 2 - check in the local storage for this project
            QString uamFileDir = QFileInfo( targetUAMFile ).absolutePath() ;

            trialPDBPath = uamFileDir+"/../../proteins/"+assumedPDBName+".pdb";
           // qDebug() << "second trial: " << trialPDBPath;
            if( QFile::exists(trialPDBPath) ){
            this->findChild<QLineEdit *>("loadPDBBox")->setText(trialPDBPath);
            }
        }
        this->statusBar()->showMessage("Loading UAM");
        double minEnergy = 0;
        double maxEnergy = 1;
        double foundRadius = -1;
        //load in the data
        while(!fileIn.atEnd()){
            std::string lineIn = fileIn.readLine().toStdString();
            //extra processing for metadata contained in the comments section
            if(lineIn.substr(0,1) == "#"){
              if(lineIn.substr(0,8) == "#RADIUS:"){
              foundRadius = std::stod( lineIn.substr(9) );
              //qDebug() << "extracted radius: " << foundRadius << "\n";

              if(foundRadius > 0){
                  this->findChild<QSpinBox *>("npViewRadius")->setValue( (int)round(foundRadius) );
              }

              }
            }

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

        double comVal = -1.0; //use this as an error value to tell the plotting routine it needs to estimate an offset
        if( results.size() > 11){
            comVal = std::stod( results[10]);
            currentFileOmega = std::stod( results[8] ) * fPi/180.0;
        }
        comData[phiIndex][thetaIndex] = comVal;


        if(energyVal < minEnergy){
            minEnergy = energyVal;
            minPhiVal = phiVal + 2.5;
            minThetaVal = thetaVal + 2.5;
        }
        //minEnergy = std::min(minEnergy,energyVal);
        maxEnergy = std::max(maxEnergy, energyVal);
            }
        }





      this->findChild<QSpinBox *>("phiInputBox")->setValue(minPhiVal);
        this->findChild<QSpinBox *>("thetaInputBox")->setValue(minThetaVal);

      //update the boxes
        this->findChild<QLineEdit *>("boltzOutBox")->setText( QString::number(boltzAverageNum/boltzAverageDenom)) ;
        this->findChild<QLineEdit *>("simpleOutBox")->setText( QString::number(simpleAverageNum/simpleAverageDenom)) ;
        updateEnergyBox();

        //plot the results
       updateHeatmapPlot();


        //this->findChild<QGraphicsView *>("heatmapView")->addPixmap(this->heatmapPixmap);
      this->statusBar()->showMessage("UAM load complete");

    }
}


}







void MainWindow::updateHeatmapPlot()
{
    QGraphicsView* heatmapBox = this->findChild<QGraphicsView *>("heatmapView") ;
    int gbWidth = heatmapBox->width() ;
    int gbHeight = heatmapBox->height() ;

    QGraphicsView* heatmapScaleBox = this->findChild<QGraphicsView *>("heatmapScaleBar") ;
    int scaleBarWidth = heatmapScaleBox->width() ;
    int scaleBarHeight = heatmapScaleBox->height() ;


    uamBoxHeight = gbHeight;
    uamBoxWidth = gbWidth;


    int currentPhi = this->findChild<QSpinBox *>("phiInputBox")->value();
     int currentTheta =  this->findChild<QSpinBox *>("thetaInputBox")->value();

        double minEnergy = 0;
        double maxEnergy = 1;
        //first loop: get min and max energies
       for( int phiIndex = 0; phiIndex < 72; ++phiIndex){
       for(int thetaIndex = 0; thetaIndex< 36; ++thetaIndex){
        double energyVal = energyData[phiIndex][thetaIndex] ;

        if(energyVal < minEnergy){
            minEnergy = energyVal;
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

        //then paint a marker at the current phi, theta
        QPen arrowPen(QColor(0,80,0),2);
        QPen arrowPen2(QColor(0,80,0),2);
        paint->setPen(arrowPen);
        paint->setBrush(QColor(0,80,0));

        double t1 = 180-currentTheta;
        double lineLength = 6;
        double arrowLength = 3;

        double t2 = t1 +lineLength;
        double t3 =  t1 + arrowLength;
        if(currentTheta <90){
        t2 = t1 - lineLength;
                t3 = t1 - arrowLength;
        }
        double p2 = currentPhi+lineLength;
        double p3 = currentPhi+arrowLength;
        if(currentPhi > 180){
            p2 = currentPhi-lineLength;
            p3 = currentPhi-arrowLength;
        }

        paint->drawLine( QLineF(currentPhi,  t1,p2,t2));

paint->setPen(arrowPen2);
        paint->drawLine( QLineF(currentPhi,  t1,p3,t1));
        paint->drawLine( QLineF(currentPhi,  t1,currentPhi,t3));
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

        //update graphics

 // scene.addPixmap(heatmapPixmap );
         scene.addPixmap(heatmapPixmap.scaled(gbWidth*0.95,gbHeight*0.95));
       heatmapBox->setScene(&scene);
       heatmapBarScene.addPixmap(heatmapScalePixmap.scaled(scaleBarWidth*0.95,scaleBarHeight*0.95)   );
       heatmapScaleBox->setScene(&heatmapBarScene);

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
    double foundCOM =  comData[closestI][closestJ];
    if(foundCOM > 0){
       this->findChild<QLineEdit *>("comDistOutBox")->setText(QString::number(foundCOM)) ;
    }
    else{
         this->findChild<QLineEdit *>("comDistOutBox")->setText("?") ;
    }

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
    loadMaterials(targetMaterialFile, true);
}

void MainWindow::loadMaterials( QString materialFile, bool doAlert ){

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

    if(doAlert == true){
    QMessageBox msgBox;
    msgBox.setText("Loaded " + QString::number(numFound) + " materials");
       msgBox.exec();
    }

}

}


void MainWindow::on_pdbTargetButton_clicked()
{
    QString targetPDB;
    bool doCorona = this->findChild<QCheckBox *>("npcpModeBox")->isChecked();
    if(doCorona == true){
        targetPDB = QFileDialog::getOpenFileName(this, tr("Biomolecule Medium file"),  this->uaGlobalPath,  tr("CSV (*.csv)"));
    }
    else{
    targetPDB = QFileDialog::getOpenFileName(this, tr("PDB Target"),  this->uaGlobalPath,  tr("PDB (*.pdb)"));

    }

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
this->findChild<QPlainTextEdit *>("uaOutputBox")->appendPlainText("---------------------\n");
this->findChild<QPlainTextEdit *>("uaOutputBox")->appendPlainText("---------------------\n");
//this->processHandle.close();
}

void MainWindow::on_runUAButton_clicked()
{

    QString targetPDB = this->findChild<QLineEdit *>("pdbTargetLine")->text();
    QString targetMaterial = this->findChild<QComboBox *>("materialDropdown")->currentText() ;
    int targetRadius = this->findChild<QSpinBox *>("radiusSpinBox")->value() ;
    int targetZetaMV = this->findChild<QSpinBox *>("zetaSpinBox")->value() ;
    QString outputFolder = this->findChild<QLineEdit *>("resultFolderBox")->text();


    bool doCorona = this->findChild<QCheckBox *>("npcpModeBox")->isChecked();

    bool localBoltzMode = this->findChild<QCheckBox *>("boltzModeCheckBox")->isChecked();

    double targetZeta = targetZetaMV/1000.0;

    double jitterMag = this->findChild<QDoubleSpinBox *>("jitterSpinBox")->value();

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
    //QString argString = "--operation-type=pdb -p "+targetPDB+" -r "+QString::number(targetRadius)+" -z "+QString::number(targetZeta)+" -m "+targetMaterial+" -o "+outputFolder;
    //qDebug() << argString << "\n";




    QStringList commandArgs;

    if(doCorona == false){
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
    commandArgs << "-j";
    commandArgs << QString::number(jitterMag);

    commandArgs << "-S";
    commandArgs << QString::number(1 + this->findChild<QComboBox *>("npTargetShapeOverride")->currentIndex());


    if(localBoltzMode == true){
        commandArgs << "-B";
        commandArgs << "1";
    }



    }
    else{
        //set up NPCoronaPredict

        int autorunSetting = this->findChild<QComboBox *>("npcpModeOptions")->currentIndex() ;

        commandArgs <<"NPCoronaPredict.py";
        //commandArgs << "--steady" ; 
        commandArgs <<"-r";
        commandArgs<< QString::number(targetRadius) ;
        commandArgs <<"-z";
        commandArgs <<QString::number(targetZeta) ;
        commandArgs <<"-m";
        commandArgs <<targetMaterial ;
        commandArgs << "-a";
        commandArgs << QString::number(autorunSetting);
        commandArgs << "-j";
        commandArgs << QString::number(jitterMag);
        commandArgs << "-S";
        commandArgs << QString::number(1 + this->findChild<QComboBox *>("npTargetShapeOverride")->currentIndex());
        if(localBoltzMode == true){
            commandArgs << "-B";
            commandArgs << "1";
        }

        
        commandArgs << "-o";
        commandArgs << targetPDB; //targetPDB was overloaded with otherproteins file if this mode is requested
        commandArgs << "--steady" ;

        //construct a project name
        QString autoProjectName = "";
QString npname = "";
QString npUAName= "";
        if(this->findChild<QCheckBox *>("autoNPBox")->isChecked() == true){

            npname =  targetMaterial + "_" + QString::number(targetRadius) + "_" + QString::number(targetZeta) ;
            autoProjectName = npname;
            npUAName = "np1R_" + QString::number(targetRadius) +"_ZP_" + QString::number(targetZeta)  ;
        }
        else{
            QString npTarget = this->findChild<QLineEdit *>("npTargetBox")->text();
           QStringList npTargetParts = npTarget.split("/");
            npname = npTargetParts.at( npTargetParts.size()-1);
             npname.chop(3);
            autoProjectName = npname ;
            npUAName = npname;
        }



        QStringList pdbTargetParts = targetPDB.split("/");
        QString mediumName = pdbTargetParts.at( pdbTargetParts.size()-1) ;
        mediumName.chop(4);
        autoProjectName = autoProjectName + "_" +  mediumName;


          this->findChild<QLineEdit *>("resultFolderBox")->setText(uaGlobalPath+"/CoronaPredictionProjects/"+autoProjectName+"/results/"+npUAName) ;

       commandArgs << "-p";
        commandArgs << autoProjectName;
    qDebug() << "Auto project name: " << autoProjectName ;
    }

    if(this->findChild<QCheckBox *>("autoNPBox")->isChecked() == false){
        QString npTarget = this->findChild<QLineEdit *>("npTargetBox")->text();
       commandArgs << "-N";
       commandArgs << npTarget;

    }
   
this->findChild<QPlainTextEdit *>("uaOutputBox")->appendPlainText("---------------------\n");
    this->findChild<QPlainTextEdit *>("uaOutputBox")->appendPlainText(commandArgs.join(" ")+"\n");
    this->findChild<QPlainTextEdit *>("uaOutputBox")->appendPlainText("---------------------\n");
    this->processHandle.setWorkingDirectory(uaGlobalPath );
    processHandle.setProcessChannelMode(QProcess::MergedChannels) ;

    // connect(&processHandle, SIGNAL ( readyReadStandardOutput() )  , this,  SLOT( updateUABox() )        );
   //  connect(&processHandle, SIGNAL(finished(int, QProcess::ExitStatus) ), this, SLOT(uaDoneAlert(int, QProcess::ExitStatus)));
     this->findChild<QPushButton *>("runUAButton")->setDisabled(true);
    this->findChild<QPushButton *>("cancelRunButton")->setDisabled(false);
    this->processHandle.start("python3", commandArgs) ;

    }
}


void MainWindow::on_phiInputBox_valueChanged(int arg1)
{
    updateEnergyBox() ;
    updateMoleculeBox();
    updateHeatmapPlot();
}


void MainWindow::on_thetaInputBox_valueChanged(int arg1)
{
    updateEnergyBox() ;
    updateMoleculeBox();
    updateHeatmapPlot();
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
                double radiusVal = 0.5;
                double chargeVal = 0.0;
                auto beadSearch = beadTypeMap.find(nameIn);
                if(  beadSearch!=beadTypeMap.end()   ){
                    radiusVal = beadSearch->second.radius;
                    chargeVal = beadSearch->second.charge;
                }

                  atomList.emplace_back( Atom(nameIn,x,y,z,false,radiusVal,chargeVal)) ;
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
             atom.orientationPhi =   2*fPi -   acos( atom.x0/(sqrt(atom.x0*atom.x0 + atom.y0*atom.y0))) ;
            }

         //   qDebug() << " x0 now: " << QString::number( atom.x0 ) << "\n";
        }
     //   qDebug() << QString::number( atomList[0].x0 ) << "\n";

      double shrinkwrapRadius = 0.2;
      double npRadius = 5;
        //then build the shrinkwrap layer for when this is shown
     for(int j=0; j<36;++j){
         for(int i=0; i<71;++i){
             double phiVal = (2.5 +i*5)* fPi/180.0;
             double thetaVal =(2.5 + j*5)* fPi/180.0;
             double radius = 0.5;

             double x0i = radius*cos(phiVal)*sin(thetaVal);
             double y0i = radius*sin(phiVal)*sin(thetaVal);
             double z0i = radius*cos(thetaVal);

             //loop over all existing non-shrinkwrap atoms, calculate the place where the bead centre would need to be to be outside
             for(auto& atom: atomList){
                 if(atom.isShrinkWrap==false){

                  double ri = sqrt( atom.x0*atom.x0 + atom.y0*atom.y0 + atom.z0 * atom.z0);
                  double dotProduct = x0i*atom.x0 + y0i*atom.y0 + z0i*atom.z0 ;
                  if(ri < 0.001){
                      ri = 0.001;
                      dotProduct = x0i*atom.x0 + y0i*atom.y0 + z0i*(atom.z0-0.01) ;
                  }
                  double inclAngle = acos(dotProduct/((ri)*(0.5) + 0.001 ));
                  double rb = atom.radius  +shrinkwrapRadius;
                  if(ri * sin(inclAngle) < rb  ){
                      double radius1 = ri * cos(inclAngle) + sqrt(rb*rb - (ri * sin(inclAngle))*(ri * sin(inclAngle)) ) ;
                      radius = std::max(radius,radius1 -  shrinkwrapRadius  );
                  }
                 }
             }

             double x = radius*cos(phiVal)*sin(thetaVal);
             double y = radius*sin(phiVal)*sin(thetaVal);
             double z = radius*cos(thetaVal);
             Atom newAtom = Atom("sw",x,y,z,true,shrinkwrapRadius);
             newAtom.orientationTheta =  thetaVal;
             newAtom.orientationPhi = phiVal ;
             atomList.emplace_back( newAtom) ;
         }
     }


}
}
    this->statusBar()->showMessage("Finished loading PDB");
updateMoleculeBox();
}

void MainWindow::calcBeadDistances(bool doRotate = false){
    //for each bead, calculate its distance to the surface of the NP in the current phi, theta configuration

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
       double targetPhiVal =( this->findChild<QSpinBox *>("phiInputBox")->value()   );
       double targetThetaVal = this->findChild<QSpinBox *>("thetaInputBox")->value();
       double deltaVal = 5.0;

    int closestI = (int)floor( (targetPhiVal)/deltaVal ) ;
    int closestJ =(int)floor( (targetThetaVal)/deltaVal ) ;
    if(closestI > 71){
        closestI = 0;
    }
    if(closestJ > 35){
        closestJ = 35;
    }
    int i = closestI;
    int j = closestJ;

    double currentPhi = (  i*5 + 2.5)* fPi/180.0;
    double currentTheta = (j*5 + 2.5)* fPi/180.0;



    double energyVal = energyData[i][j];
    double comVal = comData[i][j];


    double p = -1.0*currentPhi;
    double t = fPi- currentTheta;
    double omega = currentFileOmega;
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
        if(comVal > -1.0){
            zOffset = atom.zc2 + comVal;
        }
     int npShapeIndex = this->findChild<QComboBox *>("npShapeBox")->currentIndex();

     int xFactor = 1;
     if(npShapeIndex == 1){
         xFactor = 0;
     }

     double beadDist = sqrt( xFactor*atom.xc2*atom.xc2   + atom.yc2*atom.yc2 + zOffset*zOffset ) - npRadius;
     atom.dAvn += beadDist * exp(-energyVal)*std::sin(currentTheta);
     atom.dAvd += exp(-energyVal)*std::sin(currentTheta);

    }




    //final pass: set the average distance

    for( auto& atom: atomList){
     atom.dAv = atom.dAvn / (atom.dAvd + 0.01);
     atom.colourParam = atom.dAv;
    }
 this->statusBar()->showMessage("Distances calculated");
    updateMoleculeBox();
}


void MainWindow::calcBeadBoltzDistances(bool doRotate = false){
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



    double currentPhi = (  i*5 + 2.5)* fPi/180.0;
    double currentTheta = (j*5 + 2.5)* fPi/180.0;
    double energyVal = energyData[i][j];
    double comVal = comData[i][j];


    double p = -1.0*currentPhi;
    double t = fPi - currentTheta;
    double omega = this->findChild<QDial *>("omegaDial")->value() * fPi/180.0;
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
    int npShapeIndex = this->findChild<QComboBox *>("npShapeBox")->currentIndex();

    int xFactor = 1;
    if(npShapeIndex == 1){
        xFactor = 0;
    }
    for( auto& atom: atomList){


        double zOffset = atom.zc2 -minz + beadDelta;
        if(comVal > -1.0){
            zOffset = atom.zc2 + comVal;
        }

     double beadDist = sqrt(xFactor* atom.xc2*atom.xc2   + atom.yc2*atom.yc2 + zOffset*zOffset ) - npRadius;

     if(npShapeIndex == 2){
         beadDist = std::sqrt(  std::pow(std::max(0.0, std::abs(atom.xc2) - npRadius )  ,2) +  std::pow(std::max(0.0, std::abs(atom.yc2 ) - npRadius )  ,2) +   std::pow(std::max(0.0, std::abs(zOffset ) - npRadius )  ,2)   ) ;;
     }
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
double cylinderHalfLength = 10;
    QGraphicsView* graphicsBox = this->findChild<QGraphicsView *>("pdbView") ;
    int gbWidth = graphicsBox->width() ;
    int gbHeight = graphicsBox->height() ;
    int scaleBasis = gbWidth;
    if(gbHeight < gbWidth){
        scaleBasis = gbHeight;
    }

    std::vector< std::vector<double> > plotObjects;
    std::vector< std::vector<double> > shrinkwrapplotObjects;

    double currentPhi = this->findChild<QSpinBox *>("phiInputBox")->value() * fPi/180.0;
    double currentTheta = this->findChild<QSpinBox *>("thetaInputBox")->value() * fPi/180.0;

    double p = -1.0*currentPhi;
    double t = fPi - currentTheta;
    double screenRotate = this->findChild<QDial *>("omegaDial")->value() * fPi/180.0  ;
    double omega = screenRotate + currentFileOmega;
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
int maxNPCircles =7;

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
         cylinderHalfLength = std::max( cylinderHalfLength, sqrt(xc*xc + yc*yc)*1.1 );


    }
    double proteinRadius = extent;

   double beadDelta = 0.5 + 0.2 + npRadius ; //offset to apply to beads to displace their center from the NP center at (0,0) by one bead radius + a separation + the NP radius;

   int phiIndex =(int)( floor( (currentPhi*180.0/fPi )/5.0));
   int thetaIndex = (int)( floor((currentTheta*180.0/fPi)/5.0));
   phiIndex = std::max(0,phiIndex);
   thetaIndex = std::max(0,thetaIndex);
   phiIndex = std::min(71,phiIndex);
   thetaIndex = std::min(35,thetaIndex);


   double comVal = comData[phiIndex][thetaIndex];
   double plotComVal = comVal;
   //qDebug() << phiIndex << " " << thetaIndex << " found comval " << comVal << "\n";
   if(comVal > -1.0){
       beadDelta = comVal;
   }
    else{
       beadDelta = beadDelta - minz;
       plotComVal = beadDelta;
   }



   double zMax = 0;
    //qDebug() << " z offset applied: " << beadDelta << "\n";
double colourParamMin = 0;
double colourParamMax = 0.01;

bool make3D = false;

   //final pass: get the set of circles to plot, translating such that in the world-frame the molecule is above the NP and has its minimum point defined relative to the NP center at (0,0,0)
    for( auto& atom: atomList){

        if( atom.isShrinkWrap == false || showEnergyShrinkWrap == true){
        double colourParam = atom.colourParam ;
        colourParamMin = std::min(colourParam,colourParamMin);
        colourParamMax = std::max(colourParam,colourParamMax);
        double doOutLine = 1.0;
        if(atom.isShrinkWrap == true){
            doOutLine = 0.0;
        }


        std::vector<double> plotCircle{ atom.xc,atom.yc,(atom.zc ) + beadDelta ,atom.radius, colourParam, doOutLine, atom.charge} ;
       // qDebug() << QString::number(atom.xc) << " " << QString::number(atom.yc) << QString::number(atom.zc) ;
        zMax = std::max( zMax, (atom.zc ) + beadDelta ); //get the maximum z-coordinate used to redefine the uppermost point for transformation to graphics co-ords
        plotObjects.emplace_back(  plotCircle );
        int numLines = 8;
        if(make3D==true){
            for(int i =1; i<numLines; i++){
                  double yDelta = i * atom.radius/(numLines + 1);
                  double localRadius = std::sqrt( atom.radius*atom.radius - yDelta*yDelta);
                  std::vector<double> plotCircleN{ atom.xc,atom.yc + yDelta,(atom.zc ) + beadDelta , localRadius, colourParam, 0.0, atom.charge} ;



                  plotObjects.emplace_back( plotCircleN);
            }
        }

    }

    }

    //add in NP objects, remembering to get these in the correct position
    int numHalfCylinderSegments = 4;
    int npShapeIndex = this->findChild<QComboBox *>("npShapeBox")->currentIndex();
    if(npShapeIndex == 1){

        int circlesNeeded =  1 ;//+floor(   cylinderHalfLength/1.0 );
         double cylinderFaceWidth = std::abs( 2*npRadius * std::sin(screenRotate) ) ;
         double sideSegmentWidth = std::abs( 1.0 * std::cos(screenRotate) ) ;
         for(int i= 1; i<=circlesNeeded; i++){
             double xc1 = std::cos(screenRotate) * cylinderHalfLength*i/circlesNeeded;
             double yc1 = std::sin(screenRotate)*cylinderHalfLength*i/circlesNeeded;
             double xc3 = -xc1;
             double yc3 = -yc1;
             std::vector<double> plotCircle1{xc1    ,yc1,zMax- npRadius ,npRadius, -1, 1, cylinderFaceWidth/2.0, 0} ;
             std::vector<double> plotCircle3{xc3,yc3,zMax- npRadius ,npRadius, -1, 1,  cylinderFaceWidth/2.0 , 0} ;
                  plotObjects.emplace_back(  plotCircle1 );
                  plotObjects.emplace_back(  plotCircle3 );

         }

         //pc[3] stores the x offset for the lower points, pc[6] stores the z offset for lower points

         int numFaces = 16;

         double deltaPhi = 1.0/numFaces *  2 * fPi;
         double dl=cylinderHalfLength/(numHalfCylinderSegments);
         for(int k=-numHalfCylinderSegments+1; k<numHalfCylinderSegments+1; k++){
         for( int i =0; i< numFaces; i++){
             double phiAngle = (1.0*i)/numFaces * 2 *  fPi ;

             double zoffset = npRadius * std::sin(phiAngle) ;
             double zoffset2 = npRadius*std::sin(phiAngle + deltaPhi);

             double xoffset = k*dl * std::cos(screenRotate) - npRadius* std::cos(phiAngle) * std::sin(screenRotate) ;
             double xoffset2 =  k*dl * std::cos(screenRotate) - npRadius* std::cos(phiAngle+deltaPhi) * std::sin(screenRotate) ;



             // qDebug() << k << " " <<  k*cylinderHalfLength/(numHalfCylinderSegments) * std::cos(screenRotate);
             //  double ycentral = npRadius * std::cos(screenRotate)*std::cos(phiAngle) + k*cylinderHalfLength*std::sin(screenRotate);
              double ycentral = 0.5* (npRadius*std::cos(screenRotate)*( std::cos(phiAngle) + std::cos(phiAngle+deltaPhi)  ) +(dl+k*dl)*std::sin(screenRotate)   );

              double x0 = 0;
              double z0 = zMax- npRadius ;
              //qDebug() << phiAngle << " " << deltaPhi <<" " << xoffset << " " << xoffset2<<" " << zoffset << " " << zoffset2;
         std::vector<double> plotRect1{   xoffset  ,ycentral, z0 + zoffset +npRadius, (xoffset2-xoffset), -1, 1, (zoffset2-zoffset), 1} ;
            plotObjects.emplace_back(plotRect1);
        }


         }
        std::vector<double> plotCircle2{0,0,zMax- npRadius ,npRadius, -1, 1,  cylinderFaceWidth/2.0, 0} ;


  //   plotObjects.emplace_back(  plotCircle2 );



    }
    else if(npShapeIndex == 2){ //make a cube from -R to +R
        int numSegments = 4;
        for(int i = -numSegments; i<numSegments; i++){
           double iFrac = (double)i / (double)numSegments;
             double i1Frac = (i + 1.0) / (double)numSegments;
            double x0 =  npRadius*iFrac;
            double x1 = npRadius*i1Frac ;
            double y0 = npRadius;
            double y1 = npRadius;
            double x0p = std::cos(screenRotate)*x0 - std::sin(screenRotate)*y0 ;
            double x1p = std::cos(screenRotate)*x1 - std::sin(screenRotate)*y1 ;
            double y0p = std::sin(screenRotate)*x0 + std::cos(screenRotate)*y0 ;
            double y1p = std::sin(screenRotate)*x1 + std::cos(screenRotate)*y1 ;
            double yAv = 0.5*(y0p+y1p);
            double z0 = zMax - npRadius;
            double faceColourVal = 160 + 20*yAv/npRadius;

        std::vector<double> plotRect1{   x0p  ,yAv, z0, x1p, -faceColourVal, 1, z0+npRadius*2, 1} ;
         plotObjects.emplace_back(plotRect1);

         y0 = -npRadius;
         y1 = -npRadius;
           x0p = std::cos(screenRotate)*x0 - std::sin(screenRotate)*y0 ;
           x1p = std::cos(screenRotate)*x1 - std::sin(screenRotate)*y1 ;
           y0p = std::sin(screenRotate)*x0 + std::cos(screenRotate)*y0 ;
           y1p = std::sin(screenRotate)*x1 + std::cos(screenRotate)*y1 ;
           yAv = 0.5*(y0p+y1p);

         double faceColourVal2 = 160 + 20*yAv/npRadius;
     std::vector<double> plotRect2{   x0p  ,yAv, z0, x1p, -faceColourVal2, 1, z0+npRadius*2, 1} ;
      plotObjects.emplace_back(plotRect2);


      y0 = npRadius*iFrac;
      y1 = npRadius*i1Frac ;
      x0 = -npRadius;
      x1 = -npRadius;

        x0p = std::cos(screenRotate)*x0 - std::sin(screenRotate)*y0 ;
        x1p = std::cos(screenRotate)*x1 - std::sin(screenRotate)*y1 ;
        y0p = std::sin(screenRotate)*x0 + std::cos(screenRotate)*y0 ;
         y1p = std::sin(screenRotate)*x1 + std::cos(screenRotate)*y1 ;
       yAv = 0.5*(y0p+y1p);
      double faceColourVal3 = 160 + 20*yAv/npRadius;
  std::vector<double> plotRect3{   x0p  ,yAv, z0, x1p, -faceColourVal3, 1, z0+npRadius*2, 1} ;
   plotObjects.emplace_back(plotRect3);




   y0 = npRadius*iFrac;
   y1 = npRadius*i1Frac ;
   x0 = npRadius;
   x1 = npRadius;

     x0p = std::cos(screenRotate)*x0 - std::sin(screenRotate)*y0 ;
     x1p = std::cos(screenRotate)*x1 - std::sin(screenRotate)*y1 ;
     y0p = std::sin(screenRotate)*x0 + std::cos(screenRotate)*y0 ;
      y1p = std::sin(screenRotate)*x1 + std::cos(screenRotate)*y1 ;
    yAv = 0.5*(y0p+y1p);
   double faceColourVal4 = 160 + 20*yAv/npRadius;

std::vector<double> plotRect4{   x0p  ,yAv, z0, x1p, -faceColourVal4, 1, z0+npRadius*2, 1} ;
plotObjects.emplace_back(plotRect4);



        }




    }
    else{
        std::vector<double> plotCircle{ 0,0,zMax - npRadius,npRadius, -1, 1, npRadius, 0} ;
    plotObjects.emplace_back(  plotCircle );
    for(int i =1; i<maxNPCircles; i++){
    double npradiusSmall = npRadius  * std::cos(0.5 *  fPi * i/maxNPCircles )  ;
    double colourParamLocal = -1 -i;
    std::vector<double> plotCircle2{ 0,sqrt(npRadius*npRadius - npradiusSmall*npradiusSmall),zMax - npradiusSmall,npradiusSmall, colourParamLocal, 1, npradiusSmall, 0} ;
    plotObjects.emplace_back(  plotCircle2 );
    }

    }


 //sort by y-coord to get an approximation of depth
    std::sort(plotObjects.begin(), plotObjects.end(),
              [](const std::vector<double>& a, const std::vector<double>& b) {
      return a[1] < b[1];
    });


QPen outlinePen(Qt::black);
QPen outlineRed(Qt::red);
QPen outlineBlue(Qt::blue);
QPen currentOutline;

QBrush whiteFill(Qt::white);
QBrush greyFill(Qt::gray);

extent = std::max(1.0, zMax);
extent = std::max(extent, xAbsMax);
extent = std::max(extent, yAbsMax);

double scaleFactor = 0.5*scaleBasis/(extent + 0.01);
 double boundSize = scaleFactor*( 1 + 1.5*proteinRadius);
outlinePen.setWidth( boundSize/100);

double npLeft = -   npRadius;
double npUp = zMax -( 0 ) - npRadius;
bool plottedNP = false;

int alphaVal = this->findChild<QSlider *>("opacitySlider")->value() ;

double cylinderFaceLength =    std::cos(screenRotate)  * cylinderHalfLength/numHalfCylinderSegments;


for(const auto& pc: plotObjects){
 //qDebug() << " plotting bead at max-y-value " << pc[2] << "\n";
 double xleft = pc[0] - pc[3];
 double zup =  zMax -pc[2]  - pc[3];


 double colourVal = 0;

 //colour parameter is greater than zero: plot a regular bead
 if(pc[4] > 0){
 colourVal = ( colourParamMax - pc[4]   )/(0.01 +  colourParamMax - colourParamMin ) ;
 QBrush colourFill(  QColor(255*colourVal,0,0,alphaVal)  );
 QPen colourPen(  QColor(255*colourVal,0,0)) ;
 if(pc[5] > 0.5){
     currentOutline = outlinePen;
     if( this->findChild<QCheckBox *>("showChargeBox")->isChecked() == true  ){
     if(pc[6] > 0.5){
         currentOutline = outlineRed;
     }
     else if(pc[6] < -0.5){
         currentOutline = outlineBlue;
     }
     }
  pdbScene.addEllipse(xleft*scaleFactor, zup*scaleFactor, 2*pc[3]*scaleFactor , 2*pc[3]*scaleFactor  ,currentOutline , colourFill);


 }
 }
else{
     //interpret as a special NP plot item with half-width pc[6]; note that these interpret the pc[2] argument as the upper corner and not the centre.
     zup = pc[2];
     xleft = pc[0] - pc[6];
     // colourparam = -1 should give 160


     double colourValNP =  160 + 60 * std::sin( 0.5* (-1   -pc[4] )*fPi/maxNPCircles  );
     QBrush colourFillNP(  QColor(colourValNP,colourValNP,colourValNP)  );
     QPen colourOutlineNP(QColor(colourValNP,colourValNP,colourValNP));



     if(pc[4]>-1.1){
         colourOutlineNP = outlinePen;
     }
     if(pc[7] < 0.5){
      pdbScene.addEllipse(xleft*scaleFactor, zup*scaleFactor, 2*pc[6]*scaleFactor , 2*pc[3]*scaleFactor  ,colourOutlineNP ,colourFillNP);
     }
     else{
         //add a cylinder face  xleft*scaleFactor, zup*scaleFactor,  pc[3] stores the x offset for the lower points, pc[6] stores the z offset for lower points
         if(npShapeIndex==2){
          QPolygonF face;
          double colourValNP =  -pc[4];
          QBrush colourFillNP(  QColor(colourValNP,colourValNP,colourValNP)  );
          QPen colourOutlineNP(QColor(colourValNP,colourValNP,colourValNP));

          face << QPointF(pc[0]*scaleFactor, pc[2]*scaleFactor) << QPointF(pc[3]*scaleFactor, pc[2]*scaleFactor) << QPointF(pc[3]*scaleFactor, pc[6]*scaleFactor)<< QPointF(pc[0]*scaleFactor, pc[6]*scaleFactor) ;
          pdbScene.addPolygon(face, colourOutlineNP,colourFillNP);
         }
         else{
         QPolygonF cylinderFace;
         xleft = pc[0] - cylinderFaceLength;
         zup = pc[2];
         double cdx = pc[3];
         double cdz = pc[6];
         cylinderFace << QPointF(xleft*scaleFactor   , zup*scaleFactor) << QPointF( (xleft+cylinderFaceLength)*scaleFactor, zup*scaleFactor  ) << QPointF( (xleft+cylinderFaceLength+cdx)*scaleFactor, (zup+cdz)*scaleFactor ) << QPointF((xleft+cdx)*scaleFactor, (zup+cdz)*scaleFactor) ;
         pdbScene.addPolygon( cylinderFace ,colourOutlineNP ,colourFillNP);
         }
     }

 }

/*
 if(pc[1] > 0 && plottedNP == false){

     if(npShapeIndex == 1){ //plot a cylinder
         double cylinderFaceWidth = std::abs( 2*npRadius * std::sin(screenRotate) ) ;
        // npLeft = -cylinderFaceWidth;
         double cylinderLeft = -cylinderFaceWidth/2.0;

         int lowerSign = -1;

         if(screenRotate > 3.1415){
             lowerSign = 1;
         }

         pdbScene.addEllipse(  (cylinderLeft + lowerSign * cylinderHalfLength*std::cos(screenRotate))*scaleFactor, npUp*scaleFactor, cylinderFaceWidth*scaleFactor , 2*npRadius*scaleFactor  ,outlinePen ,greyFill);
         pdbScene.addEllipse(  cylinderLeft*scaleFactor, npUp*scaleFactor, cylinderFaceWidth*scaleFactor , 2*npRadius*scaleFactor  ,outlinePen ,greyFill);
         pdbScene.addEllipse(  (cylinderLeft - lowerSign * cylinderHalfLength*std::cos(screenRotate))*scaleFactor, npUp*scaleFactor, cylinderFaceWidth*scaleFactor , 2*npRadius*scaleFactor  ,outlinePen ,greyFill);
     }
     else{
     pdbScene.addEllipse( npLeft*scaleFactor, npUp*scaleFactor, 2*npRadius*scaleFactor , 2*npRadius*scaleFactor  ,outlinePen ,greyFill);
     }
     plottedNP = true;
 }
*/



}

//replot the NP's outline dashed
QPen dashedPen(Qt::blue);
QVector<qreal> dashes;
dashes << 4 << 4;
dashedPen.setDashPattern(dashes );
QBrush noFillBrush(Qt::white);
QBrush blueFill(Qt::blue);
QPen bluePen(Qt::blue);
QBrush greenFill(Qt::green);
QPen greenPen(Qt::green);
noFillBrush.setStyle(Qt::BrushStyle( Qt::NoBrush));
if(npShapeIndex == 1){
    double cylinderOutLineLeft = - cylinderHalfLength * std::abs( std::cos(screenRotate) ) ;
    pdbScene.addRect( cylinderOutLineLeft*scaleFactor, npUp*scaleFactor, -2*cylinderOutLineLeft*scaleFactor, 2*npRadius*scaleFactor  ,dashedPen, noFillBrush );
}
else if(npShapeIndex == 2){
    pdbScene.addRect(   -npRadius*scaleFactor, npUp*scaleFactor, 2*npRadius*scaleFactor, 2*npRadius*scaleFactor  ,dashedPen, noFillBrush  );
}
else{
pdbScene.addEllipse( npLeft*scaleFactor, npUp*scaleFactor, 2*npRadius*scaleFactor , 2*npRadius*scaleFactor  ,dashedPen, noFillBrush );
//add extra wireframe
//double wfsPhi = 45.0 * 3.1414 / 180.0;
//double effectiveRadius = npRadius * std::cos(wfsPhi);
//pdbScene.addEllipse( -effectiveRadius*scaleFactor, npUp*scaleFactor, 2*effectiveRadius*scaleFactor , 2*npRadius*scaleFactor  ,dashedPen, noFillBrush );

}


//plot COMs

double proteinCOMSize = 0.1;
double xleft = 0 - proteinCOMSize;
double zup = zMax - proteinCOMSize;

if( this->findChild<QCheckBox *>("showCOMBox")->isChecked()  == true ){
pdbScene.addEllipse(xleft*scaleFactor, (-1.0*plotComVal+zup)*scaleFactor, 2*proteinCOMSize*scaleFactor , 2*proteinCOMSize*scaleFactor  ,greenPen , greenFill) ;
pdbScene.addEllipse( (-0.1)*scaleFactor, (zMax-0.1)*scaleFactor, 2*0.1*scaleFactor , 2*0.1*scaleFactor  ,bluePen, blueFill );
}

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
    loadMaterials( baseMaterialSet,false);
  }

  QString baseBeadSet = QDir::cleanPath(  uaGlobalPath+"/pmfp-beadsetdef/PMFP-BeadSet.csv" );
 if(QFile::exists(baseBeadSet ) ){
   loadBeadSetFile( baseBeadSet );
 }

}

void MainWindow::on_colourBoltz_clicked()
{
    showEnergyShrinkWrap = false;
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



    double currentPhi = (  i*5 + 2.5)* fPi/180.0;
    double currentTheta = (j*5 + 2.5)* fPi/180.0;
    double energyVal = energyData[i][j];
    double p = -1.0*currentPhi;
    double t = fPi - currentTheta;
    double omega = this->findChild<QDial *>("omegaDial")->value() * fPi/180.0;
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
    double minEnergy =0;
    double maxEnergy = 0.1;
    for(int i=0;i<72;++i){
        for(int j=0;j<36;++j){
            minEnergy =std::min( energyData[i][j], minEnergy);
            maxEnergy = std::max(energyData[i][j],maxEnergy);
        }
    }

    showEnergyShrinkWrap = true;
    for( auto& atom: atomList){
        int i =(int)( floor( (atom.orientationPhi*180.0/fPi )/5.0));
        int j = (int)( floor((atom.orientationTheta*180.0/fPi)/5.0));
        i = std::max(0,i);
        j = std::max(0,j);
        i = std::min(71,i);
        j = std::min(35,j);
        //qDebug() << QString::number(atom.orientationPhi*180.0/3.1415) << " mapped to " << QString::number(i) << "\n";
        atom.colourParam = 1.0 - (energyData[i][j] -maxEnergy)/(minEnergy - maxEnergy);

            qDebug() << QString::number(energyData[i][j]) << " " <<QString::number(atom.colourParam) << "\n";

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


void MainWindow::on_npcpModeBox_stateChanged(int arg1)
{
     // qDebug() << "Checkbox new state: " << QString::number(arg1) << "\n";

    if (this->findChild<QCheckBox *>("npcpModeBox")->isChecked() == true){
      //swap text over as needed
         this->findChild<QPushButton *>("pdbTargetButton")->setText("Medium file");
        this->findChild<QComboBox *>("npcpModeOptions")->setEnabled(true) ;
        this->findChild<QLineEdit *>("resultFolderBox")->setDisabled(true);

    }
    else{
    this->findChild<QPushButton *>("pdbTargetButton")->setText("PDB Target");
this->findChild<QComboBox *>("npcpModeOptions")->setEnabled(false) ;
         this->findChild<QLineEdit *>("resultFolderBox")->setDisabled(false);
    }


}


void MainWindow::on_mediumEditTable_customContextMenuRequested(const QPoint &pos)
{
    //QModelIndex index = this->findChild<QTableWidget *>("mediumEditTable")->indexAt(pos);
    tableMenu.popup( this->findChild<QTableWidget *>("mediumEditTable")->viewport()->mapToGlobal(pos)    ) ;
    //qDebug() << index << "\n";
}


void MainWindow::addMoleculeToMedium(){
QTableWidget *tableWidget = this->findChild<QTableWidget *>("mediumEditTable") ;
  blockMediumColouring=true;
        tableWidget->insertRow ( tableWidget->rowCount() );
        tableWidget->setItem   ( tableWidget->rowCount()-1,      0,  new QTableWidgetItem("NewMolecule"));
                tableWidget->setItem   ( tableWidget->rowCount()-1,      1,  new QTableWidgetItem("0.0"));
                  blockMediumColouring=false;

    colourStructureRow(tableWidget->rowCount()-1)          ;

}
void MainWindow::removeMoleculeFromMedium(){

     QTableWidget *tableWidget = this->findChild<QTableWidget *>("mediumEditTable") ;
   // int currentRow =this->findChild<QTableWidget *>("mediumEditTable")->currentRow();
     int currentRow = tableWidget->currentRow();
  blockMediumColouring=true;
   // qDebug() << currentRow ;
   // qDebug() << tableWidget->rowCount();
    if(currentRow > -1 ){
      //qDebug() << QString(currentRow) << " can be removed";
      tableWidget->removeRow(currentRow);
    }
      blockMediumColouring=false;


}

void MainWindow::on_mediumNewButton_clicked()
{
    QTableWidget *tableWidget = this->findChild<QTableWidget *>("mediumEditTable") ;
    blockMediumColouring=true;
    tableWidget->setRowCount(0);
blockMediumColouring=false;
}


void MainWindow::on_mediumSaveButton_clicked()
{
    QString outputFile =   QFileDialog::getSaveFileName(this, tr("Medium CSV File"),uaGlobalPath ,  tr("Comma separated variable file (*.csv) "));
if( !outputFile.endsWith(".csv")){
 outputFile.append(".csv");
}
    QFile file(outputFile);
    QTableWidget *tableWidget = this->findChild<QTableWidget *>("mediumEditTable") ;


       if(file.open(QIODevice::WriteOnly | QIODevice::Text)) {
           QTextStream fileOut(&file);
           fileOut << "#MoleculeID,Concentration[mol/L]\n";
          int numRows =  tableWidget->rowCount();
         for(int i=0; i<numRows; ++i){


             fileOut << tableWidget->item(i,0)->text() << "," << tableWidget->item(i,1)->text() << "\n";
         }

       }
    file.close();
}


void MainWindow::on_mediumLoadButton_clicked()
{
    QString targetMediumFile = QFileDialog::getOpenFileName(this, tr("Medium File"),  this->uaGlobalPath,  tr("CSV (*.csv)"));
     QTableWidget *tableWidget = this->findChild<QTableWidget *>("mediumEditTable") ;
     blockMediumColouring=true;
     tableWidget->setRowCount(0);
    QFile file(targetMediumFile);
    if(file.open(QIODevice::ReadOnly | QIODevice::Text)) {
    QTextStream fileIn(&file);

    while(!fileIn.atEnd()){
        std::string lineIn = fileIn.readLine().toStdString();


        if(lineIn.substr(0, 1) != "#") {
            std::vector<std::string> results;
             boost::split(results, lineIn, [](char c){return c == ',';});
             //silicaquartz,surface/SiO2-Quartz,hamaker/SiO2_Quartz.dat,1
             //qDebug() << QString::fromStdString(results[0]) << " " << QString::fromStdString(results[1]) << " " <<  QString::fromStdString(results[2] )<< "\n";
             if(results.size()==2){
            // materialTypes.emplace_back( MaterialType( QString::fromStdString(results[0])   , QString::fromStdString(results[1]) ,QString::fromStdString(results[2])     ) ) ;
             //addBeadTypeMaterialBox->addItem( QString::fromStdString(results[0]) ,   QList<QVariant>() <<  QString::fromStdString(results[1]) <<  QString::fromStdString(results[2])) ;
             //this->findChild<QComboBox *>("materialDropdown")->addItem( QString::fromStdString(results[0]) ,   QList<QVariant>() <<  QString::fromStdString(results[1]) <<  QString::fromStdString(results[2]))  ;

                 tableWidget->insertRow ( tableWidget->rowCount() );
                 tableWidget->setItem   ( tableWidget->rowCount()-1,      0,  new QTableWidgetItem(QString::fromStdString(results[0])));
                 tableWidget->setItem   ( tableWidget->rowCount()-1,      1,  new QTableWidgetItem(QString::fromStdString(results[1])));


             }

        }
    }


}
colourStructures();
blockMediumColouring=false;
}

void MainWindow::on_cancelRunButton_clicked()
{

 processHandle.kill();
    this->findChild<QPushButton *>("runUAButton")->setDisabled(false);
   this->findChild<QPushButton *>("cancelRunButton")->setDisabled(true);

}


void MainWindow::on_checkStructureButton_clicked()
{

    QTableWidget *tableWidget = this->findChild<QTableWidget *>("mediumEditTable") ;
    int numRows =  tableWidget->rowCount();
    for(int i=0; i<numRows; ++i){

        QString pdbName = tableWidget->item(i,0)->text() ;

        bool isPDBDB = false;
        bool isAFDB = false;

        if(pdbName.left(4) =="PDB-"  ){
        isPDBDB = true ;
        }
        if(pdbName.left(5) == "AFDB-"){
        isAFDB = true ;
        }
         if(  !QFile::exists( uaGlobalPath+"/all_proteins/"+pdbName+".pdb")){


             if( isPDBDB == true){
                 //do PDB testing



                 auto doPDBTest = QMessageBox::question(this, "Find PDB structure?", "No structure found for "+pdbName+". Check PDB?");

                 if(doPDBTest == QMessageBox::Yes){
                     //try to get from PDB, if so set foundStructure = true

                     QUrl afPath("https://files.rcsb.org/download/" +pdbName.mid(4)+ ".pdb");
                      MainWindow::startDownload( afPath) ;

                 }



             }
             else if(isAFDB == true){
                   auto doAFTest = QMessageBox::question(this, "Find AlphaFold structure?", "No structure found for "+pdbName+". Check AlphaFoldDB?");
                   if(doAFTest == QMessageBox::Yes){
                   //try to get from alphafold, if so set foundStructure = true
                   QUrl afPath("https://alphafold.ebi.ac.uk/files/AF-" +pdbName.mid(5)+ "-F1-model_v4.pdb");
                    MainWindow::startDownload( afPath) ;
                    }
             }






         }

       colourStructureRow(i);
    }


       // colourStructures();

}

void MainWindow::startDownload(QUrl targetURL){


networkManager.get(QNetworkRequest(targetURL));
}


void MainWindow::downloadReplyFinished(QNetworkReply *reply){
//make all_proteins if needed

    QString targetFile = reply->url().fileName();

    QDir baseDir = QDir(uaGlobalPath);
    baseDir.mkpath( "all_proteins");

    if(reply->error()){
    QMessageBox::warning(this, tr("NPCoronaPredict-GUI"),  "Failed to find structure for: "+targetFile+" \n"   );
    }
    else{


    if( targetFile.size()!= 8  ){   //assume XXXX.pdb is from PDB and anything else is from AF 
        targetFile = "AFDB-"+ targetFile.split("-")[1]+".pdb";   //set AFDB name
    }
    else{
       targetFile= "PDB-"+targetFile;
    }
    //save to file
    QFile *saveFile = new QFile(uaGlobalPath+"/all_proteins/"+targetFile);
    if(saveFile->open(QFile::Append)){
        saveFile->write(reply->readAll());
        saveFile->flush();
        saveFile->close();
    }
      colourStructures();
      QMessageBox::warning(this, tr("NPCoronaPredict-GUI"),  "Found structure for: "+targetFile+" \n"   );
   delete saveFile;


}


    reply->deleteLater();
}


void MainWindow::colourStructures()
{
        QTableWidget *tableWidget = this->findChild<QTableWidget *>("mediumEditTable") ;
        int numRows =  tableWidget->rowCount();
        for(int i=0; i<numRows; ++i){
           colourStructureRow(i);
        }

}


void MainWindow::colourStructureRow(int row)
{
        QTableWidget *tableWidget = this->findChild<QTableWidget *>("mediumEditTable") ;


            QString pdbName = tableWidget->item(row,0)->text() ;
             if(  QFile::exists( uaGlobalPath+"/all_proteins/"+pdbName+".pdb")){
                 //qDebug() << " found pdb \n";
                 tableWidget->item(row,0)->setData(Qt::BackgroundRole, QColor(175,255,175));
                 tableWidget->item(row,1)->setData(Qt::BackgroundRole, QColor(175,255,175));
             }
             else{
                 //qDebug() << " did not find pdb \n";
                 tableWidget->item(row,0)->setData(Qt::BackgroundRole, QColor(255,175,175));
                 tableWidget->item(row,1)->setData(Qt::BackgroundRole, QColor(255,175,175));

             }


}


void MainWindow::on_mediumEditTable_cellChanged(int row, int column)
{
    if(blockMediumColouring==false){
         colourStructureRow(row);
    }

}

void MainWindow::on_mediumEditTable_currentCellChanged(int currentRow, int currentColumn, int previousRow, int previousColumn)
{
    if(blockMediumColouring==false){
         colourStructureRow(currentRow);
    }
}


void MainWindow::on_showCOMBox_stateChanged(int arg1)
{
    updateMoleculeBox();
}


void MainWindow::on_opacitySlider_sliderMoved(int position)
{
  //  updateMoleculeBox();
}


void MainWindow::on_opacitySlider_valueChanged(int value)
{
    updateMoleculeBox();
}


void MainWindow::on_loadBeadmapButton_clicked()
{
    QString targetBeadSetFile = QFileDialog::getOpenFileName(this, tr("Bead definition file"),  this->uaGlobalPath , tr("CSV-file (*.csv)"));
    if(targetBeadSetFile!=""){
        loadBeadSetFile(targetBeadSetFile);
    }

}

void MainWindow::loadBeadSetFile(QString targetFile){
    QFile file(targetFile);
    if(file.open(QIODevice::ReadOnly | QIODevice::Text)) {
    QTextStream fileIn(&file);
    while(!fileIn.atEnd()){
        std::string lineIn = fileIn.readLine().toStdString();
        if(lineIn.substr(0,1) != "#"   and lineIn.length()>5){
            std::vector<std::string> results;
             boost::split(results, lineIn, [](char c){return c == ',';}  , boost::token_compress_on   );
             if(results.size() == 3){
                 //qDebug() << "loaded " << QString::fromStdString(results[0] );
                 //BeadType newBead(  "ZZZ", 0.5, -1.0);
                 beadTypeMap[ results[0] ] = BeadType(   results[0],std::stod(results[2]), std::stod(results[1])) ;
             }
        }
    }
    }

    //loop over atoms and update

    for(auto& atom: atomList){
        auto beadSearch = beadTypeMap.find(atom.atomName);
        double radiusVal = 0.5;
        double chargeVal = 0.0;
        if(  beadSearch!=beadTypeMap.end()   ){
            radiusVal = beadSearch->second.radius;
            chargeVal = beadSearch->second.charge;
        }
        atom.radius = radiusVal;
        atom.charge = chargeVal;
    }

updateMoleculeBox();
}

void MainWindow::on_showChargeBox_stateChanged(int arg1)
{
    updateMoleculeBox();
}


void MainWindow::on_npShapeBox_currentIndexChanged(int index)
{
    updateMoleculeBox();
}


void MainWindow::on_materialDropdown_currentIndexChanged(int index)
{
    //make sure the selected NP shape is compatible with the chosen material - if the material is a CNT set the NP to the correct shape, if not default to sphere

    QString currentNPShape = this->findChild<QComboBox *>("npTargetShapeOverride")->currentText();

    QString targetMaterial = this->findChild<QComboBox *>("materialDropdown")->currentText() ;

    if( targetMaterial.contains("mwcnt")){
        this->findChild<QComboBox *>("npTargetShapeOverride")->setCurrentIndex(4);
    }
    else if(targetMaterial.contains("cnt")){
        this->findChild<QComboBox *>("npTargetShapeOverride")->setCurrentIndex(3);
    }
    else{
        if(currentNPShape=="MWCNT" || currentNPShape=="SWCNT"){
            this->findChild<QComboBox *>("npTargetShapeOverride")->setCurrentIndex(0);
        }
    }

}


void MainWindow::on_npTargetShapeOverride_currentIndexChanged(int index)
{
    if(index == 0){
        this->findChild<QComboBox *>("npShapeBox")->setCurrentIndex(  0) ; //set view to sphere
    }
    else if(index == 1 || index == 3 || index == 4){
        this->findChild<QComboBox *>("npShapeBox")->setCurrentIndex(  1); //set view to cylinder
    }
    else if(index == 2){
        this->findChild<QComboBox *>("npShapeBox")->setCurrentIndex(  2) ; //set view to cube
    }
    else{
       this->findChild<QComboBox *>("npShapeBox")->setCurrentIndex(  0) ; //set view to sphere??
    }
}

