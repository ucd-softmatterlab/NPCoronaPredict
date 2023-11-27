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

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    //initialise vectors for storage of bead types and added beads
    std::default_random_engine generator;
    addBead = new AddBead(this);
    addBeadType = new AddBeadType(this);
    addBrush = new AddBrush(this);
        addShell = new AddShell(this);
        tipsWindowI = new tipsWindow(this);
        //scene = new QGraphicsScene();

        this->findChild<QTableWidget *>("beadTypeTable")->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
        this->findChild<QTableWidget *>("beadTable")->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);


    connect(addBead, &AddBead::sendNewNPBead,this , &MainWindow::recieveNewNPBead);
connect(addBeadType, &AddBeadType::sendNewNPBeadType,this , &MainWindow::recieveNewNPBeadTypeAuto);
connect(addShell, &AddShell::sendNewShell,this , &MainWindow::recieveNewShell);
connect(addBrush, &AddBrush::sendNewBrush,this , &MainWindow::recieveNewBrush);
QString currentPath = QDir::currentPath() ;
currentPath = QDir::cleanPath(currentPath+"/../");

QAction *removeSelectedBeads = new QAction("Remove multiple", this);
npBeadTableMenu.addAction( removeSelectedBeads);
connect(removeSelectedBeads,SIGNAL(triggered()), this, SLOT( removeMultipleNPBeads()));

this->uaGlobalPath = currentPath;
this->findChild<QLineEdit *>("uaPath")->setText(this->uaGlobalPath);



this->materialTypes.emplace_back( MaterialType( "Custom"   , "",""     ) ) ;
addBeadType->findChild<QComboBox *>("materialTypeBox")->addItem( "Custom" ,   QList<QVariant>() <<  QString::fromStdString("") <<  QString::fromStdString("")) ;
addShell->findChild<QComboBox *>("materialTypeBox")->addItem( "Custom" ,   QList<QVariant>() <<  QString::fromStdString("") <<  QString::fromStdString("")) ;

 }

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_newBeadTypeButton_clicked()
{
    //open newBeadType dialogue
    addBeadType->searchPath = uaGlobalPath;



    addBeadType->exec();
}


void MainWindow::on_newBeadButton_clicked()
{

    QComboBox *brushComboTypes = addBead->findChild<QComboBox *>("beadTypeID") ;
    brushComboTypes->clear();
    for(int i = 0; i < (int)beadTypes.size(); ++i){
     brushComboTypes->addItem( QString::number(i) ) ;
    }

    addBead->exec();
}

void MainWindow::on_newBrushButton_clicked()
{
    QComboBox *brushComboTypes = addBrush->findChild<QComboBox *>("brushBeadID") ;
    brushComboTypes->clear();


         //   qDebug() << "Updating bead types \n";
    for(int i = 0; i < (int)beadTypes.size(); ++i){
     //   qDebug() << i << "\n";
     brushComboTypes->addItem( QString::number(i) ) ;
//qDebug() << i << "\n";
    }
   // qDebug() << "Clearing existing bead radii \n";
            addBrush->beadRadii.erase(addBrush->beadRadii.begin(), addBrush->beadRadii.end()) ;
      //  qDebug() << "Updating  bead radii \n";
     for( const auto& np: beadTypes){
        //   qDebug() << "Adding bead radii"   << np.radius  << " done \n";
     addBrush->beadRadii.emplace_back( np.radius   ) ;

    }
    addBrush->currentOuterRadius = npOuterRadius;
    addBrush->updateDensityBox();
    addBrush->exec();
}



void MainWindow::on_newShellButton_clicked()
{
    addShell->searchPath = uaGlobalPath;
    addShell->exec();
}

void MainWindow::recieveNewNPBead( double x, double y, double z, int beadTypeID, bool doUpdate=false){
    NPBead newNPBead(beadTypeID,x,y,z);
    npBeads.emplace_back(newNPBead);
    //qDebug() << "New bead: " << x << " " << y << " " << z << " type: " << beadTypeID << "\n";
    //qDebug() << "All beads: \n";
    //for( const auto& np: npBeads){
    //qDebug() <<  np.x << " " << np.y << " " << np.z << " type: " << np.beadTypeID << "\n";
    //}

    if(doUpdate==true){
    updateBeadTable() ;
    }


}

void MainWindow::recieveNewNPBeadTypeAuto( std::string hamakerFileIn, std::string surfaceDirIn, float   radiusIn, float   surfacePotentialIn, float surfFactorIn, float coreFactorIn, float ljCutoffIn, int correctionOverrideIn){
    recieveNewNPBeadType(hamakerFileIn,surfaceDirIn,radiusIn,surfacePotentialIn,surfFactorIn,coreFactorIn,ljCutoffIn,correctionOverrideIn,-1,true);
}

void MainWindow::recieveNewNPBeadType( std::string hamakerFileIn, std::string surfaceDirIn, float   radiusIn, float   surfacePotentialIn, float surfFactorIn, float coreFactorIn, float ljCutoffIn, int correctionOverrideIn, int beadTypeIDIn, bool doUpdate){
    int newBeadID;
    if(beadTypeIDIn == -1){
    newBeadID = beadTypes.size();
    }
    else{
        newBeadID = beadTypeIDIn;
    }
    QDir uaDir(uaGlobalPath);
    /*
    std::string surfaceRelPath  = uaDir.relativeFilePath(QString::fromStdString(surfaceDirIn)).toStdString()  ;
    std::string hamakerRelPath = uaDir.relativeFilePath(QString::fromStdString(hamakerFileIn)).toStdString()  ;

    */
    std::string surfaceRelPath  =  surfaceDirIn  ;
    std::string hamakerRelPath = hamakerFileIn  ;

    BeadType newNPBead(newBeadID  , hamakerRelPath, surfaceRelPath, radiusIn,  surfacePotentialIn, surfFactorIn, coreFactorIn, ljCutoffIn, correctionOverrideIn   );
    beadTypes.emplace_back(newNPBead);
    if(doUpdate==true){
    updateBeadTypeTable( ) ;
    }
}

void MainWindow::recieveNewShell(std::string hamakerFileIn, std::string surfaceDirIn, float   innerRadius, float outerRadius, float   surfacePotentialIn,  double ljCutoff){
     int newBeadTypeID = beadTypes.size();
     recieveNewNPBeadTypeAuto( hamakerFileIn , surfaceDirIn , outerRadius, surfacePotentialIn, 1.0, 1.0, ljCutoff, 1   ) ;
     recieveNewNPBead(0, 0, 0  , newBeadTypeID);
     recieveNewNPBeadTypeAuto( hamakerFileIn , surfaceDirIn , innerRadius, 0.0, -1.0, -1.0, ljCutoff, 1  ) ;
     recieveNewNPBead(0, 0, 0  , newBeadTypeID+1);
}




double MainWindow::updateM(double mj, int numPoints){
    double piVal = 3.1415;
    //DLMF/Boost convention: ellipticE takes k as a parameter, k^2 = m
    //Mathematica/Scipy convention: ellipticE takes m as a parameter
    //Kaoy's algotihm uses  the mathematica/scipy convention for an argument m = -mj^2 where mj is real
    //thus for Boost input we have k^2 = - mj^2 => k = +- I mj

    //us

    //boost is defined only for -1 ... k ... 1 whereas mj can take arbitrary values, but because of the imaginary term this works out
    //EllipticE is ellint_2, ellipticK is ellint_1

    double ellipeVal = sqrt(1 + mj*mj) *boost::math::ellint_2( mj/sqrt(1 + mj*mj) ) ;
    double ellipkVal = boost::math::ellint_1(  mj/sqrt(1 + mj*mj) )/sqrt(1 + mj*mj) ;

    //double kVal = -m;
   // qDebug() << "input m: " << m << " k: " << kVal << "\n";
    //boost uses k = m^2, scipy = m . which here is -mInput*mInput.
    //double ellipeVal = scipyellipe( -mj*mj);//  ;// scspec.ellipe(-m*m);
    //double ellipkVal = scipyellipk( -mj*mj); // boost::math::ellint_1(kVal);//scspec.ellipk(-m*m);

    return (mj*piVal*numPoints*(2*ellipeVal - ellipkVal))/(numPoints*piVal*ellipeVal - numPoints*piVal*ellipkVal + mj * pow(ellipeVal,2)     );
}
double MainWindow::updateTheta(double theta, double m, double j){
    double piVal = 3.1415;
    return theta + ((2*j-1)*piVal - m*((2*j-1)*piVal/m))/(m * sqrt(  (1 + pow(m,2) *2* pow(sin(theta),2) ) ) );
}




void MainWindow::recieveNewBrush(double brushOccupancy,double brushRadialDist, int beadTypeID, bool forceAttach){
//qDebug() << "Starting brush: \n";
    double beadRadius = 1;

    //std::default_random_engine generator;
    std::uniform_real_distribution<double> uniformDist(0.0,1.0) ;
    //bool forceAttach = true;
    std::vector<NPBead> trialBeads;

    if ( std::find(knownBeadTypes.begin(), knownBeadTypes.end(), beadTypeID) != knownBeadTypes.end() ){
        //qDebug() << "Checking bead type " << beadTypeID <<  " in array of size " << beadTypes.size() <<  "\n";
        beadRadius = beadTypes[beadTypeID].radius ;
        //brushRadialDist is the bead centres, so this corresponds to npRadius+beadRadius in the Python version
        //Koay's algorithm
        int pointsNeeded = (int)( 4 *  brushRadialDist*brushRadialDist  / (1.1 * beadRadius*beadRadius) );
        //qDebug() << "Requires " << pointsNeeded << " beads \n";
        double mVal = sqrt( pointsNeeded * 3.1415);
        for(int j = 0; j<10; ++j){
            mVal = updateM(mVal, pointsNeeded);
          //  qDebug() << " updated m: " << mVal << "\n" ;
        }

        for(int i = 0; i < pointsNeeded; ++i){

        double thetaVal =acos( 1-  ( (1.0+2.0*i )/pointsNeeded   )  );
        //qDebug() << "theta value: " << thetaVal << " from " << 1-  ( (1.0+2.0*i )/pointsNeeded   ) << "\n";
        for( int j = 0; j< 21; ++j){
            thetaVal = updateTheta(thetaVal, mVal, 1);
        }

        double phiVal = thetaVal * mVal;


         double bx = brushRadialDist*cos(phiVal)*sin(thetaVal);
          double by = brushRadialDist*sin(phiVal)*sin(thetaVal);
           double bz = brushRadialDist*cos(thetaVal);

           if(forceAttach == false){
        if(   uniformDist(generator) < brushOccupancy){
         recieveNewNPBead( bx, by, bz, beadTypeID);
          }
         }
           else{
             trialBeads.emplace_back( NPBead( beadTypeID, bx, by, bz) ) ;
           }

         }

    if(forceAttach == true){
        //attach the trial beads to reach the target occupancy under the condition that they are in contact with another bead

       int totalNumBeadsAdded = 0;
       int numBeadsThisIteration = 0;
       int targetNumBeads = round( brushOccupancy * pointsNeeded );
       while(totalNumBeadsAdded < targetNumBeads){
           numBeadsThisIteration = 0;
       std::shuffle(trialBeads.begin(),trialBeads.end(), generator);
       std::vector<int> acceptedBeadIDs;
       double distanceEps = 0.2; //beads are touching if their surface to surface distance is less than this

       for(int i = 0; i < (int)trialBeads.size(); ++i){

           for(int j =0; j < (int)npBeads.size(); ++j){
           double ccd = sqrt( pow( trialBeads[i].x - npBeads[j].x,2) + pow( trialBeads[i].y - npBeads[j].y,2) +  pow( trialBeads[i].z - npBeads[j].z,2)  );
            if(ccd - beadRadius - beadTypes[npBeads[j].beadTypeID].radius  < distanceEps )
            {
                //accept bead
                acceptedBeadIDs.emplace_back(i);
                numBeadsThisIteration += 1;
                totalNumBeadsAdded += 1;
                //qDebug() << " bead " << i << " accepted due to contact with bead " << j << "\n";
                break;
            }
           }

           if(totalNumBeadsAdded >= targetNumBeads){
               break;
           }

       }
       //loop over all beads that were accepted, add these to the NP, remove from the list of candidates
       //reverse the ordering to make sure we can pop elements from the end without causing errors
       std::reverse(acceptedBeadIDs.begin(),acceptedBeadIDs.end());
       for(int i = 0; i<(int)acceptedBeadIDs.size(); ++i){
        int beadID = acceptedBeadIDs[i];
        NPBead acceptedBead = trialBeads[ beadID];
        recieveNewNPBead( acceptedBead.x, acceptedBead.y, acceptedBead.z, acceptedBead.beadTypeID) ;
        trialBeads.erase( trialBeads.begin() +  beadID );
       }

       // qDebug() << "Iteration complete. Beads accepted this pass: " << numBeadsThisIteration << " total beads accepted: " << totalNumBeadsAdded << "\n";
       if(numBeadsThisIteration==0){

           QMessageBox msgBox;
           msgBox.setText("No further beads can be added, check if a substrate exists.");
              msgBox.exec();
           break;
       }


       }

    }


    }
    else{
        qDebug() << "Bead ID not found, no beads will be generated";
        QMessageBox msgBox;
        msgBox.setText("Bead ID does not correspond to a known bead, skipping generation");
        msgBox.exec();
    }
updateBeadTable();
}


void MainWindow::updateBeadTypeTable( ){


    QTableWidget *beadTable = this->findChild<QTableWidget *>("beadTypeTable");
    beadTable->setRowCount(0);
    beadTable->setColumnCount(9);
    for( const auto& np: beadTypes){
    int currentRow = beadTable->rowCount();
    beadTable->insertRow( currentRow);
    beadTable->setItem(currentRow , 0, new QTableWidgetItem( QString::number(np.beadID ) ));

    beadTable->item(currentRow, 0)->setFlags(beadTable->item(currentRow, 0)->flags() &  ~Qt::ItemIsEditable);

    beadTable->setItem(currentRow , 1, new QTableWidgetItem( QString::fromStdString(np.surfaceDir)));
   beadTable->setItem(currentRow , 2, new QTableWidgetItem( QString::fromStdString(np.hamakerFile)));
    beadTable->setItem(currentRow , 3, new QTableWidgetItem( QString::number(np.radius  )));
    beadTable->setItem(currentRow , 4, new QTableWidgetItem( QString::number(np.surfacePotential  )));
     beadTable->setItem(currentRow ,5, new QTableWidgetItem( QString::number(np.surfFactor  )));
     beadTable->setItem(currentRow ,6, new QTableWidgetItem( QString::number(np.coreFactor  )));
     beadTable->setItem(currentRow ,7, new QTableWidgetItem( QString::number(np.ljCutoffVal  )));
     beadTable->setItem(currentRow ,8, new QTableWidgetItem( QString::number(np.correctionOverride  )));
    }

    if( beadTypes.size() > 0){
        //enable the bead and brush buttons
        this->findChild<QPushButton *>("newBeadButton")->setDisabled(false);
        this->findChild<QPushButton *>("newBrushButton")->setDisabled(false);
    }
    else{
        //disable bead and brush
        this->findChild<QPushButton *>("newBeadButton")->setDisabled(true);
        this->findChild<QPushButton *>("newBrushButton")->setDisabled(true);
    }
    rebuildBeadTypeList();
    updateBindingRadii();
    updateGraphicsWindow();
}

void MainWindow::updateBeadTable(){
    QTableWidget *beadTable = this->findChild<QTableWidget *>("beadTable");
    beadTable->setRowCount(0);
    beadTable->setColumnCount(5);
    for( const auto& np: npBeads){
    //qDebug() <<  np.beadTypeID << "\n";
    int currentRow = beadTable->rowCount();
    beadTable->insertRow( currentRow);
    beadTable->setItem(currentRow , 0, new QTableWidgetItem( QString::number(np.beadTypeID ) ));
    beadTable->setItem(currentRow , 1, new QTableWidgetItem( QString::number(np.x) ));
    beadTable->setItem(currentRow , 2, new QTableWidgetItem( QString::number(np.y ) ));
    beadTable->setItem(currentRow , 3, new QTableWidgetItem( QString::number(np.z  )));
    // beadTable->setItem(currentRow , 4, new QTableWidgetItem( QString::number(np.z  )));

    QTableWidgetItem *deleteBox = new QTableWidgetItem();
    deleteBox->setCheckState(Qt::Unchecked);
    beadTable->setItem(currentRow , 4, deleteBox );
    }

    updateBindingRadii();
    updateGraphicsWindow();

}


void MainWindow::rebuildBeadTypeList(){
    knownBeadTypes.erase( knownBeadTypes.begin(), knownBeadTypes.end());
    for(const auto& beadType : beadTypes){
        knownBeadTypes.push_back( beadType.beadID);
    }
}

void MainWindow::on_actionQuit_triggered()
{
   QCoreApplication::quit();
}


void MainWindow::on_findUADir_clicked()
{
    QString uaDir = QFileDialog::getExistingDirectory(this, tr("UA Install Directory"), this->uaGlobalPath, QFileDialog::ShowDirsOnly);
    if(uaDir != ""){
        this->findChild<QLineEdit *>("uaPath")->setText(uaDir);
        this->uaGlobalPath = uaDir;
    }
}


void MainWindow::saveNP( QString filename, bool doRotate = false, bool isPDBFile = false){

    double rotateMatrix[3][3];
    if(doRotate==true){
       // qDebug() << "Applying random rotation \n";


        std::uniform_real_distribution<double> uniformDist(0.0,1.0) ;

        //generate arvo matrix
         double piVal = 3.1415;
        double u1= uniformDist(generator);
        double u2= uniformDist(generator);
        double u3= uniformDist(generator);

        rotateMatrix[0][0] = -(cos(2*piVal*u1)*(1 - 2*u3*pow(cos(2*piVal*u2),2))) - 2*u3*cos(2*piVal*u2)*sin(2*piVal*u1)*sin(2*piVal*u2);
        rotateMatrix[0][1] = -((1 - 2*u3*pow(cos(2*piVal*u2),2))*sin(2*piVal*u1)) +    2*u3*cos(2*piVal*u1)*cos(2*piVal*u2)*sin(2*piVal*u2);
        rotateMatrix[0][2] = 2*sqrt(1 - u3)*sqrt(u3)*cos(2*piVal*u2);
        rotateMatrix[1][0] = 2*u3*cos(2*piVal*u1)*cos(2*piVal*u2)*sin(2*piVal*u2) + sin(2*piVal*u1)*(1 - 2*u3*pow(sin(2*piVal*u2),2));
        rotateMatrix[1][1] = 2*u3*cos(2*piVal*u2)*sin(2*piVal*u1)*sin(2*piVal*u2) - cos(2*piVal*u1)*(1 - 2*u3*pow(sin(2*piVal*u2),2));
        rotateMatrix[1][2] = 2*sqrt(1 - u3)*sqrt(u3)*sin(2*piVal*u2);
        rotateMatrix[2][0] = 2*sqrt(1 - u3)*sqrt(u3)*cos(2*piVal*u1)*cos(2*piVal*u2) -    2*sqrt(1 - u3)*sqrt(u3)*sin(2*piVal*u1)*sin(2*piVal*u2);
        rotateMatrix[2][1] = 2*sqrt(1 - u3)*sqrt(u3)*cos(2*piVal*u2)*sin(2*piVal*u1) +    2*sqrt(1 - u3)*sqrt(u3)*cos(2*piVal*u1)*sin(2*piVal*u2);
        rotateMatrix[2][2] = -1 + 2*(1 - u3);


    }
    else{
        rotateMatrix[0][0] = 1;
        rotateMatrix[0][1] = 0;
        rotateMatrix[0][2] = 0;
        rotateMatrix[1][0] = 0;
        rotateMatrix[1][1] = 1;
        rotateMatrix[1][2] = 0;
        rotateMatrix[2][0] = 0;
        rotateMatrix[2][1] = 0;
        rotateMatrix[2][2] = 1;
    }


    //open file

   QFile file(filename);



   if(file.open(QIODevice::WriteOnly | QIODevice::Text)) {
       QTextStream fileOut(&file);

       if(isPDBFile == false){
       fileOut <<  "#NPDesigner NP file \n" ;
       fileOut << "#NP-parameters \n";
       fileOut << "INNERBOUND," << npInnerRadius << "\n";
       fileOut << "OUTERBOUND," << npOuterRadius << "\n";
       fileOut << "#Bead definitions \n" ;
       fileOut << "#Radius,zeta,corefactor,surfacefactor,shape,hamakerFile,surfaceDir,pmfCutoff,correctionType \n";


       //write out NP bead types
       for( const auto& bt: beadTypes){
           QString beadTypeBaseStr="%1,%2,%3";
          // beadTypeBaseStr = beadTypeBaseStr.arg(bt.radius).arg(bt.surfacePotential).arg(bt.coreFactor).arg(bt.surfFactor).arg(1).arg(bt.hamakerFile).arg(bt.surfaceDir).arg(bt.ljCutoffVal).arg( bt.correctionOverride)  ;
        fileOut << "TYPE," << bt.radius<<"," << bt.surfacePotential/1000.0 << "," << bt.coreFactor << "," << bt.surfFactor << ",1,";
        fileOut <<      QString::fromStdString(bt.hamakerFile ) << "," << QString::fromStdString(bt.surfaceDir) << "," << bt.ljCutoffVal << "," << bt.correctionOverride << "\n";
       }

       fileOut << "#Beads\n" ;
       fileOut << "#type,x,y,z\n" ;

       //write out NP bead locations
       for(const auto& np: npBeads){

           double xOut = rotateMatrix[0][0] * np.x + rotateMatrix[0][1] * np.y + rotateMatrix[0][2] * np.z;
           double yOut = rotateMatrix[1][0] * np.x + rotateMatrix[1][1] * np.y + rotateMatrix[1][2] * np.z;
           double zOut = rotateMatrix[2][0] * np.x + rotateMatrix[2][1] * np.y + rotateMatrix[2][2] * np.z;



      fileOut << "BEAD," <<  np.beadTypeID << "," << xOut << "," << yOut <<  "," << zOut  <<  "\n"  ;


       }

       }
       else{
           fileOut << "TITLE NPDesigner NP file\n";
            int beadTypesSeen = 0;
           //write out NP bead types
           for( const auto& bt: beadTypes){
               QString beadTypeBaseStr="%1,%2,%3";
              // beadTypeBaseStr = beadTypeBaseStr.arg(bt.radius).arg(bt.surfacePotential).arg(bt.coreFactor).arg(bt.surfFactor).arg(1).arg(bt.hamakerFile).arg(bt.surfaceDir).arg(bt.ljCutoffVal).arg( bt.correctionOverride)  ;
            fileOut << "REMARK TYPE " << beadTypesSeen << " " << bt.radius<<"," << bt.surfacePotential << "," << bt.coreFactor << "," << bt.surfFactor << ",1,";
            beadTypesSeen++;
            fileOut <<      QString::fromStdString(bt.hamakerFile ) << "," << QString::fromStdString(bt.surfaceDir) << "," << bt.ljCutoffVal << "," << bt.correctionOverride << "\n";
           }

           int npBeadsSeen = 0;
           //1 = serial number = right justified three digits
           //ATOM  SSSSS  N   RRR ASSSS    XXXXXXXXYYYYYYYYZZZZZZZZOOOOOOBBBBBB           N
           QString PDBOutLine = "ATOM  %1  N   N%2 A%3    %4%5%6%7%8           N\n";
           for(const auto& np: npBeads){

               double xOut = rotateMatrix[0][0] * np.x + rotateMatrix[0][1] * np.y + rotateMatrix[0][2] * np.z;
               double yOut = rotateMatrix[1][0] * np.x + rotateMatrix[1][1] * np.y + rotateMatrix[1][2] * np.z;
               double zOut = rotateMatrix[2][0] * np.x + rotateMatrix[2][1] * np.y + rotateMatrix[2][2] * np.z;


          npBeadsSeen++;
          //fileOut << "ATOM" <<  np.beadTypeID << "," << xOut*10 << "," << yOut*10 <<  "," << zOut*10  <<  "\n"  ;
          fileOut<< PDBOutLine.arg(npBeadsSeen,  5).arg(np.beadTypeID,2,10,QLatin1Char('0')).arg(npBeadsSeen,4).arg( xOut*10, 8,  'f', 3).arg(yOut*10, 8,   'f', 3).arg(zOut*10, 8,   'f', 3).arg(1.00,  6, 'f', 2).arg(0.00,  6, 'f', 2) ;

       }


   }
   //close file
   file.close();

}


}

void MainWindow::on_saveNPButton_clicked()
{
    QString baseFileName = QFileDialog::getSaveFileName(this, tr("NP File"),uaGlobalPath ,  tr("NP-File (*.np);;PDB (*.pdb)"));



    //qDebug() << " Saving to " << baseFileName << "\n";
    bool isPDB = false;
    int chopChars = 3;
    QString outFileName = baseFileName;
    if( !outFileName.endsWith(".np")){

     if(outFileName.endsWith(".pdb")){
         isPDB = true;
         chopChars = 4;
     }
    else{
     outFileName.append(".np");
 }
    }

    saveNP( outFileName , false, isPDB);

    int numOrientations = this->findChild<QSpinBox *>("numNPOrients")->value();
    if(numOrientations > 1){
     for(int i = 1; i <= numOrientations; ++i){


         QString rotFileName = baseFileName;
          if( rotFileName.endsWith(".np") ||  rotFileName.endsWith(".pdb")){
               rotFileName.chop(chopChars);

          }
              rotFileName.append("_");
              rotFileName.append(QString::number(i) );
              if(isPDB){
                  rotFileName.append(".pdb");
              }
              else{
              rotFileName.append(".np");
}

         saveNP(rotFileName,true , isPDB);
     }
    }

}


void MainWindow::on_actionTips_triggered()
{
     tipsWindowI->exec();
}

void MainWindow::updateGraphicsWindow(){

    QGraphicsView* graphicsBox = this->findChild<QGraphicsView *>("graphicsView") ;
    int gbWidth = graphicsBox->width() ;
    int gbHeight = graphicsBox->height() ;
    int scaleBasis = gbWidth;
    if(gbHeight < gbWidth){
        scaleBasis = gbHeight;
    }

    std::vector< std::vector<double> > plotObjects;
    double innerBound = npInnerRadius;
    double delta=2;
    double outerBound = npOuterRadius;
    double maxRadius = 0;
        QPen outlinePen(Qt::black);


    //first pass: find the inner/outer limits and figure out plot order. we map z to up and down the screen, x to left and right, y to back-to-front.
    for( const auto& np: npBeads){
    //  qDebug() << beadTypes.size() << " found bead type ID: " << np.beadTypeID << "\n";
      //  if( (int)beadTypes.size() > np.beadTypeID){
           if(std::find(knownBeadTypes.begin(), knownBeadTypes.end(), np.beadTypeID) != knownBeadTypes.end() ){
            //recognised bead, add to plot
            auto npBeadType = beadTypes[ np.beadTypeID ] ;
            //qDebug() << "Found np of radius: " << npBeadType.radius << "\n";
            std::vector<double> plotCircle{ np.x, np.z, np.y + npBeadType.radius ,npBeadType.radius } ;

            plotObjects.emplace_back(  plotCircle );
            if(maxRadius < npBeadType.radius){
                maxRadius = npBeadType.radius;
            }

        }

    }



    double unbindingBound = outerBound + delta;


    double scaleFactor =0.5*scaleBasis/( unbindingBound + 1.0 );
    scene.clear( );


    //scene.addEllipse(xleft*scaleFactor, -zup*scaleFactor, 2*pc[3]*scaleFactor , 2*pc[3]*scaleFactor  ,outlinePen , whiteFill);

    //second pass: plot everything in order
    //QGraphicsScene scene;
    std::sort(plotObjects.begin(), plotObjects.end(),
              [](const std::vector<double>& a, const std::vector<double>& b) {
      return a[2] < b[2];
    });

        QBrush whiteFill(Qt::white);

    for(const auto& pc: plotObjects){
     //qDebug() << " plotting bead at max-y-value " << pc[2] << "\n";
     double xleft = pc[0] - pc[3];
     double zup = -pc[1]  - pc[3];
     scene.addEllipse(xleft*scaleFactor, zup*scaleFactor, 2*pc[3]*scaleFactor , 2*pc[3]*scaleFactor  ,outlinePen , whiteFill);
    }

    QPen dashPen(Qt::DashLine );
    QBrush redFill(Qt::red);
    outlinePen.setWidth(2);

        graphicsBox->setScene( &scene );
    //qDebug() << "updating graphics view \n" ;
    //qDebug() << "inner radius: " << innerBound << "\n";
    //qDebug() << "outer radius: " << unbindingBound << "\n";
     dashPen.setBrush(Qt::blue);
     innerBound = innerBound+0.1;
     scene.addEllipse(-innerBound*scaleFactor ,-innerBound*scaleFactor,2*innerBound*scaleFactor , 2*innerBound*scaleFactor  ,dashPen );
 dashPen.setBrush(Qt::red);
    scene.addEllipse(-unbindingBound*scaleFactor,-unbindingBound*scaleFactor,2*unbindingBound*scaleFactor , 2*unbindingBound*scaleFactor  ,dashPen );
    //"graphicsView"

    dashPen.setBrush(Qt::black);
    QBrush blackBrush(Qt::black);
    double originPointRadius = 0.2;
    scene.addEllipse(-originPointRadius*scaleFactor,-originPointRadius*scaleFactor,2.0*scaleFactor*originPointRadius, 2.0*scaleFactor*originPointRadius ,dashPen, blackBrush );
    double scaleBarLength = maxRadius;

    if(scaleBarLength < 5){
        scaleBarLength = 5;
    }
    else{
        scaleBarLength=   5 * round( scaleBarLength/5.0) ;
    }

    scene.addRect( -(1+unbindingBound)*scaleFactor,(1+unbindingBound)*scaleFactor , scaleBarLength*scaleFactor, 5, dashPen, blackBrush );
    QGraphicsTextItem *scaleLabel = scene.addText(QString::number(scaleBarLength) + " nm" ) ;
    scaleLabel->setPos( -(1+unbindingBound)*scaleFactor,(1+unbindingBound)*scaleFactor );
    graphicsBox->show();



}




void MainWindow::on_updateTables_clicked()
{
//reset bead types, re-initialise from defined beads
    QTableWidget *beadTable = this->findChild<QTableWidget *>("beadTypeTable");
     QTableWidget *npBeadTable = this->findChild<QTableWidget *>("beadTable");
    float minRadius = 0.1;

    /*
     for(int i = 0; i < (int)beadTypes.size(); ++i){

     beadTypes[i].surfaceDir = beadTable->item(i, 1)->data(0).toString().toStdString();
     beadTypes[i].hamakerFile = beadTable->item(i, 2)->data(0).toString().toStdString();

     beadTypes[i].radius = std::max(minRadius, beadTable->item(i, 3)->data(0).toFloat() );
     beadTypes[i].surfacePotential = beadTable->item(i, 4)->data(0).toFloat();
      beadTypes[i].surfFactor = beadTable->item(i, 5)->data(0).toFloat();
     beadTypes[i].coreFactor = beadTable->item(i, 6)->data(0).toFloat();
   beadTypes[i].ljCutoffVal = beadTable->item(i, 7)->data(0).toFloat();
   beadTypes[i].correctionOverride = beadTable->item(i, 8)->data(0).toInt();
     }


    //reset beads, re-initialise
    std::vector<int> popIndices;
     //QTableWidget *npBeadTable = this->findChild<QTableWidget *>("beadTable");
      for(int i = 0; i < (int)npBeads.size(); ++i){


      npBeads[i].beadTypeID = npBeadTable ->item(i, 0)->data(0).toInt();
      npBeads[i].x =npBeadTable ->item(i, 1)->data(0).toFloat();
      npBeads[i].y = npBeadTable ->item(i, 2)->data(0).toFloat();
      npBeads[i].z =npBeadTable ->item(i, 3)->data(0).toFloat();
      if(npBeadTable->item(i,4)->checkState() == Qt::Checked  ){
          popIndices.insert(popIndices.begin(), i );
      }
     }
      //delete beads which were marked as checked. Since we have the array in reverse order this removes from the end of the list first and keeps numbering intact
      for(const auto& delIndex: popIndices){
       npBeads.erase( npBeads.begin()+  delIndex);
       npBeadTable->removeRow(delIndex);


      }
     */
     beadTypes.clear();
     npBeads.clear();
    std::vector<int> popIndices;
    int numBeadTypeTableRows = beadTable->rowCount();
    int numNPBeadRows = npBeadTable->rowCount();

     for(int i = 0; i < numBeadTypeTableRows; ++i){

      int newBeadTypeID = beadTable->item(i, 0)->data(0).toInt();
     std::string newSurfaceDir = beadTable->item(i, 1)->data(0).toString().toStdString();
     std::string newHamakerFile = beadTable->item(i, 2)->data(0).toString().toStdString();

     float newRadius = std::max(minRadius, beadTable->item(i, 3)->data(0).toFloat() );
     float newSurfacePotential = beadTable->item(i, 4)->data(0).toFloat();
      float newSurfFactor = beadTable->item(i, 5)->data(0).toFloat();
     float newCoreFactor = beadTable->item(i, 6)->data(0).toFloat();
   float newLjCutoffVal = beadTable->item(i, 7)->data(0).toFloat();
   int newCorrectionOverride = beadTable->item(i, 8)->data(0).toInt();
   recieveNewNPBeadType(newHamakerFile, newSurfaceDir, newRadius, newSurfacePotential, newSurfFactor,newCoreFactor,newLjCutoffVal, newCorrectionOverride,newBeadTypeID,false);

     }

     for(int i = 0; i < numNPBeadRows; ++i){

     if(npBeadTable->item(i,4)->checkState() == Qt::Checked  ){
         popIndices.insert(popIndices.begin(), i );
     }
     else{
         recieveNewNPBead(  npBeadTable ->item(i, 1)->data(0).toFloat(), npBeadTable ->item(i, 2)->data(0).toFloat(), npBeadTable ->item(i, 3)->data(0).toFloat(), npBeadTable ->item(i, 0)->data(0).toInt(), false  );
     }
    }

     for(const auto& delIndex: popIndices){
      //npBeads.erase( npBeads.begin()+  delIndex);
      npBeadTable->removeRow(delIndex);


     }
     updateBindingRadii();
     //reset the view
     updateGraphicsWindow();
}



void MainWindow::on_beadTypeTable_cellChanged(int row, int column)
{
    //QTableWidget *beadTable = this->findChild<QTableWidget *>("beadTypeTable");
    //beadTable->row()
    //beadTypes[row].correctionOverride =  ;
   // qDebug() << " bead type table has changed \n";
}



void MainWindow::on_recenterNP_clicked()
{

    double xc=0;
    double yc=0;
    double zc=0;
    double totalMass = 0.001;

    //pass 1: get center of mass
    for( const auto& np: npBeads){
    //  qDebug() << beadTypes.size() << " found bead type ID: " << np.beadTypeID << "\n";
        if( std::find(knownBeadTypes.begin(), knownBeadTypes.end(), np.beadTypeID) != knownBeadTypes.end() ){
            //recognised bead, add to plot
            auto npBeadType = beadTypes[ np.beadTypeID ] ;
            //qDebug() << "Found np of radius: " << npBeadType.radius << "\n";
            double beadMass = pow( npBeadType.radius ,3) ;
             xc += np.x * beadMass;
             yc += np.y * beadMass;
             zc += np.z * beadMass;
             totalMass += beadMass;
        }

    }
    xc = xc/totalMass;
    yc = yc/totalMass;
    zc = zc/totalMass;

    //pass 2: shift each component

    for( int i = 0; i < (int)npBeads.size(); ++i){
    //  qDebug() << beadTypes.size() << " found bead type ID: " << np.beadTypeID << "\n";

             npBeads[i].x = npBeads[i].x - xc;
             npBeads[i].y = npBeads[i].y - yc;
             npBeads[i].z = npBeads[i].z -  zc;


    }

updateBeadTable();
}


void MainWindow::updateBindingRadii(){

    double innerBound  = 0;
    double outerBound = 1;
      if(this->findChild<QCheckBox *>("autoBounds")->isChecked() == true){

    for( const auto& np: npBeads){
    //  qDebug() << beadTypes.size() << " found bead type ID: " << np.beadTypeID << "\n";
        if(std::find(knownBeadTypes.begin(), knownBeadTypes.end(), np.beadTypeID) != knownBeadTypes.end()){
            //recognised bead, add to plot
            auto npBeadType = beadTypes[ np.beadTypeID ] ;
            //qDebug() << "Found np of radius: " << npBeadType.radius << "\n";
            if( innerBound < npBeadType.radius){
                innerBound = npBeadType.radius;
            }
            double npOuterDist = sqrt( pow(np.x,2) + pow(np.y,2) + pow(np.z,2) ) + npBeadType.radius ;
            if(outerBound < npOuterDist){
                outerBound = npOuterDist;
            }

        }

    }


        npInnerRadius = innerBound;
        npOuterRadius = outerBound;
        this->findChild<QLineEdit *>("lowerBoundLine")->setText(QString::number(innerBound));
         this->findChild<QLineEdit *>("upperBoundLine")->setText(QString::number(outerBound));
    }
    else{
         // qDebug() << "Manually updating bounding radii \n";
          npInnerRadius = this->findChild<QLineEdit *>("lowerBoundLine")->text().toFloat();
          npOuterRadius = this->findChild<QLineEdit *>("upperBoundLine")->text().toFloat();
          if(npOuterRadius < npInnerRadius){
              npOuterRadius = npInnerRadius;
          }
        //  qDebug() << " new inner: " << npInnerRadius << "\n";
        //  qDebug() << "new outer: "<< npOuterRadius << "\n";
      }
updateGraphicsWindow() ;

}

void MainWindow::on_autoBounds_stateChanged(int arg1)
{

    if(this->findChild<QCheckBox *>("autoBounds")->isChecked() == false){
        //automatic updating of the bounding radius is disabled - these fields can now be edited
        this->findChild<QLineEdit *>("lowerBoundLine")->setEnabled(true)  ;
        this->findChild<QLineEdit *>("upperBoundLine")->setEnabled(true)  ;
    }
    else{
        //automatic updating of the bounding radius is enabled - recalculate and disable these fields
        updateBindingRadii();
        this->findChild<QLineEdit *>("lowerBoundLine")->setEnabled(false)  ;
        this->findChild<QLineEdit *>("upperBoundLine")->setEnabled(false)  ;
    }
}


void MainWindow::on_lowerBoundLine_textEdited(const QString &arg1)
{
     if(this->findChild<QCheckBox *>("autoBounds")->isChecked() == false){
         updateBindingRadii();
     }

}


void MainWindow::on_upperBoundLine_textEdited(const QString &arg1)
{
    if(this->findChild<QCheckBox *>("autoBounds")->isChecked() == false){
        updateBindingRadii();
    }
}


void MainWindow::on_actionNew_triggered()
{


    auto confirm = QMessageBox::question(this, "NPDesigner", "Discard current NP?");
    if(confirm == QMessageBox::Yes){
        //delete all NP beads and then types
        beadTypes.clear();
        npBeads.clear();

        updateBeadTypeTable();
        updateBeadTable();
    }

}


void MainWindow::on_actionLoad_triggered()
{

    auto confirm = QMessageBox::question(this, "NPDesigner", "Discard current NP and load a new one?");
    if(confirm == QMessageBox::Yes){
        //delete all NP beads and then types
        beadTypes.clear();
        npBeads.clear();


            QString targetNPFile = QFileDialog::getOpenFileName(this, tr("NP File"),  this->uaGlobalPath,  tr("NP (*.np)"));

            QFile file(targetNPFile);
            if(file.open(QIODevice::ReadOnly | QIODevice::Text)) {
            QTextStream fileIn(&file);

            while(!fileIn.atEnd()){
                std::string lineIn = fileIn.readLine().toStdString();


                if(lineIn.size() > 3 && lineIn.substr(0, 1) != "#") {
                           try {



                              std::vector<std::string> results;
                               boost::split(results, lineIn, [](char c){return c == ',';});

                              if( results[0] == "TYPE"){
                              //register a new type
                               double newradiusval = std::stod(results[1]);
                               double newzetaval = std::stod(results[2])*1000;

                              double newcoreFactor = std::stod(results[3]);
                              double newsurfFactor = std::stod(results[4]);
                              //shape.emplace_back(std::stoi(results[5]);
                              std::string newhamakerFile = results[6];
                              std::string newpmfFile = (results[7]);
                              double newpmfCutoff = (std::stod(results[8]));
                              double newcorrectionType= (std::stoi(results[9]));


                              BeadType newNPBead(beadTypes.size()  , newhamakerFile , newpmfFile,newradiusval, newzetaval, newsurfFactor,newcoreFactor, newpmfCutoff, newcorrectionType   );
                              beadTypes.emplace_back(newNPBead);



                              }
                              else if( results[0] == "INNERBOUND"){
                              double fileInnerBound = std::stod(results[1]);
                              }
                              else if( results[0] == "OUTERBOUND"){
                              double fileOuterBound = std::stod(results[1]);
                              }
                              else if( results[0] == "BEAD"){
                              //place an NP bead
                              int    npBeadTypeVal = std::stoi(results[1]);
                              double xval = std::stod(results[2]);
                              double yval = std::stod(results[3]);
                              double zval = std::stod(results[4]);
                              NPBead newBead(npBeadTypeVal,xval,yval,zval);
                              npBeads.emplace_back( newBead  ) ;
                              }
                              else{
                              qDebug() << "NP input line " <<  QString::fromStdString(lineIn) << " \n not recognised, will attempt to parse using legacy system. WARNING: Mixing this with new-style will lead to misassigned beads. \n";

                              double xval = std::stod(results[0]);
                              double yval = std::stod(results[1]);
                              double zval = std::stod(results[2]);

                              double newradiusval = std::stod(results[3]);
                               double newzetaval = std::stod(results[4])*1000;
                               double newcoreFactor = std::stod(results[5]);
                               double newsurfFactor= std::stod(results[6]);
                               int shapeVal = std::stoi(results[7]) ;
                               std::string newhamakerFile  = results[8] ;
                               std::string newpmfFile  = results[9];
                               double newpmfCutoff = std::stod(results[10]) ;
                               int newcorrectionType = std::stoi(results[11]) ;
                               int foundBeadType = 0;
                               int npBeadTypeVal = 0;
                               int j = 0;
                               //test to see if the bead matches any beads which have already been seen - replace this with object based testing
                               for(j=0 ; j < (int)beadTypes.size(); ++j){
                                BeadType testBead = beadTypes[j];
                                if( abs(testBead.radius -newradiusval ) < 0.1 &&  abs( testBead.surfacePotential  -newzetaval) < 0.1
                                        && abs( testBead.coreFactor - newcoreFactor ) < 0.01 && abs( testBead.surfFactor - newsurfFactor) < 0.01
                                         && testBead.hamakerFile == newhamakerFile && testBead.surfaceDir == newpmfFile
                                        &&  abs( testBead.ljCutoffVal - newpmfCutoff ) < 0.01  && abs( testBead.correctionOverride - newcorrectionType) < 0.1 ){
                                foundBeadType = 1;
                                npBeadTypeVal = j;
                                }

                               }

                               if(foundBeadType == 0){
                               //register a new bead type
                                   npBeadTypeVal = beadTypes.size();
                                   BeadType newNPBead(npBeadTypeVal  , newhamakerFile , newpmfFile,newradiusval, newzetaval, newsurfFactor,newcoreFactor, newpmfCutoff, newcorrectionType   );
                                   beadTypes.emplace_back(newNPBead);

                       //   std::cout << "Legacy reader: Defining bead type " << npBeadTypeVal << " radius " << radius[npBeadTypeVal] << "\n";

                              }
                              //then add the NP bead
                              //std::cout<<"Legacy reader: Adding bead of type " << npBeadTypeVal << "\n";
                               NPBead newBead(npBeadTypeVal,xval,yval,zval);
                               npBeads.emplace_back( newBead  ) ;
                              }






                           //previous expression: sqrt( xval*xval + yval*yval + zval*zval   ) + radiusval;






                           }
                           catch (const std::invalid_argument& ia) {
                               qDebug() << "Error: Failed to parse NP file '" << targetNPFile << "'\n";
                               qDebug() << "Error: Failed on line: " << QString::fromStdString(lineIn) << "\n";
                               std::exit(1);
                           }

                       }













            }

            }

        updateBeadTypeTable();
        updateBeadTable();
    }




}




void MainWindow::on_loadMaterialSet_clicked()
{
    QString targetMaterialFile = QFileDialog::getOpenFileName(this, tr("Material File"),  this->uaGlobalPath,  tr("CSV (*.csv)"));
    //Load in the file, iterate over contents, update the material set

    QFile file(targetMaterialFile);
    if(file.open(QIODevice::ReadOnly | QIODevice::Text)) {
    QTextStream fileIn(&file);
    int numFound = 0;
    QComboBox *addBeadTypeMaterialBox = addBeadType->findChild<QComboBox *>("materialTypeBox") ;
    QComboBox *addShellMaterialBox = addShell->findChild<QComboBox *>("materialTypeBox") ;

    while(!fileIn.atEnd()){
        std::string lineIn = fileIn.readLine().toStdString();


        if(lineIn.substr(0, 1) != "#") {
            std::vector<std::string> results;
             boost::split(results, lineIn, [](char c){return c == ',';});
             //silicaquartz,surface/SiO2-Quartz,hamaker/SiO2_Quartz.dat,1
             //qDebug() << QString::fromStdString(results[0]) << " " << QString::fromStdString(results[1]) << " " <<  QString::fromStdString(results[2] )<< "\n";
             if(results.size()==4){
             materialTypes.emplace_back( MaterialType( QString::fromStdString(results[0])   , QString::fromStdString(results[1]) ,QString::fromStdString(results[2])     ) ) ;
             addBeadTypeMaterialBox->addItem( QString::fromStdString(results[0]) ,   QList<QVariant>() <<  QString::fromStdString(results[1]) <<  QString::fromStdString(results[2])) ;
             addShellMaterialBox->addItem( QString::fromStdString(results[0]) ,   QList<QVariant>() <<  QString::fromStdString(results[1]) <<  QString::fromStdString(results[2])) ;

             numFound += 1;
             }

        }
    }

    QMessageBox msgBox;
    msgBox.setText("Loaded " + QString::number(numFound));
       msgBox.exec();

}
}

void MainWindow::on_uaPath_textEdited(const QString &uaDir)
{
    this->uaGlobalPath = uaDir;
}


void MainWindow::on_beadTable_customContextMenuRequested(const QPoint &pos)
{
    npBeadTableMenu.popup( this->findChild<QTableWidget *>("beadTable")->viewport()->mapToGlobal(pos)    ) ;
}

void MainWindow::removeMultipleNPBeads(){
    qDebug() << "removing multiple beads";
    QTableWidget *npBeadTable =  this->findChild<QTableWidget *>("beadTable");
    for(auto selectedItem: npBeadTable->selectedItems() ){
     int targetRow = selectedItem->row();
     if(targetRow > -1){
        npBeadTable->item( targetRow, 4)-> setCheckState(Qt::Checked);
     }
    }

    on_updateTables_clicked();
}

