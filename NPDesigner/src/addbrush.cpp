#include "addbrush.h"
#include "ui_addbrush.h"
#include <QDebug>
#include <QMessageBox>
AddBrush::AddBrush(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::AddBrush)
{
    ui->setupUi(this);
}

AddBrush::~AddBrush()
{
    delete ui;
}

void AddBrush::on_buttonBox_accepted()
{

    //qDebug() << " starting brush \n";
    bool forceAttach = false;

    if(   this->findChild<QCheckBox *>("forceAttachBox")->isChecked() == false  ){
        forceAttach = true;
    }

    double brushOccupancy = this->findChild<QDoubleSpinBox *>("brushOccupancy")->value() ;
    // qDebug() << " occupancy: " << brushOccupancy << "\n";
    QString ystring = this->findChild<QLineEdit *>("brushRadialDist")->text()   ;
    double brushRadialDist = ystring.toFloat();
     //    qDebug() << " dist: " << brushRadialDist<< "\n";
 //   int beadTypeID = (this->findChild<QLineEdit *>("brushBeadID")->text() ).toInt() ;
     int beadTypeID =  this->findChild<QComboBox  *>("brushBeadID")->currentIndex( );
//qDebug() << " bead type: " << beadTypeID<< "\n";
   // qDebug() << " sending brush \n";
   lastBeadIndex = beadTypeID;
     if(beadTypeID == -1){
         qDebug() << "No valid brush bead selected \n";

         QMessageBox::warning(this, tr("NPDesigner"),       tr("No valid brush bead selected. Please register a bead first.\n" ) );

     }
     else{
    emit( sendNewBrush(brushOccupancy,brushRadialDist,beadTypeID,forceAttach)) ;
     }


}


void AddBrush::on_brushOccupancy_valueChanged(double arg1)
{
    //qDebug() << arg1 << "\n";
  updateDensityBox();
}

void AddBrush::updateDensityBox(){
double brushOccupancy = this->findChild<QDoubleSpinBox *>("brushOccupancy")->value() ;
int beadTypeID =  this->findChild<QComboBox  *>("brushBeadID")->currentIndex( );
//QString ystring = this->findChild<QLineEdit *>("brushRadialDist")->text()   ;
double brushRadialDist = this->findChild<QLineEdit *>("brushRadialDist")->text().toFloat();


double beadRadius = 5;//this->beadRadii[beadTypeID];
if( (int)this->beadRadii.size() > beadTypeID){
    beadRadius = beadRadii[beadTypeID] ;
}

double estNumBeads = brushOccupancy * (4 * brushRadialDist*brushRadialDist)/(1.1*  beadRadius*beadRadius) ;
double estDensity = estNumBeads / (4 * brushRadialDist*brushRadialDist * 3.1415 ) ;
this->findChild<QLineEdit *>("densityOutBox")->setText( QString::number(estDensity)    ) ;
this->findChild<QLineEdit *>("numbeadsOutBox")->setText( QString::number(estNumBeads)    ) ;
double suggestRadius = beadRadius+currentOuterRadius;
this->findChild<QLineEdit *>("outerLayerTarget")->setText( QString::number(suggestRadius)    ) ;



}


void AddBrush::setRadialToSuggest(){
    double suggestRadius = this->findChild<QLineEdit *>("outerLayerTarget")->text().toFloat();
    this->findChild<QLineEdit *>("brushRadialDist")->setText( QString::number(suggestRadius)    ) ;
}

void AddBrush::on_brushBeadID_currentIndexChanged(int index)
{
    updateDensityBox();
    setRadialToSuggest();
}


void AddBrush::on_brushRadialDist_editingFinished()
{
    updateDensityBox();
    setRadialToSuggest();
}


