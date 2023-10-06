//Define activities for the AddBeadType window interface for adding bead types

#include "addbeadtype.h"
#include "ui_addbeadtype.h"
#include <QtDebug>
#include <QFileDialog>
#include <QComboBox>
AddBeadType::AddBeadType(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::AddBeadType)
{
    ui->setupUi(this);
}

AddBeadType::~AddBeadType()
{
    delete ui;
}

void AddBeadType::on_buttonBox_accepted()
{
   QString surfDirString = this->findChild<QPlainTextEdit *>("surfDir")->toPlainText()   ;
   QString hamFileString = this->findChild<QPlainTextEdit *>("hamakerFile")->toPlainText()   ;

   double surfFactor = (this->findChild<QLineEdit *>("surfFactor")->text()).toFloat()   ;
   double coreFactor = (this->findChild<QLineEdit *>("coreFactor")->text()).toFloat()   ;
   double zetaVal =  (this->findChild<QLineEdit *>("zeta")->text()).toFloat()   ;
  double ljCutoffVal =  (this->findChild<QLineEdit *>("ljCutoff")->text()).toFloat()   ;
  double radius  =  (this->findChild<QLineEdit *>("radius")->text()).toFloat()   ;
  int correctionOverride = 1;

   //qDebug() << "New bead: " << xval << " " << yval << " " << zval << " type: " << beadTypeID << "\n";
   emit( sendNewNPBeadType( hamFileString.toStdString() , surfDirString.toStdString() , radius, zetaVal, surfFactor, coreFactor, ljCutoffVal, correctionOverride   )) ;


}


void AddBeadType::on_surfFindDirButton_clicked()
{


    QString surfDir = QFileDialog::getExistingDirectory(this, tr("Surface Directory"), this->searchPath, QFileDialog::ShowDirsOnly);
    if(surfDir != ""){
         QDir uaDir(this->searchPath);
        QString localSurfPath = uaDir.relativeFilePath(surfDir);

        this->findChild<QPlainTextEdit *>("surfDir")->setPlainText(localSurfPath);
    }
}


void AddBeadType::on_hamFindButton_clicked()
{
    QString hamFile = QFileDialog::getOpenFileName(this, tr("Hamaker File"),  this->searchPath,  tr("Hamaker (*.dat)"));
    if(hamFile != ""){
        QDir uaDir(this->searchPath);
       QString localHamPath = uaDir.relativeFilePath(hamFile);

        this->findChild<QPlainTextEdit *>("hamakerFile")->setPlainText(localHamPath);
    }
}


void AddBeadType::on_materialTypeBox_currentIndexChanged(int index)
{
    //get new material type chosen
   // QVariant currentData = this->findChild<QComboBox *>("materialTypeBox")->currentData();
    auto currentData =this->findChild<QComboBox *>("materialTypeBox")->currentData().value<QList<QVariant>>();
   // qDebug() << currentData[0].toString() << " " << currentData[1].toString() << "\n";
     this->findChild<QPlainTextEdit *>("surfDir")->setPlainText(currentData[0].toString());
     this->findChild<QPlainTextEdit *>("hamakerFile")->setPlainText(currentData[1].toString());

}

