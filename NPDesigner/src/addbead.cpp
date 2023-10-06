#include "addbead.h"
#include "ui_addbead.h"
#include <QtDebug>
#include <QMessageBox>
AddBead::AddBead(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::AddBead)
{
    ui->setupUi(this);
}

AddBead::~AddBead()
{
    delete ui;
}



void AddBead::on_buttonBox_accepted()
{
   QString xstring = this->findChild<QLineEdit *>("newBeadX")->text()   ;
   double xval = xstring.toFloat();
   QString ystring = this->findChild<QLineEdit *>("newBeadY")->text()   ;
   double yval = ystring.toFloat();
   QString zstring = this->findChild<QLineEdit *>("newBeadZ")->text()   ;
   double zval = zstring.toFloat();
   //qDebug () << "Looking for bead type ID \n";
   int beadTypeID = this->findChild<QComboBox *>("beadTypeID")->currentIndex() ;
  // qDebug() << "New bead: " << xval << " " << yval << " " << zval << " type: " << beadTypeID << "\n";

   if(beadTypeID == -1){
       qDebug() << "No valid bead selected \n";

       QMessageBox::warning(this, tr("NPDesigner"),       tr("No valid bead selected. Please register a bead first.\n" ) );

   }
   else{
   emit( sendNewNPBead(xval,yval,zval,beadTypeID, true)) ;
   }

}

