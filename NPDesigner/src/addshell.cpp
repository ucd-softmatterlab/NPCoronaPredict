#include "addshell.h"
#include "ui_addshell.h"
#include <QFileDialog>


AddShell::AddShell(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::AddShell)
{
    ui->setupUi(this);
}

AddShell::~AddShell()
{
    delete ui;
}

void AddShell::on_buttonBox_accepted()
{
    std::string surfDir = this->findChild<QPlainTextEdit *>("surfDirIn")->toPlainText().toStdString()   ;
    std::string hamakerFile = this->findChild<QPlainTextEdit *>("hamakerIn")->toPlainText().toStdString()   ;
    double outerRadius = (this->findChild<QLineEdit *>("outerRadius")->text()).toFloat()   ;
     double innerRadius = (this->findChild<QLineEdit *>("innerRadius")->text()).toFloat()   ;
     double zeta = (this->findChild<QLineEdit *>("zeta")->text()).toFloat()   ;
     double ljCutoffVal = (this->findChild<QLineEdit *>("ljCutoff")->text()).toFloat()   ;

       emit( sendNewShell(hamakerFile, surfDir, innerRadius, outerRadius, zeta, ljCutoffVal ));
//emit( sendNewNPBeadType( hamakerFile , surfDir , outerRadius, zeta, 1.0, 1.0, ljCutoffVal, 0   )) ;
//emit( sendNewNPBeadType( hamakerFile , surfDir , innerRadius, 0.0, -1.0, -1.0, ljCutoffVal, 0  )) ;

}


void AddShell::on_findHamaker_clicked()
{
    //to do - add relative pathing in, use QDir dir("BasePath") + dir.relativeFilePath
    //  QString UABaseDir =  qobject_cast <MainWindow *>(this->parent()  )->uaGlobalPath;

      QString hamFile = QFileDialog::getOpenFileName(this, tr("Hamaker File"),  this->searchPath,  tr("Hamaker (*.dat)"));
    if(hamFile != ""){
        this->findChild<QPlainTextEdit *>("hamakerIn")->setPlainText(hamFile);
    }
}


void AddShell::on_findSurf_clicked()
{

    QString surfDir = QFileDialog::getExistingDirectory(this, tr("Surface Directory"),  this->searchPath, QFileDialog::ShowDirsOnly);
    if(surfDir != ""){
        this->findChild<QPlainTextEdit *>("surfDirIn")->setPlainText(surfDir);
    }
}

