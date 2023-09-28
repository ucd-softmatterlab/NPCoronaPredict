#ifndef ADDSHELL_H
#define ADDSHELL_H

#include <QDialog>


namespace Ui {

class AddShell;
}

class AddShell : public QDialog
{
    Q_OBJECT

public:
    explicit AddShell(QWidget *parent = nullptr);
    ~AddShell();
    QString searchPath = "";
signals:
    void  sendNewShell(std::string hamakerFileIn, std::string surfaceDirIn, float   innerRadius, float outerRadius, float   surfacePotentialIn,  double ljCutoff);
private slots:
    void on_buttonBox_accepted();

    void on_findHamaker_clicked();

    void on_findSurf_clicked();

private:
    Ui::AddShell *ui;
};

#endif // ADDSHELL_H
