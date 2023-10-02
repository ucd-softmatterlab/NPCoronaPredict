#ifndef ADDBEADTYPE_H
#define ADDBEADTYPE_H

#include <QDialog>

namespace Ui {
class AddBeadType;
}

class AddBeadType : public QDialog
{
    Q_OBJECT

public:
    explicit AddBeadType(QWidget *parent = nullptr);
    ~AddBeadType();
    QString searchPath = "";
signals:
    void  sendNewNPBeadType(std::string hamakerFileIn, std::string surfaceDirIn, float   radiusIn, float   surfacePotentialIn, float surfFactorIn, float coreFactorIn, float ljCutoffIn, int correctionOverrideIn);
private slots:
    void on_buttonBox_accepted();

    void on_surfFindDirButton_clicked();

    void on_hamFindButton_clicked();

    void on_materialTypeBox_currentIndexChanged(int index);

private:
    Ui::AddBeadType *ui;
};

#endif // ADDBEADTYPE_H
