#ifndef ADDBEAD_H
#define ADDBEAD_H

#include <QDialog>

namespace Ui {
class AddBead;
}

class AddBead : public QDialog
{
    Q_OBJECT

public:
    explicit AddBead(QWidget *parent = nullptr);
    ~AddBead();
signals:
    void  sendNewNPBead(double x, double y, double z, int beadID, bool doUpdate);
private slots:
    void on_buttonBox_accepted();

private:
    Ui::AddBead *ui;
};

#endif // ADDBEAD_H
