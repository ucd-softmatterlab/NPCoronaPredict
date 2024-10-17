#ifndef ADDBRUSH_H
#define ADDBRUSH_H

#include <QDialog>

namespace Ui {
class AddBrush;
}

class AddBrush : public QDialog
{
    Q_OBJECT

public:
    explicit AddBrush(QWidget *parent = nullptr);
    ~AddBrush();
    void updateDensityBox();
    void setRadialToSuggest();
    std::vector<double> beadRadii;
    double currentOuterRadius;
    int lastBeadIndex;
signals:
    void  sendNewBrush(double brushOccupancy, double brushRadialDist, int beadID, bool forceAttach);
private slots:
    void on_buttonBox_accepted();

    void on_brushOccupancy_valueChanged(double arg1);

    void on_brushBeadID_currentIndexChanged(int index);

    void on_brushRadialDist_editingFinished();

private:
    Ui::AddBrush *ui;
};

#endif // ADDBRUSH_H
