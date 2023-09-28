#include "tipswindow.h"
#include "ui_tipswindow.h"

tipsWindow::tipsWindow(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::tipsWindow)
{
    ui->setupUi(this);
}

tipsWindow::~tipsWindow()
{
    delete ui;
}
