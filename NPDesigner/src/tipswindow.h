#ifndef TIPSWINDOW_H
#define TIPSWINDOW_H

#include <QDialog>

namespace Ui {
class tipsWindow;
}

class tipsWindow : public QDialog
{
    Q_OBJECT

public:
    explicit tipsWindow(QWidget *parent = nullptr);
    ~tipsWindow();

private:
    Ui::tipsWindow *ui;
};

#endif // TIPSWINDOW_H
