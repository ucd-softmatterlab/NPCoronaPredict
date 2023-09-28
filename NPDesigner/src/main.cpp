#include "mainwindow.h"
#include <QApplication>





int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;





    //register the buttons on the main window
    //auto newBeadTypeBtn = w.findChild<QPushButton*>("newBeadTypeButton");
   // w.connect( newBeadTypeBtn, SIGNAL (clicked()), QApplication::instance(), SLOT (openAddBead()));



    w.show();
    return a.exec();
}
