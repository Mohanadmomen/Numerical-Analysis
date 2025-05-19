#include "mainwindow.h"
#include <QIcon>
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;

    w.setWindowTitle("Numerical Analysis Tool");

    // Set the window icon
    w.setWindowIcon(QIcon(":Numerical/Icon.png"));

    w.show();
    return a.exec();
}
