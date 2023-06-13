#include "stdafx.h"
#include "tonemapper1.h"
#include "ImageViewer.h"
#include <QtWidgets/QApplication>
#include"qt_test_gl.h"
int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    ImageViewer w;
    //qt_test_gl w;
    w.show();
    return a.exec();

}
