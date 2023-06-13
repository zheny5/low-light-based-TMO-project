// qt_test.h
#pragma once

#include <QtWidgets/QMainWindow>
#include <Qlabel>
#include <QPushButton>
#include <Qopenglwidget>
#include "ui_tonemapper1.h"
#include "GlWidget.h"

class qt_test_gl : public QMainWindow
{
    Q_OBJECT

public:
    qt_test_gl(QWidget* parent = Q_NULLPTR);
    void onButtonClick();

private:
    Ui::Tonemapper1Class ui;
    QLabel* lab;
    QPushButton* button;

    GlWidget* opengl_widget;
};
