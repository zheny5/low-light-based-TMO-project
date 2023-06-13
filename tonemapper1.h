#pragma once

#include <QtWidgets/QWidget>
#include "ui_tonemapper1.h"

class Tonemapper1 : public QWidget
{
    Q_OBJECT

public:
    Tonemapper1(QWidget *parent = nullptr);
    ~Tonemapper1();

private:
    Ui::Tonemapper1Class ui;
};
