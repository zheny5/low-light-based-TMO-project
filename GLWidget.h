#pragma once

#ifndef __GL_WIDGET_H__
#define __GL_WIDGET_H__

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QTimer>
#include <opencv2/opencv.hpp>
//#include <gl/GL.h>
//#include<gl/freeglut.h>
//#include<gl/GLU.h>
class GlWidget : public QOpenGLWidget, protected QOpenGLFunctions
{
    Q_OBJECT

public:
    GlWidget(QWidget* parent = 0);
    ~GlWidget();

protected:
    void initializeGL() override;
    void resizeGL(int w, int h) override;
    void paintGL() override;

private:
    QTimer* timer;
};

#endif // WIDGET_H
