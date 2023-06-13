#include"stdafx.h"
#include "GlWidget.h"

GlWidget::GlWidget(QWidget * parent)
	: QOpenGLWidget(parent)
{
}

GlWidget::~GlWidget()
{

}

void GlWidget::initializeGL()
{
	this->initializeOpenGLFunctions(); // Initialize OpenGL functions for the current context
	timer = new QTimer();
	timer->setInterval(33);
	connect(timer, SIGNAL(timeout()), this, SLOT(update()));
	timer->start();
}

void GlWidget::resizeGL(int w, int h)
{
	glViewport(0, 0, w, h); // Adjust the viewport
}

void GlWidget::paintGL()
{
	static float r = 0.0, g = 0.0, b = 0.0;
	r += 0.01;
	g += 0.02;
	b += 0.03;
	if (r > 1) r = 0;
	if (g > 1) g = 0;
	if (b > 1) b = 0;
	glClearColor(r, g, b, 1.0f); // Clear the screen
	glClear(GL_COLOR_BUFFER_BIT); // Clear the color buffer

	//cv::Mat image = cv::imread("D:/Research/ToneMap/sihdr/input/clip_95/002.png");



}