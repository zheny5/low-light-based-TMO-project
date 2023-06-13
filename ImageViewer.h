#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_ImageViewer.h"
#include <opencv2/opencv.hpp>
#include <QMap>
#include <opencv2/photo/photo.hpp>

//#include<QScrollArea>

class ImageViewer : public QMainWindow
{
    Q_OBJECT

public:
    ImageViewer(QWidget* parent = Q_NULLPTR);
    //private slots:
    void loadFile();
    void saveAs(); // Save the image
    void zoomOut(); // Zoom out
    void zoomIn(); // Zoom in
    void normalSize(); // Restore to normal size
    void fitToWindow(); // Fit to window
    void ImageViewer::displayTonemapped(QString fullPath, std::string method);
    void debugOpenCVimage(cv::Mat image);
private:
    Ui::ImageViewerClass ui;
    QImage image;
    QImage tonemappedImage;
    double scaleFactor; // Scale factor
    QImage mat2qim(cv::Mat& mat);//cv::mat to qimage
    QMap<QPushButton*, QLabel*> buttonLabelMap;//for mapping the label to the button
    QMap<QLabel*, qreal> labelScaleFactorMap;//for mapping the scale to the label
    QLabel* getLabel();
    void qtdisplayImgTest();
    void qtdisplayImgResult();
public:
};
