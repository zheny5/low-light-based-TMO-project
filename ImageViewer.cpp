/*************************************************************
* File: ImageViewer.cpp
* Description: Implementation of the ImageViewer class, an image viewer application.
**************************************************************/

#include "stdafx.h"
#include "ImageViewer.h"
#include <QFileDialog>
#include <QScrollBar>
#include <QImageWriter>
#include <iostream>
#include <QDebug>
#include "Percep_low_light.h"
/**
 * Constructor for ImageViewer class.
 * Initializes the UI and connects signals to slots.
 * @param parent The parent widget.
 */
ImageViewer::ImageViewer(QWidget* parent)
    : QMainWindow(parent)
{
    ui.setupUi(this);

    // Map buttons to labels
    buttonLabelMap.insert(ui.pushButton, ui.label);
    buttonLabelMap.insert(ui.pushButton_2, ui.label);
    buttonLabelMap.insert(ui.pushButton_3, ui.label);
    buttonLabelMap.insert(ui.pushButton_4, ui.label);
    buttonLabelMap.insert(ui.pushButton_5, ui.label_2);
    buttonLabelMap.insert(ui.pushButton_6, ui.label_2);
    buttonLabelMap.insert(ui.pushButton_7, ui.label_2);
    buttonLabelMap.insert(ui.pushButton_8, ui.label_2);

    // Initialize scale factors for labels
    labelScaleFactorMap.insert(ui.label, 1.0);
    labelScaleFactorMap.insert(ui.label_2, 1.0);

    // Connect signals to slots
    connect(ui.actionOpen, &QAction::triggered, this, &ImageViewer::loadFile); // Custom function
    connect(ui.actionExit, &QAction::triggered, this, &QWidget::close); // System function
    connect(ui.actionSaveAs, &QAction::triggered, this, &ImageViewer::saveAs);
    connect(ui.pushButton, &QPushButton::clicked, this, &ImageViewer::zoomOut);
    connect(ui.pushButton_2, &QPushButton::clicked, this, &ImageViewer::zoomIn);
    connect(ui.pushButton_3, &QPushButton::clicked, this, &ImageViewer::normalSize);
    connect(ui.pushButton_4, &QPushButton::clicked, this, &ImageViewer::fitToWindow);
    connect(ui.pushButton_5, &QPushButton::clicked, this, &ImageViewer::fitToWindow);
    connect(ui.pushButton_6, &QPushButton::clicked, this, &ImageViewer::normalSize);
    connect(ui.pushButton_7, &QPushButton::clicked, this, &ImageViewer::zoomOut);
    connect(ui.pushButton_8, &QPushButton::clicked, this, &ImageViewer::zoomIn);
}

/**
 * Loads an image file and displays it in the UI.
 * Supports PNG, BMP, JPG, and HDR (EXR) image formats.
 * If the loaded file is not an image, it is processed using tonemapping.
 */
void ImageViewer::loadFile()
{
    QFileDialog dialog(this, "Open Image", NULL, "Image Files (*.png) \n Image Files (*.bmp *.jpg) \n HDR Image Files (*.exr)"); // Priority on bmp format
    if (dialog.exec() != QDialog::Accepted)
    {
        return;
    }
    QString fullPath = dialog.selectedFiles().first(); // Get the full path of the image
    if (fullPath.endsWith("png", Qt::CaseInsensitive) || fullPath.endsWith("jpg", Qt::CaseInsensitive) || fullPath.endsWith("jpeg", Qt::CaseInsensitive))
    {
        // Load and display the image
        image.load(fullPath);
        QFileInfo fi(fullPath); // Get the file name
        const QString title = QString("%1 - %2").arg(fi.fileName()).arg(windowTitle());
        setWindowTitle(title);
        const QString msg = QString("Resolution: %1×%2, Bit Depth: %3").arg(image.width()).arg(image.height()).arg(image.depth());
        ui.statusBar->showMessage(msg);
        ui.label->setPixmap(QPixmap::fromImage(image));
        ui.label->adjustSize();
        ui.scrollArea->setWidget(ui.label);
        ui.scrollArea_6->setWidget(ui.label_2);
        qtdisplayImgTest();
    }
    else 
    {
        // Process the file using tonemapping
        displayTonemapped(fullPath, "percep_low_light");
        qtdisplayImgResult();
    }
    return;
}

/**
 * Saves the currently displayed image to a file.
 * Supports BMP, PNG, and JPG image formats.
 */
void ImageViewer::saveAs()
{
    QFileDialog dialog(this, "Save Image", NULL, "Image Files (*.bmp) \n Image Files (*.png *.jpg)");
    dialog.setAcceptMode(QFileDialog::AcceptSave);
    if (dialog.exec() != QDialog::Accepted)
    {
        return;
    }
    QImageWriter writer(dialog.selectedFiles().first());
    writer.write(image);
}

/**
 * Zooms out the displayed image by reducing the scale factor.
 */
void ImageViewer::zoomOut() {
    QLabel* label = getLabel();
    if (label) {
        labelScaleFactorMap[label] *= 0.8;
        label->setScaledContents(true);
        label->resize(labelScaleFactorMap[label] * label->pixmap(Qt::ReturnByValue).size());
    }
}

/**
 * Zooms in the displayed image by increasing the scale factor.
 */
void ImageViewer::zoomIn() {
    QLabel* label = getLabel();
    if (label) {
        labelScaleFactorMap[label] *= 1.2;
        label->setScaledContents(true);
        label->resize(labelScaleFactorMap[label] * label->pixmap(Qt::ReturnByValue).size());
    }
}

/**
 * Resets the displayed image to its original size.
 */
void ImageViewer::normalSize() {
    QLabel* label = getLabel();
    if (label) {
        labelScaleFactorMap[label] = 1.0;
        label->setScaledContents(true);
        label->resize(labelScaleFactorMap[label] * label->pixmap(Qt::ReturnByValue).size());
    }
}

/**
 * Adjusts the displayed image to fit the window size.
 */
void ImageViewer::fitToWindow() {
    QLabel* label = getLabel();
    if (label) {
        label->setScaledContents(true);
        label->resize(ui.scrollArea->width() - 1, ui.scrollArea->height() - 1);
    }
}

/**
 * Returns the label associated with the clicked button.
 * @return The label widget or nullptr if not found.
 */
QLabel* ImageViewer::getLabel() {
    QPushButton* button = qobject_cast<QPushButton*>(sender());
    if (button) {
        return buttonLabelMap.value(button);
    }
    else {
        return nullptr;
    }
}

/**
 * Displays a tone-mapped version of the input image using the specified method.
 * Currently, only the "percep_low_light" method is implemented.
 * @param fullPath The full path of the input image file.
 * @param method The tonemapping method to use.
 */
void ImageViewer::displayTonemapped(QString fullPath, std::string method)
{
    if (method == "percep_low_light")
    {
        Percep_low_light pllight;
        auto exrimage = pllight.converto_Mat(fullPath);
        pllight.HDRToLMSR(pllight.tonemapMats_.hdrImage, pllight.tonemapMats_.LMtx, pllight.tonemapMats_.MMtx, pllight.tonemapMats_.SMtx, pllight.tonemapMats_.RMtx);
        pllight.purkinje(pllight.tonemapMats_.LMtx, pllight.tonemapMats_.MMtx, pllight.tonemapMats_.SMtx, pllight.tonemapMats_.RMtx, pllight.tonemapMats_.lms, pllight.tonemapMats_.blend, &pllight.tonemapVals_);
        pllight.lms2display(pllight.tonemapMats_.lms, pllight.tonemapMats_.lms2displayOutput);
        pllight.tonemapping4lms(pllight.tonemapMats_.lms2displayOutput);
        cv::Mat reduceRangeInput = cv::imread("AfterDurandTonemapping.ppm", 1);
        //cv::Mat reduceRangeInput = cv::imread("D:/Research/ToneMap/Code/Perceptually based tone mapping for low-light conditions/a42-kirk/LOWLIGHT_LINUX/LOWLIGHT_LINUX/src/images/AfterDurandTonemapping.ppm", 1);
        pllight.reduceRange(&reduceRangeInput, pllight.tonemapMats_.blend, reduceRangeInput.rows, reduceRangeInput.cols, "./output.png", *pllight.tonemapMats_.inputImage);
    }
}

/**
 * Displays the input OpenCV image for debugging purposes.
 * @param image The OpenCV image to display.
 */
void ImageViewer::debugOpenCVimage(cv::Mat image)
{
    int type = image.type();
    if (type == CV_32FC3) cv::normalize(image, image, 0, 255, cv::NORM_MINMAX, CV_8UC3);
    cv::namedWindow("Display Image", cv::WINDOW_NORMAL);
    cv::imshow("Display Image", image);
    cv::waitKey(0);
}

/**
 * Displays a test image in the label_2 widget using Qt.
 * Used for testing purposes.
 */
void ImageViewer::qtdisplayImgTest()
{
    cv::Mat image = cv::imread("D:/Research/ToneMap/sihdr/input/clip_95/092.png");
    tonemappedImage = mat2qim(image);
    ui.label_2->setPixmap(QPixmap::fromImage(tonemappedImage));
    ui.label_2->adjustSize();
}
void ImageViewer::qtdisplayImgResult()
{
    cv::Mat image1 = cv::imread("./AfterDurandTonemapping.ppm");
    //cv::Mat image1 = cv::imread("./BeforeDurandTonemapping.exr", cv::IMREAD_ANYCOLOR | cv::IMREAD_ANYDEPTH);
    QImage image1_beforetmo = mat2qim(image1);
    ui.label->setPixmap(QPixmap::fromImage(image1_beforetmo));
    ui.label->adjustSize();
    ui.scrollArea->setWidget(ui.label);
    ui.scrollArea_6->setWidget(ui.label_2);
    cv::Mat image2 = cv::imread("./output.png");
    tonemappedImage = mat2qim(image2);
    ui.label_2->setPixmap(QPixmap::fromImage(tonemappedImage));
    ui.label_2->adjustSize();
}
/**
 * Converts an OpenCV Mat to a QImage.
 * @param mat The OpenCV Mat to convert.
 * @return The converted QImage.
 */
QImage ImageViewer::mat2qim(cv::Mat& mat)
{
    cvtColor(mat, mat, cv::COLOR_BGR2RGB);
    QImage qim((const unsigned char*)mat.data, mat.cols, mat.rows, mat.step,
        QImage::Format_RGB888);
    return qim;
}
