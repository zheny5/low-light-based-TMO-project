#include "stdafx.h"
#include "Percep_low_light.h"
#include <opencv2/opencv.hpp>
#include <iostream>
#include<opencv2/xphoto/tonemap.hpp>
#include<algorithm>
namespace
{
	void G(cv::Mat* signal, float k1, float k2, cv::Mat* val) //helper function for the purkinje function
	{
		int rows, cols;
		rows = signal->rows;
		cols = signal->cols;

		for (int r = 0; r < rows; r++)
		{
			for (int c = 0; c < cols; c++)
			{
				val->at<float>(r, c) = 1 / (pow((1 + k1 * signal->at<float>(r, c)), k2));
			}
		}
	}
}
Percep_low_light::Percep_low_light()
{
	// Initialize the Mats
	tonemapMats_.LMtx = new cv::Mat();
	tonemapMats_.MMtx = new cv::Mat();
	tonemapMats_.SMtx = new cv::Mat();
	tonemapMats_.RMtx = new cv::Mat();
	tonemapMats_.hdrImage = new cv::Mat();
	tonemapMats_.lms = new cv::Mat();
	tonemapMats_.blend = new cv::Mat();
	tonemapMats_.lms2displayOutput = new cv::Mat();
	tonemapMats_.inputImage = new cv::Mat();

}

Percep_low_light::Percep_low_light(const TonemapVals& tvals) :tonemapVals_(tvals)
{
	// Initialize the Mats
	tonemapMats_.LMtx = new cv::Mat();
	tonemapMats_.MMtx = new cv::Mat();
	tonemapMats_.SMtx = new cv::Mat();
	tonemapMats_.RMtx = new cv::Mat();
	tonemapMats_.hdrImage = new cv::Mat();
	tonemapMats_.lms = new cv::Mat();
	tonemapMats_.blend = new cv::Mat();
	tonemapMats_.lms2displayOutput = new cv::Mat();
	tonemapMats_.inputImage = new cv::Mat();
}

Percep_low_light::~Percep_low_light()
{
	// Delete Mats
	delete tonemapMats_.LMtx;
	delete tonemapMats_.MMtx;
	delete tonemapMats_.SMtx;
	delete tonemapMats_.RMtx;
	delete tonemapMats_.hdrImage;
	delete tonemapMats_.lms;
	delete tonemapMats_.blend;
	delete tonemapMats_.lms2displayOutput;
	delete tonemapMats_.inputImage;
}


cv::Mat Percep_low_light::converto_Mat(QString fullPath)
{   
    //two operations : 1.imread for qstring 2.RGB->BGR//value type:32FC3
    cv::Mat exrimage = cv::imread(fullPath.toStdString(), cv::IMREAD_ANYCOLOR | cv::IMREAD_ANYDEPTH);
	cv::Mat inputimage;
	double minVal, maxVal;
	cv::minMaxLoc(exrimage, &minVal, &maxVal); //find minimum and maximum intensities
	//exrimage.convertTo(exrimage, CV_8UC3, 1.0 / (maxVal - minVal), -minVal * 1.0 / (maxVal - minVal)); //apply linear transformation
	exrimage = (exrimage - minVal) / (maxVal - minVal);
	tonemapMats_.inputImage = new cv::Mat(exrimage.rows, exrimage.cols, CV_8UC3);
	exrimage.convertTo(*tonemapMats_.inputImage, CV_8UC3, 255.0 / 1); // Assuming your .exr image data is normalized between 0 and 1. If not, replace 1 with maximum value in your exr image.
	inputimage = *tonemapMats_.inputImage;
	inputimage.convertTo(*tonemapMats_.hdrImage, CV_32FC3, tonemapVals_.exposure, 0.0);
	//*(tonemapMats_.hdrImage) = exrimage;
	//*(tonemapMats_.inputImage) = exrimage;
    //cv::cvtColor(exrimage, *(tonemapMats_.hdrImage), cv::COLOR_BGR2RGB);
	 
	//if the input is a .png/ldr image
	/*cv::Mat inputimage = cv::imread(fullPath.toStdString(), cv::IMREAD_COLOR);
	tonemapMats_.inputImage = new cv::Mat(inputimage.rows, inputimage.cols, CV_8UC3);
	*tonemapMats_.inputImage = inputimage;
	inputimage.convertTo(*tonemapMats_.hdrImage, CV_32FC3, tonemapVals_.exposure, 0.0);*/

    // note: if we need to process tone mapping on non-exr file, it needs to be extended, to scale ldr to hdr range.
    return *(tonemapMats_.hdrImage);
}

void Percep_low_light::HDRToLMSR(cv::Mat* hdrMat, cv::Mat* LMtx, cv::Mat* MMtx, cv::Mat* SMtx, cv::Mat* RMtx)
{
	float *buffer;
	int rows, cols;
	rows = hdrMat->rows;
	cols = hdrMat->cols;
	LMtx->create(rows, cols, CV_32FC1);
	MMtx->create(rows, cols, CV_32FC1);
	SMtx->create(rows, cols, CV_32FC1);
	RMtx->create(rows, cols, CV_32FC1);
	
	//read in buffers
	buffer = new float[rows * cols * 3];
	for (int h = 0; h < rows; h++)
	{
		for (int w = 0; w < cols; w++)
		{
			int pixIdx = w + (h * cols);
			buffer[pixIdx * 3 + 0] = hdrMat->at<cv::Vec3f>(h, w)[2];
			buffer[pixIdx * 3 + 1] = hdrMat->at<cv::Vec3f>(h, w)[1];
			buffer[pixIdx * 3 + 2] = hdrMat->at<cv::Vec3f>(h, w)[0];
		}
	}
	
	//create convert matrix
	cv::Mat ConvertMtx(4, 3, CV_32FC1);

	//matrix from cvx 
	ConvertMtx.at<float>(0, 0) = 0.000167179131;
	ConvertMtx.at<float>(0, 1) = 0.021378405283;
	ConvertMtx.at<float>(0, 2) = -0.000600300420;
	ConvertMtx.at<float>(1, 0) = 0.000037285102;
	ConvertMtx.at<float>(1, 1) = 0.017335413184;
	ConvertMtx.at<float>(1, 2) = -0.000607142696;
	ConvertMtx.at<float>(2, 0) = -0.000384788819;
	ConvertMtx.at<float>(2, 1) = 0.011063773501;
	ConvertMtx.at<float>(2, 2) = -0.000583754041;
	ConvertMtx.at<float>(3, 0) = -0.000211589221;
	ConvertMtx.at<float>(3, 1) = 0.018000447150;
	ConvertMtx.at<float>(3, 2) = -0.000621909939;

	//calculate the LMSR matrix
	cv::Mat LMSRMtx(4, rows * cols, CV_32FC1);

	cv::Mat LMSR(4, 1, CV_32FC1);
	cv::Mat RGB_mtx(3, 1, CV_32FC1);

	float minL = FLT_MAX;
	float minM = FLT_MAX;
	float minS = FLT_MAX;
	float minR = FLT_MAX;

	//float exposureMultiplier = 0;

	//test
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			int pixIdx = j + (i * cols);
			RGB_mtx.at<float>(0, 0) = buffer[pixIdx * 3 + 0];
			RGB_mtx.at<float>(1, 0) = buffer[pixIdx * 3 + 1];
			RGB_mtx.at<float>(2, 0) = buffer[pixIdx * 3 + 2];

			LMSR = ConvertMtx * RGB_mtx;
			LMtx->at<float>(i, j) = LMSR.at<float>(0, 0);
			MMtx->at<float>(i, j) = LMSR.at<float>(1, 0);
			SMtx->at<float>(i, j) = LMSR.at<float>(2, 0);
			RMtx->at<float>(i, j) = LMSR.at<float>(3, 0);
			if (LMtx->at<float>(i, j) < minL)
				minL = LMtx->at<float>(i, j);
			if (MMtx->at<float>(i, j) < minM)
				minM = MMtx->at<float>(i, j);
			if (SMtx->at<float>(i, j) < minS)
				minS = SMtx->at<float>(i, j);
			if (RMtx->at<float>(i, j) < minR)
				minR = RMtx->at<float>(i, j);
		}
	}

}

void Percep_low_light::purkinje(cv::Mat* LMtx, cv::Mat* MMtx, cv::Mat* SMtx, cv::Mat* RMtx, cv::Mat* I, cv::Mat* blend, TonemapVals* tvals)
{
	float lmax, mmax, smax, k1, k2, k5, k3PC, k3KC, rw, k6, p, y;
	lmax = .637;
	mmax = .392;
	smax = 1.606;

	//TODO: make these variables (figure out the range) 
	//determines the rod contribution to the neural signals
	k1 = tvals->kappa1; //.33;
	k2 = tvals->kappa2; //.5;

	k5 = tvals->z; //.25; //z
	k6 = tvals->rho3; //.4; //rho3
	k3PC = tvals->rho1; //.8; //rho1
	k3KC = tvals->rho4; //.6; //rho4
	rw = tvals->rho2; //.139; //rho2
	p = tvals->alpha; //.6189; //alpha
	y = tvals->y; //15.0; //y
	int rows = LMtx->rows;
	int cols = LMtx->cols;
	int r, c;

	I->create(rows, cols, CV_32FC3);
	blend->create(rows, cols, CV_32FC1);

	cv::Mat LCopy = LMtx->clone();
	cv::Mat MCopy = MMtx->clone();
	cv::Mat SCopy = SMtx->clone();

	for (r = 0; r < rows; r++)
	{
		for (c = 0; c < cols; c++)
		{
			LCopy.at<float>(r, c) = LMtx->at<float>(r, c) + k5 * RMtx->at<float>(r, c);
			MCopy.at<float>(r, c) = MMtx->at<float>(r, c) + k5 * RMtx->at<float>(r, c);
			SCopy.at<float>(r, c) = SMtx->at<float>(r, c) + k6 * RMtx->at<float>(r, c);
		}
	}

	cv::Mat* GL = new cv::Mat(rows, cols, CV_32FC1);
	cv::Mat* GM = new cv::Mat(rows, cols, CV_32FC1);
	cv::Mat* GS = new cv::Mat(rows, cols, CV_32FC1);
	G(&LCopy, k1, k2, GL);
	G(&MCopy, k1, k2, GM);
	G(&SCopy, k1, k2, GS);

	//opponent colors
	cv::Mat* M_L = new cv::Mat(rows, cols, CV_32FC1);
	cv::Mat* S_L_M = new cv::Mat(rows, cols, CV_32FC1);
	cv::Mat* LtM = new cv::Mat(rows, cols, CV_32FC1);
	float a = 15 * k5; //x
	float b = 1 + rw * k3PC; //the real rho1
	float d = k3PC + rw; //the real rho2

	for (r = 0; r < rows; r++)
	{
		for (c = 0; c < cols; c++)
		{
			M_L->at<float>(r, c) = a * (b * GM->at<float>(r, c) / mmax - d * GL->at<float>(r, c) / lmax);
			S_L_M->at<float>(r, c) = y * (k6 * GS->at<float>(r, c) / smax - \
				k3KC * (p * k5 * GL->at<float>(r, c) / lmax + (1 - p) * k5 * GM->at<float>(r, c) / mmax));
			LtM->at<float>(r, c) = 5 * k5 * (p * GL->at<float>(r, c) / lmax + (1 - p) * GM->at<float>(r, c) / mmax);

			blend->at<float>(r, c) = 1 - LtM->at<float>(r, c);
			if (blend->at<float>(r, c) < 0) blend->at<float>(r, c) = 0;

			//set opponent color values
			M_L->at<float>(r, c) = RMtx->at<float>(r, c) * M_L->at<float>(r, c) + tvals->redGreenVal;
			S_L_M->at<float>(r, c) = RMtx->at<float>(r, c) * S_L_M->at<float>(r, c) + tvals->blueYellowVal;
			LtM->at<float>(r, c) = RMtx->at<float>(r, c) * LtM->at<float>(r, c) + tvals->luminanceVal;

			//calculate lms values
			float value = (LtM->at<float>(r, c) - M_L->at<float>(r, c)) / 2;
			I->at<cv::Vec3f>(r, c)[0] = value;
			I->at<cv::Vec3f>(r, c)[1] = LtM->at<float>(r, c) - value;
			I->at<cv::Vec3f>(r, c)[2] = S_L_M->at<float>(r, c) + LtM->at<float>(r, c);

			I->at<cv::Vec3f>(r, c)[0] += LMtx->at<float>(r, c);
			I->at<cv::Vec3f>(r, c)[1] += MMtx->at<float>(r, c);
			I->at<cv::Vec3f>(r, c)[2] += SMtx->at<float>(r, c);
		}
	}
}

void Percep_low_light::lms2display(cv::Mat* lms, cv::Mat* outImage)
{
	int rows = lms->rows;
	int cols = lms->cols;
	outImage->create(rows, cols, CV_32FC3);

	std::cout << "TEST: rows " << rows << " cols" << cols << std::endl;
	cv::Mat colors(3, rows * cols, CV_32FC1);
	for (int r = 0; r < rows; r++) {
		for (int c = 0; c < cols; c++)
		{
			colors.at<float>(0, c + cols * r) = lms->at<cv::Vec3f>(r, c)[0];
			colors.at<float>(1, c + cols * r) = lms->at<cv::Vec3f>(r, c)[1];
			colors.at<float>(2, c + cols * r) = lms->at<cv::Vec3f>(r, c)[2];
		}
	}

	//monitorLMS matrix based on Adam's monitor settings
	cv::Mat monitorLMS(3, 3, CV_32FC1);
	monitorLMS.at<float>(0, 0) = 3.737001334053387;
	monitorLMS.at<float>(0, 1) = 4.474289710481138;
	monitorLMS.at<float>(0, 2) = 0.446763557036099;
	monitorLMS.at<float>(1, 0) = 1.260813115018955;
	monitorLMS.at<float>(1, 1) = 4.394344052644990;
	monitorLMS.at<float>(1, 2) = 0.587779102275466;
	monitorLMS.at<float>(2, 0) = 0.046245589359130;
	monitorLMS.at<float>(2, 1) = 0.244453440088045;
	monitorLMS.at<float>(2, 2) = 1.218879359722785;

	cv::Mat color(3, rows * cols, CV_32FC1);
	color = monitorLMS.inv() * colors;

	//reshape color matrix
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			//clamp the negative values to 0
			if (color.at<float>(0, j + cols * i) < 0.0f)  color.at<float>(0, j + cols * i) = 0.0f;
			if (color.at<float>(1, j + cols * i) < 0.0f)  color.at<float>(1, j + cols * i) = 0.0f;
			if (color.at<float>(2, j + cols * i) < 0.0f)  color.at<float>(2, j + cols * i) = 0.0f;
			outImage->at<cv::Vec3f>(i, j)[0] = color.at<float>(0, j + cols * i);
			outImage->at<cv::Vec3f>(i, j)[1] = color.at<float>(1, j + cols * i);
			outImage->at<cv::Vec3f>(i, j)[2] = color.at<float>(2, j + cols * i);
		}
	}
}

void Percep_low_light::tonemapping4lms(cv::Mat* img)
{
	int width = img->cols;
	int height = img->rows;

	cv::Mat output(height, width, CV_32FC3);

	for (int h = 0; h < height; h++)
	{
		for (int w = 0; w < width; w++)
		{
			output.at<cv::Vec3f>(h, w)[0] = img->at<cv::Vec3f>(h, w)[0];
			output.at<cv::Vec3f>(h, w)[1] = img->at<cv::Vec3f>(h, w)[1];
			output.at<cv::Vec3f>(h, w)[2] = img->at<cv::Vec3f>(h, w)[2];
		}
	}

	if (!img->empty()) {
		//write output as an exr and save it and then read it back in as a float* 
		cv::Mat outputexr;
		cv::cvtColor(output, outputexr, cv::COLOR_BGR2RGB);
		cv::imwrite("./BeforeDurandTonemapping.exr", outputexr);
		//cv::Mat beforeDurandExr = cv::imread("./BeforeDurandTonemapping.exr", cv::IMREAD_ANYCOLOR | cv::IMREAD_ANYDEPTH);
		cv::Ptr<cv::xphoto::TonemapDurand> durandTMO = cv::xphoto::createTonemapDurand();
		durandTMO->setGamma(2.2f);//2.2
		durandTMO->setSigmaSpace(0.02 * std::min(output.cols, output.rows));//0.02
		durandTMO->setContrast(50.0f);
		durandTMO->setSigmaColor(0.4f);//0.4
		durandTMO->setSaturation(1.0f);
		//apply Durand tone mapping
		cv::Mat afterDurandLdr;
		durandTMO->process(output, afterDurandLdr);
		afterDurandLdr *= 255;
		afterDurandLdr.convertTo(afterDurandLdr, CV_8UC3);
		cv::imwrite("./AfterDurandTonemapping.ppm", afterDurandLdr);

	} 
}

void Percep_low_light::reduceRange(cv::Mat* img, cv::Mat* alpha, int rows, int cols, const char* outputFileName, cv::Mat& LDRImg)
{
	assert(img->type() == CV_8UC3);

	float x = 0.1f;
	cv::Mat* Idouble = new cv::Mat(rows, cols, CV_32FC3);
	cv::Mat* MaxImage = new cv::Mat(rows, cols, CV_32FC1);
	cv::Mat* alphaScaling = new cv::Mat(rows, cols, CV_32FC1);
	cv::Mat* outputMtx = new cv::Mat(rows, cols, CV_8UC3);

	//for blending with the input.
	//Mat LDRImg = imread("Result/tonemappedHDR.ppm",1);  
	float alphaSum = 0;
	for (int r = 0; r < rows; r++) {
		for (int c = 0; c < cols; c++)
		{
			Idouble->at<cv::Vec3f>(r, c)[0] = ((double)img->at<cv::Vec3b>(r, c)[0]) / 256.0;
			Idouble->at<cv::Vec3f>(r, c)[1] = ((double)img->at<cv::Vec3b>(r, c)[1]) / 256.0;
			Idouble->at<cv::Vec3f>(r, c)[2] = ((double)img->at<cv::Vec3b>(r, c)[2]) / 256.0;
			MaxImage->at<float>(r, c) = 255.0;
			alphaScaling->at<float>(r, c) = (alpha->at<float>(r, c) * (1 - x)) + x;
			MaxImage->at<float>(r, c) *= alphaScaling->at<float>(r, c);

			outputMtx->at<cv::Vec3b>(r, c)[0] = (unsigned int)(Idouble->at<cv::Vec3f>(r, c)[0] * MaxImage->at<float>(r, c));
			outputMtx->at<cv::Vec3b>(r, c)[1] = (unsigned int)(Idouble->at<cv::Vec3f>(r, c)[1] * MaxImage->at<float>(r, c));
			outputMtx->at<cv::Vec3b>(r, c)[2] = (unsigned int)(Idouble->at<cv::Vec3f>(r, c)[2] * MaxImage->at<float>(r, c));

			//blending
			float colorBlendStart = 210.0;
			float intensity = 0.0, bright_factor = 0.0;
			intensity = 0.3 * (float)LDRImg.at<cv::Vec3b>(r, c)[0] + 0.6 * (float)LDRImg.at<cv::Vec3b>(r, c)[1] + 0.1 * (float)LDRImg.at<cv::Vec3b>(r, c)[2];
			if (intensity >= 210.0) {
				bright_factor = (intensity - colorBlendStart) / (float)(255 - colorBlendStart);
			}
			float finalBlendR = (1 - bright_factor) * alpha->at<float>(r, c) + bright_factor;
			float finalBlendG = (1 - bright_factor) * alpha->at<float>(r, c) + bright_factor;
			float finalBlendB = (1 - bright_factor) * alpha->at<float>(r, c) + bright_factor;

#if 0
			float purkinjeVsFixedR = ((float)LDRImg.at<Vec3b>(r, c)[0] - colorBlendStart) / (255.0f - colorBlendStart);
			float purkinjeVsFixedG = ((float)LDRImg.at<Vec3b>(r, c)[1] - colorBlendStart) / (255.0f - colorBlendStart);
			float purkinjeVsFixedB = ((float)LDRImg.at<Vec3b>(r, c)[2] - colorBlendStart) / (255.0f - colorBlendStart);
			if (purkinjeVsFixedR < 0) purkinjeVsFixedR = 0;
			if (purkinjeVsFixedG < 0) purkinjeVsFixedG = 0;
			if (purkinjeVsFixedB < 0) purkinjeVsFixedB = 0;
			float finalBlendR = (1 - purkinjeVsFixedR) * alpha->at<float>(r, c) + purkinjeVsFixedR * 1.0f;
			float finalBlendG = (1 - purkinjeVsFixedG) * alpha->at<float>(r, c) + purkinjeVsFixedG * 1.0f;
			float finalBlendB = (1 - purkinjeVsFixedB) * alpha->at<float>(r, c) + purkinjeVsFixedB * 1.0f;
#endif   
			outputMtx->at<cv::Vec3b>(r, c)[0] = (unsigned int)(finalBlendR * LDRImg.at<cv::Vec3b>(r, c)[0] + ((1 - finalBlendR) * outputMtx->at<cv::Vec3b>(r, c)[0]));
			outputMtx->at<cv::Vec3b>(r, c)[1] = (unsigned int)(finalBlendG * LDRImg.at<cv::Vec3b>(r, c)[1] + ((1 - finalBlendG) * outputMtx->at<cv::Vec3b>(r, c)[1]));
			outputMtx->at<cv::Vec3b>(r, c)[2] = (unsigned int)(finalBlendB * LDRImg.at<cv::Vec3b>(r, c)[2] + ((1 - finalBlendB) * outputMtx->at<cv::Vec3b>(r, c)[2]));
			alphaSum += alpha->at<float>(r, c);
		}
	}

	float alphaMean = alphaSum / (rows * cols);
	printf("The average blend value is %f\n", alphaMean);
	cv::imwrite(outputFileName, *outputMtx);
}


