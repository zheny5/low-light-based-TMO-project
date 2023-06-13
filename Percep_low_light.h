#pragma once
#include <opencv2/opencv.hpp>
class Percep_low_light
{
public:
	typedef struct
	{
		float exposure = 64.0;
		float redGreenVal = 0.0;
		float blueYellowVal = 0.0;
		float luminanceVal = 0.0;
		boolean preview = 1;
		float kappa1 = .33; //rod contribution
		float kappa2 = .5; //rod contribution
		float rho1 = .8; //weight of medium cones
		float rho2 = .139; //weight of long cones
		float rho3 = .4; //weight of short cones
		float rho4 = .6; //weight of long and med cones
		float alpha = .6189; //blending ratio btw long and med cones
		float y = 15.0; //scale for blue/yellow channel
		float z = .24; //scale for luminance channel
	}TonemapVals;

	typedef struct
	{
		cv::Mat* LMtx = nullptr;
		cv::Mat* MMtx = nullptr;
		cv::Mat* SMtx = nullptr;
		cv::Mat* RMtx = nullptr;
		cv::Mat* hdrImage = nullptr;
		cv::Mat* inputImage = nullptr;
		cv::Mat* lms = nullptr;
		cv::Mat* blend = nullptr;
		cv::Mat* lms2displayOutput = nullptr;
	}TonemapMats;

	TonemapVals tonemapVals_;
	TonemapMats tonemapMats_;

public:
	Percep_low_light();
	Percep_low_light(const TonemapVals& tvals);
	~Percep_low_light();
public:
	cv::Mat converto_Mat(QString fullPath);
	void HDRToLMSR(cv::Mat* hdrMat, cv::Mat* LMtx, cv::Mat* MMtx, cv::Mat* SMtx, cv::Mat* RMtx);
	void purkinje(cv::Mat* LMtx, cv::Mat* MMtx, cv::Mat* SMtx, cv::Mat* RMtx, cv::Mat* I, cv::Mat* blend, TonemapVals* tvals);
	void lms2display(cv::Mat* lms, cv::Mat* outImage);
	void tonemapping4lms(cv::Mat* img);
	void reduceRange(cv::Mat* img, cv::Mat* alpha, int rows, int cols, const char* outputFileName, cv::Mat& LDRImg);
};
