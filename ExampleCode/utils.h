#ifndef __UTILS_H__
#define __UTILS_H__

#include <vector>
#include <cmath>
#include "loadPNM.h"

#define PI 3.14159265358979323
#define LOW_IGNORE_THRESH 0.005
#define HIGH_IGNORE_THRESH 0.920
#define TWO_STOP 4
#define GAMMA 1.5
#define N 2
#define RADIUS 255
#define DIAMETER 511
#define DIMENSION 3
#include <iostream>

/*//////////////////////////////////////////////////////////////////////////
Util functions
//////////////////////////////////////////////////////////////////////////*/

bool isInCircle(unsigned int width, unsigned int height);
float getX(int width);
float getY(int height);

//Center wighting function
float w(float x);

vector<float> getSurfaceNormal(float x, float y);
vector<float> getReflectanceVector(vector<float>& normal, vector<float>& v);
float findMaxIntensity(float* image_in, unsigned int width, unsigned int height, unsigned int numComponents);

// Tone Maps a PFM image into range 0-1
void toneMapper(float* image_in, unsigned int width, unsigned int height, unsigned int numComponents);

char doNothing(float value);
void gammaFunc(float gamma, float* image_in, unsigned int width, unsigned int height, unsigned int numComponents);

void NExposureScale(int n, float* image_in, unsigned int width, unsigned int height, unsigned int numComponents);

char reflectanceToPPM(float value);

void applyFunctionOnAllPixelsPFM(vector<float*> images_in, float* image_out, 
    unsigned int width, unsigned int height, unsigned int numComponents, float (*func)(vector<float*>, unsigned int));

void applyFunctionOnAllPixelsPPM(vector<unsigned char*> images_in, 
    unsigned char* image_out,
    unsigned int width, unsigned int height, unsigned int numComponents, 
    char (*func)(vector<unsigned char*>, unsigned int));

void applyFunctionOnAllPixelsPPMFromPFM(float* image_in, 
    unsigned char* image_out,
    unsigned int width, unsigned int height, unsigned int numComponents, 
    char (*func)(float));

char exposureGamma(float value);

void LoadPPMAndSavePFM(const char *image_in, const char *image_out);

void LoadPFMSavePPM(const char *image_in, const char *image_out, char (*func)(float));

void reflectanceCircleAndSavePFM(const char *image_out);

void reflectanceSphereSavePPM(const char *image_in, const char *image_out);

void mapLatLongToSphere(const char* reflectance, const char* latLongMap, 
    const char* image_out);

float returnHDRCoponentPFM(vector<float*> imgs_in, unsigned int index);

void processHDRAndSavePFM(vector<const char*> images_in, const char *image_out);
#endif
