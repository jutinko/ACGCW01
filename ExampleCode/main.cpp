/*//////////////////////////////////////////////////////////////////////////
Author: Abhijeet Ghosh
Year: 2013
//////////////////////////////////////////////////////////////////////////*/
#include <cstdio>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include "loadPNM.h"

#define PI 3.14159265358979323
#define LOW_IGNORE_THRESH 0.005
#define HIGH_IGNORE_THRESH 0.920
#define TWO_STOP 4
#define EXPOSUREFACTOR 10
#define GAMMA 2.2
#define RADIUS 255
#define DIAMETER 511
#define DIMENSION 3
#define uint unsigned int

#include <iostream>

using namespace std;

unsigned int width;
unsigned int height;
unsigned int numComponents;

bool isInCircle(uint width, uint height)
{
  uint x = RADIUS-width > 0 ? RADIUS-width : width-RADIUS;
  uint y = RADIUS-height > 0 ? RADIUS-height : height-RADIUS;
  return (x*x+y*y < RADIUS*RADIUS);
}

float getX(int width)
{
  return (float)(width-RADIUS)/RADIUS;
}

float getY(int height)
{
  return (float)(height-RADIUS)/RADIUS;
}

//Center wighting function
float w(float x) 
{
  //return sin(PI*x);
  if(x <= 0.5)
  {
    return 2*x;
  }
  return -2*x+2;
}

float findMaxIntensity(float* image_in, uint width, uint height, uint numComponents)
{
  float result = -1;
  float currValue;
  for(uint i = 0; i < height; ++i)
  {
    for(uint j = 0; j < width; ++j )
    {
      currValue = 0;
      for(uint k = 0; k < numComponents; ++k)
      {
        uint index = i*width*numComponents + j*numComponents + k; //index within the image
        currValue += image_in[index];
      }
      if(currValue/numComponents > result)
      {
        result = currValue/numComponents;
      }
    }
  }
  return result;
}

// Tone Maps a PFM image into range 0-1
void toneMapper(float* image_in, uint width, uint height, uint numComponents)
{
  float maxIntensity = findMaxIntensity(image_in, width, height, numComponents);
  for ( uint i = 0 ; i < height ; ++i ) // height
  {
    for ( uint j = 0 ; j < width ; ++j ) // width
    {
      for ( uint k = 0 ; k < numComponents ; ++k ) // color channels - 3 for RGB images
      {
        uint index = i*width*numComponents + j*numComponents + k; //index within the image
        image_in[index] /= maxIntensity;
      }
    }
  }
}

void NExposureScale(int n, float* image_in, uint width, uint height, uint numComponents)
{
  int stops = pow(2, n);
  for ( uint i = 0 ; i < height ; ++i ) // height
  {
    for ( uint j = 0 ; j < width ; ++j ) // width
    {
      for ( uint k = 0 ; k < numComponents ; ++k ) // color channels - 3 for RGB images
      {
        uint index = i*width*numComponents + j*numComponents + k; //index within the image
        float scaledValue = image_in[index]*stops;
        image_in[index] = scaledValue > 1.0f ? 1.0f : scaledValue;
      }
    }
  }
}

void gammaFunc(float gamma, float* image_in, uint width, uint height, uint numComponents)
{
  for ( uint i = 0 ; i < height ; ++i ) // height
  {
    for ( uint j = 0 ; j < width ; ++j ) // width
    {
      for ( uint k = 0 ; k < numComponents ; ++k ) // color channels - 3 for RGB images
      {
        uint index = i*width*numComponents + j*numComponents + k; //index within the image
        image_in[index] = pow(image_in[index], 1/gamma);
      }
    }
  }
}

void applyFunctionOnAllPixelsPFM(vector<float*> images_in, float* image_out, 
    uint width, uint height, uint numComponents, float (*func)(vector<float*>, uint))
{
  for ( uint i = 0 ; i < height ; ++i ) // height
  {
    for ( uint j = 0 ; j < width ; ++j ) // width
    {
      for ( uint k = 0 ; k < numComponents ; ++k ) // color channels - 3 for RGB images
      {
        uint index = i*width*numComponents + j*numComponents + k; //index within the image
        image_out[index] = func(images_in, index);
      }
    }
  }
}

void applyFunctionOnAllPixelsPPM(vector<unsigned char*> images_in, 
    unsigned char* image_out,
    uint width, uint height, uint numComponents, 
    char (*func)(vector<unsigned char*>, uint))
{
  for ( uint i = 0 ; i < height ; ++i ) // height
  {
    for ( uint j = 0 ; j < width ; ++j ) // width
    {
      for ( uint k = 0 ; k < numComponents ; ++k ) // color channels - 3 for RGB images
      {
        uint index = i*width*numComponents + j*numComponents + k; //index within the image
        image_out[index] = func(images_in, index);
      }
    }
  }
}

void applyFunctionOnAllPixelsPPMFromPFM(vector<float*> images_in, 
    unsigned char* image_out,
    uint width, uint height, uint numComponents, 
    char (*func)(vector<float*>, uint))
{
  for ( uint i = 0 ; i < height ; ++i ) // height
  {
    for ( uint j = 0 ; j < width ; ++j ) // width
    {
      for ( uint k = 0 ; k < numComponents ; ++k ) // color channels - 3 for RGB images
      {
        uint index = i*width*numComponents + j*numComponents + k; //index within the image
        if(isInCircle(j, i))
        {
          image_out[index] = func(images_in, index);
        }
      }
    }
  }
}

void LoadPPMAndSavePFM(const char *image_in, const char *image_out)
{
  unsigned char *img_in = loadPNM(image_in, width, height, numComponents);
  float *img_out = new float [width*height*numComponents];

  for ( uint i = 0 ; i < height ; ++i ) // height
  {
    for ( uint j = 0 ; j < width ; ++j ) // width
    {
      for ( uint k = 0 ; k < numComponents ; ++k ) // color channels - 3 for RGB images
      {
        uint index = i*width*numComponents + j*numComponents + k; //index within the image

        //typecast 0 - RADIUS values to the 0.0f -> 1.0f range 
        img_out[index] = static_cast<float>(img_in[index])/255.0f; //typecast all color channels of each pixel
                
      }
    }
  }
  WritePFM(image_out, width, height, numComponents, img_out);
}

void LoadPFMAndSavePPM(const char *image_in, const char *image_out)
{
  float* img_in = loadPFM(image_in, width, height, numComponents);
  //toneMapper(img_in, width, height, numComponents);
  //NExposureScale(EXPOSUREFACTOR, img_in, width, height, numComponents);
  //gammaFunc(GAMMA, img_in, width, height, numComponents);
  unsigned char *img_out = new unsigned char [width*height*numComponents];

  for ( uint i = 0 ; i < height ; ++i ) // height
  {
    for ( uint j = 0 ; j < width ; ++j ) // width
    {
      for ( uint k = 0 ; k < numComponents ; ++k ) // color channels - 3 for RGB images
      {
        uint index = i*width*numComponents + j*numComponents + k; //index within the image
        //typecast 0.0f -> 1.0f values to the 0 - 255 range 
        img_out[index] = static_cast<unsigned char>(img_in[index]*255.0f); //typecast all color channels of each pixel
                
      }
    }
  }
  WritePNM(image_out, width, height, numComponents, img_out);
}

vector<float> getSurfaceNormal(float x, float y)
{
  vector<float> normal;
  float z = sqrt(1-(x*x+y*y));
  normal.push_back(x);
  normal.push_back(y);
  normal.push_back(z);
  return normal;
}

vector<float> getReflectanceVector(vector<float> normal, vector<float> v)
{
  float nDotV = 0;
  for(int i = 0; i < DIMENSION; ++i)
  {
    nDotV += normal[i]*v[i];
  }

  vector<float> result;
  for(int i = 0; i < DIMENSION; ++i)
  {
    result.push_back(2*nDotV*normal[i]-v[i]);
  }
  return result;
}

void CreateAndSavePFM(const char *image_out)
{
  width = DIAMETER; // set size of image to DIAMETERxDIAMETER pixels
  height = DIAMETER;
  numComponents = 3;
  
  float *img_out = new float [width*height*numComponents];

  vector<float> normal;
  vector<float> r;
  vector<float> v;
  v.push_back(0.0f);
  v.push_back(0.0f);
  v.push_back(1.0f);
  for ( uint i = 0 ; i < height ; ++i ) // height
  {
    for ( uint j = 0 ; j < width ; ++j ) // width
    {
      normal = getSurfaceNormal(getX(j), getY(i));
      r = getReflectanceVector(normal, v);
      for ( uint k = 0 ; k < numComponents ; ++k ) // color channels - 3 for RGB images
      {
        uint index = i*width*numComponents + j*numComponents + k; //index within the image
        img_out[index] = isInCircle(i, j) ? r[k] : 0.0f;
      }
    }
  }
  WritePFM(image_out, width, height, numComponents, img_out);
  delete img_out;
}

char returnSameValuePPM(vector<unsigned char*> imgs_in, uint index)
{
  return imgs_in[0][index];
}

char reflectanceToPPM(vector<float*> imgs_in, uint index)
{
  return (imgs_in[0][index]+1.0f)/2*255;
}

void reflectanceSphereSavePPM(const char *image_in, const char *image_out)
{
  float* img_in = loadPFM(image_in, width, height, numComponents);
  unsigned char *img_out = new unsigned char [width*height*numComponents];
  vector<float*> imgs_in;
  imgs_in.push_back(img_in);
  applyFunctionOnAllPixelsPPMFromPFM(imgs_in, img_out, width, height, numComponents, reflectanceToPPM);
  WritePNM(image_out, width, height, numComponents, img_out);
  delete img_out;
}

float returnSameValuePFM(vector<float*> imgs_in, uint index)
{
  return imgs_in[0][index];
}

void mapLatLongToSphere(const char* reflectance, const char* latLongMap, 
    const char* image_out)
{
  uint widthRef, heightRef, numComponentsRef;
  uint widthLat, heightLat, numComponentsLat;

  float* reflectance_f = loadPFM(reflectance, widthRef, heightRef, numComponentsRef);
  float* latLongMap_f = loadPFM(latLongMap, widthLat, heightLat, numComponentsLat);
  float* img_out = new float[widthRef*heightRef*numComponentsRef];

  for ( uint i = 0 ; i < heightRef ; ++i ) // height
  {
    for ( uint j = 0 ; j < widthRef ; ++j ) // width
    {
      if(isInCircle(j, i))
      {
        vector<float> xyz;
        for ( uint k = 0 ; k < numComponentsRef ; ++k ) // color channels - 3 for RGB images
        {
          uint index = i*widthRef*numComponentsRef + j*numComponentsRef + k; 
          xyz.push_back(reflectance_f[index]);
        }
        float theta = acos(xyz[1])/PI;
        float phi = (atan2(xyz[2], xyz[0])+PI)/(2*PI);
        uint mappedW = phi*widthLat;
        uint mappedH = theta*heightLat;
        float value;
        for ( uint k = 0 ; k < numComponentsRef ; ++k ) 
        {
          uint index = i*widthRef*numComponentsRef + j*numComponentsRef + k; 

          value = latLongMap_f[mappedH*widthLat*numComponentsLat+
            mappedW*numComponentsLat+k];
          img_out[index] = value > 1.0f ? 1.0f : value;
        }
      }
    }
  }
  cout << "befoer triple forloop" << endl;

  WritePFM(image_out, widthRef, heightRef, numComponentsRef, img_out);
  delete img_out;
}

//This method assumes that the order of images in the given vector
//is in ascedning order of exposure time
float returnHDRCoponentPFM(vector<float*> imgs_in, uint index)
{
  float result = 0;
  float exponent = 0;
  float sumWeightedValue = 0;
  float currValue;
  float deltaT = 1;
  vector<float*>::iterator i = imgs_in.begin();
  while(i != imgs_in.end())
  {
    currValue = (*i)[index];
    if(currValue < LOW_IGNORE_THRESH || currValue > HIGH_IGNORE_THRESH)
    {
      ++i;
      deltaT *= TWO_STOP;
      continue;
    }
    exponent += log(currValue/deltaT)*w(currValue);
    deltaT *= TWO_STOP;
    sumWeightedValue += w(currValue);
    ++i;
  }

  result = exp(exponent/sumWeightedValue);
  return result;
}

void processHDRAndSavePFM(vector<const char*> images_in, const char *image_out)
{
  vector<float*> imgs_in;
  for(vector<const char*>::iterator i = images_in.begin(); i != images_in.end();
      ++i)
  {
    imgs_in.push_back(loadPFM(*i, width, height, numComponents));
  }
  float *img_out = new float [width*height*numComponents];
  applyFunctionOnAllPixelsPFM(imgs_in, img_out, width, height, numComponents, returnHDRCoponentPFM);

  WritePFM(image_out, width, height, numComponents, img_out);

  delete img_out;
}

int main(int argc, char** argv)
{
  cerr<<"main invoked: arguments - <image_out (.pfm)> "<<endl;
  cerr<<"main invoked: arguments - <image_in (.ppm)> <image_out (.ppm)> "<<endl;
  
  if(argc < 2)
  {
    cout<<"Too few arguments! .... exiting."<<endl;
    return 0;
  }        

  //if(argc == 2)
  //{  
    //mapLatLongToSphere(argv[1], argv[2], argv[3]); //Creates and saves a PFM
    LoadPFMAndSavePPM(argv[1], argv[2]);
  //} 
 // else 
 // {  
 //   vector<const char*> imgs_in;
 //   //The last parameter is the output file
 //   for(int i = 1; i < argc-1; ++i)
 //   {
 //     imgs_in.push_back(argv[i]);
 //   }
 //   processHDRAndSavePFM(imgs_in, argv[argc-1]); //Loads and saves a PFM file
 //   const char* outName = "ppmout.ppm";
 //   LoadPFMAndSavePPM(argv[argc-1], outName);
 // }  
  return 0;
}
