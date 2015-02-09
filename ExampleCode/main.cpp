/*//////////////////////////////////////////////////////////////////////////
Author: Abhijeet Ghosh
Year: 2013
//////////////////////////////////////////////////////////////////////////*/
#include <cstdio>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include "loadPNM.h"

#define PI 3.14159265358979323
#define LOW_IGNORE_THRESH 0.005
#define HIGH_IGNORE_THRESH 0.920
#define TWO_STOP 4
#define uint unsigned int

#include <iostream>

using namespace std;

unsigned int width;
unsigned int height;
unsigned int numComponents;

//Center wighting function
float w(float x) 
{
  //return sin(PI*x);
  if(x <= 0.5)
    return 2*x;
  return -2*x+2;
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

void applyFunctionOnAllPixelsPPM(vector<unsigned char*> images_in, unsigned char* image_out, 
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

				//typecast 0 - 255 values to the 0.0f -> 1.0f range 
				img_out[index] = static_cast<float>(img_in[index])/255.0f; //typecast all color channels of each pixel
								
			}
		}
	}
	WritePFM(image_out, width, height, numComponents, img_out);
}

void LoadPFMAndSavePPM(const char *image_in, const char *image_out)
{
	float *img_in = loadPFM(image_in, width, height, numComponents);
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

void CreateAndSavePFM(const char *image_out)
{
	width = 511; // set size of image to 511x511 pixels
	height = 511;
	numComponents = 3;
	
	float *img_out = new float [width*height*numComponents];

	for ( uint i = 0 ; i < height ; ++i ) // height
  {
		for ( uint j = 0 ; j < width ; ++j ) // width
		{
			for ( uint k = 0 ; k < numComponents ; ++k ) // color channels - 3 for RGB images
			{
				uint index = i*width*numComponents + j*numComponents + k; //index within the image

				//set image to white
				img_out[index] = 1.0f; //RGB all set to white
			}
		}
	}
	WritePFM(image_out, width, height, numComponents, img_out);
}

char returnSameValuePPM(vector<unsigned char*> imgs_in, uint index)
{
  return imgs_in[0][index];
}

void LoadAndSavePPM(const char *image_in, const char *image_out)
{
	unsigned char *img_in = loadPNM(image_in, width, height, numComponents);
	unsigned char *img_out = new unsigned char [width*height*numComponents];
  vector<unsigned char*> imgs_in;
  imgs_in.push_back(img_in);
  applyFunctionOnAllPixelsPPM(imgs_in, img_out, width, height, numComponents, returnSameValuePPM);
	WritePNM(image_out, width, height, numComponents, img_out);
  delete img_out;
}

float returnSameValuePFM(vector<float*> imgs_in, uint index)
{
  return imgs_in[0][index];
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

void LoadAndSavePFM(vector<const char*> images_in, const char *image_out)
{
  vector<float*> imgs_in;
  for(vector<const char*>::iterator i = images_in.begin(); i != images_in.end();
      ++i)
  {
    imgs_in.push_back(loadPFM(*i, width, height, numComponents));
  }
	//float *img_in = loadPFM(image_in, width, height, numComponents);
	float *img_out = new float [width*height*numComponents];
  //imgs_in.push_back(img_in);
  //applyFunctionOnAllPixelsPFM(imgs_in, img_out, width, height, numComponents, returnSameValuePFM);
      cout << "befeore apply fdsjak;e" << endl;
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

  if(argc == 2)
  {  
    CreateAndSavePFM(argv[1]); //Creates and saves a PFM
  } else 
  {  
    vector<const char*> imgs_in;
    //The last parameter is the output file
    for(int i = 1; i < argc-1; ++i)
    {
      imgs_in.push_back(argv[i]);
    }
    LoadAndSavePFM(imgs_in, argv[argc-1]); //Loads and saves a PFM file
  }  
  return 0;
}
