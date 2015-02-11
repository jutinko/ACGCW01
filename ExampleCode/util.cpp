#include "utils.h"

bool isInCircle(unsigned int width, unsigned int height)
{
  unsigned int x = RADIUS-width > 0 ? RADIUS-width : width-RADIUS;
  unsigned int y = RADIUS-height > 0 ? RADIUS-height : height-RADIUS;
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
  if(x <= 0.5)
  {
    return 2*x;
  }
  return -2*x+2;
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

float findMaxIntensity(float* image_in, unsigned int width, unsigned int height, unsigned int numComponents)
{
  float result = -1;
  float currValue;
  for(unsigned int i = 0; i < height; ++i)
  {
    for(unsigned int j = 0; j < width; ++j)
    {
      currValue = 0;
      for(unsigned int k = 0; k < numComponents; ++k)
      {
        unsigned int index = i*width*numComponents + j*numComponents + k;
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
void toneMapper(float* image_in, unsigned int width, unsigned int height, unsigned int numComponents)
{
  float maxIntensity = findMaxIntensity(image_in, width, height, numComponents);
  for(unsigned int i = 0; i < height; ++i)
  {
    for(unsigned int j = 0; j < width; ++j)
    {
      for(unsigned int k = 0; k < numComponents; ++k)
      {
        unsigned int index = i*width*numComponents + j*numComponents + k;
        image_in[index] /= maxIntensity;
      }
    }
  }
}

float gammaFunc(float value)
{
  return pow(value, 1/GAMMA);
}

float NExposureScale(float value)
{
  int stops = pow(2, N);
  float result = value*stops;
  return  result> 1.0f ? 1.0f : result;
}

void applyFunctionOnAllPixelsPFM(vector<float*> images_in, float* image_out, 
    unsigned int width, unsigned int height, unsigned int numComponents, float (*func)(vector<float*>, unsigned int))
{
  for(unsigned int i = 0; i < height; ++i)
  {
    for(unsigned int j = 0; j < width; ++j)
    {
      for(unsigned int k = 0; k < numComponents; ++k)
      {
        unsigned int index = i*width*numComponents + j*numComponents + k;
        image_out[index] = func(images_in, index);
      }
    }
  }
}

void applyFunctionOnAllPixelsPPM(vector<unsigned char*> images_in, 
    unsigned char* image_out,
    unsigned int width, unsigned int height, unsigned int numComponents, 
    char (*func)(vector<unsigned char*>, unsigned int))
{
  for(unsigned int i = 0; i < height; ++i)
  {
    for(unsigned int j = 0; j < width; ++j)
    {
      for(unsigned int k = 0; k < numComponents; ++k)
      {
        unsigned int index = i*width*numComponents + j*numComponents + k;
        image_out[index] = func(images_in, index);
      }
    }
  }
}

void applyFunctionOnAllPixelsPPMFromPFM(float* image_in, 
    unsigned char* image_out,
    unsigned int width, unsigned int height, unsigned int numComponents, 
    char (*func)(float))
{
  for(unsigned int i = 0; i < height; ++i)
  {
    for(unsigned int j = 0; j < width; ++j)
    {
      for(unsigned int k = 0; k < numComponents; ++k)
      {
        unsigned int index = i*width*numComponents + j*numComponents + k;
        image_out[index] = func(image_in[index]);
      }
    }
  }
}

char exposureGamma(float value)
{
  return gammaFunc(NExposureScale(value))*255.0f;
}

char doNothing(float value)
{
  return value*255.0f;
}

char reflectanceToPPM(float value)
{
  return (value+1.0f)/2*255;
}

void LoadPFMSavePPM(const char *image_in, const char *image_out, char (*func)(float))
{
  unsigned int width, height, numComponents;
  float* img_in = loadPFM(image_in, width, height, numComponents);
  //toneMapper(img_in, width, height, numComponents);
  unsigned char *img_out = new unsigned char [width*height*numComponents];
  applyFunctionOnAllPixelsPPMFromPFM(img_in, img_out, width, height, numComponents, func);
  WritePNM(image_out, width, height, numComponents, img_out);
  delete img_out;
}

void reflectanceCircleAndSavePFM(const char *image_out)
{
  unsigned int width, height, numComponents;

  width = DIAMETER; 
  height = DIAMETER;
  numComponents = DIMENSION;
  
  float *img_out = new float [width*height*numComponents];

  vector<float> normal;
  vector<float> r;
  vector<float> v;
  v.push_back(0.0f);
  v.push_back(0.0f);
  v.push_back(1.0f);
  for(unsigned int i = 0; i < height; ++i)
  {
    for(unsigned int j = 0; j < width; ++j)
    {
      normal = getSurfaceNormal(getX(j), getY(i));
      r = getReflectanceVector(normal, v);
      for(unsigned int k = 0; k < numComponents; ++k)
      {
        unsigned int index = i*width*numComponents + j*numComponents + k;
        img_out[index] = isInCircle(i, j) ? r[k] : 0.0f;
      }
    }
  }
  WritePFM(image_out, width, height, numComponents, img_out);
  delete img_out;
}


void mapLatLongToSphere(const char* reflectance, const char* latLongMap, 
    const char* image_out)
{
  unsigned int widthRef, heightRef, numComponentsRef;
  unsigned int widthLat, heightLat, numComponentsLat;

  float* reflectance_f = loadPFM(reflectance, widthRef, heightRef, numComponentsRef);
  float* latLongMap_f = loadPFM(latLongMap, widthLat, heightLat, numComponentsLat);
  float* img_out = new float[widthRef*heightRef*numComponentsRef];

  for(unsigned int i = 0; i < heightRef; ++i)
  {
    for(unsigned int j = 0; j < widthRef; ++j)
    {
      if(isInCircle(j, i))
      {
        vector<float> xyz;
        for(unsigned int k = 0; k < numComponentsRef; ++k)
        {
          unsigned int index = i*widthRef*numComponentsRef + j*numComponentsRef + k; 
          xyz.push_back(reflectance_f[index]);
        }
        float theta = acos(xyz[1])/PI;
        float phi = (atan2(xyz[2], xyz[0])+PI)/(2*PI);
        unsigned int mappedW = phi*widthLat;
        unsigned int mappedH = theta*heightLat;
        float value;
        for(unsigned int k = 0; k < numComponentsRef; ++k) 
        {
          unsigned int index = i*widthRef*numComponentsRef + j*numComponentsRef + k; 

          value = latLongMap_f[mappedH*widthLat*numComponentsLat+
            mappedW*numComponentsLat+k];
          img_out[index] = value > 1.0f ? 1.0f : value;
        }
      }
    }
  }

  WritePFM(image_out, widthRef, heightRef, numComponentsRef, img_out);
  delete img_out;
}

//This method assumes that the order of images in the given vector
//is in ascedning order of exposure time
float returnHDRCoponentPFM(vector<float*> imgs_in, unsigned int index)
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
  unsigned int width, height, numComponents;
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
