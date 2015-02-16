/*//////////////////////////////////////////////////////////////////////////
Author: Abhijeet Ghosh
Year: 2013
//////////////////////////////////////////////////////////////////////////*/
#include "utils.h"

using namespace std;

int main(int argc, char** argv)
{
  if(argc < 2)
  {
    cout<<"Too few arguments! .... exiting."<<endl;
    return 0;
  }        

  //reflectanceCircleAndSavePFM(argv[1]);
  //cout << "her" << endl;
 //mapLatLongToSphere(argv[1], argv[2], argv[3]); //Creates and saves a PFM
  LoadPFMSavePPM(argv[1], argv[2], exposureGamma);
  //  vector<const char*> imgs_in;
  //  //The last parameter is the output file
  //  for(int i = 1; i < argc-1; ++i)
  //  {
  //    cout << argv[i] << endl;
  //    imgs_in.push_back(argv[i]);
  //  }
  //  processHDRAndSavePFM(imgs_in, argv[argc-1]); //Loads and saves a PFM file
 //   const char* outName = "ppmout.ppm";
 //   LoadPFMAndSavePPM(argv[argc-1], outName);
 // }  
  return 0;
}
