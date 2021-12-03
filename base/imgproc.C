
#include <cmath>
#include "imgproc.h"

#include <iostream>

#include <OpenImageIO/imageio.h>
OIIO_NAMESPACE_USING

using namespace img;


ImgProc::ImgProc() :
  Nx (0),
  Ny (0),
  Nc (0),
  Nsize (0),
  img_data (nullptr)
{}


ImgProc::~ImgProc()
{
   clear();
}
   

void ImgProc::clear()
{
   if( img_data != nullptr ){ delete[] img_data; img_data = nullptr;}
   Nx = 0;
   Ny = 0;
   Nc = 0;
   Nsize = 0;
}


void ImgProc::clear(int nX, int nY, int nC)
{
   clear();
   Nx = nX;
   Ny = nY;
   Nc = nC;
   Nsize = (long)Nx * (long)Ny * (long)Nc;
   img_data = new float[Nsize];
#pragma omp parallel for
   for(long i=0;i<Nsize;i++){ img_data[i] = 0.0; }
}


bool ImgProc::load( const std::string& filename )
{
   auto in = ImageInput::create (filename);
   if (!in) {return false;}
   ImageSpec spec;
   in->open (filename, spec);
   clear();
   Nx = spec.width;
   Ny = spec.height;
   Nc = spec.nchannels;
   Nsize = (long)Nx * (long)Ny * (long)Nc;
   img_data = new float[Nsize];
   in->read_image(TypeDesc::FLOAT, img_data);
   in->close ();
   return true;
}


void ImgProc::value( int i, int j, std::vector<float>& pixel) const
{
   pixel.clear();
   if( img_data == nullptr ){ return; }
   if( i<0 || i>=Nx ){ return; }
   if( j<0 || j>=Ny ){ return; }
   pixel.resize(Nc);
   for( int c=0;c<Nc;c++ )
   {
      pixel[c] = img_data[index(i,j,c)];
   }
   return;
}


void ImgProc::set_value( int i, int j, const std::vector<float>& pixel)
{
   if( img_data == nullptr ){ return; }
   if( i<0 || i>=Nx ){ return; }
   if( j<0 || j>=Ny ){ return; }
   if( Nc > (int)pixel.size() ){ return; }
#pragma omp parallel for
   for( int c=0;c<Nc;c++ )
   {
      img_data[index(i,j,c)] = pixel[c];
   }
   return;
}


ImgProc::ImgProc(const ImgProc& v) :
  Nx (v.Nx),
  Ny (v.Ny),
  Nc (v.Nc),
  Nsize (v.Nsize)
{
   img_data = new float[Nsize];
#pragma omp parallel for
   for( long i=0;i<Nsize;i++){ img_data[i] = v.img_data[i]; }
}


ImgProc& ImgProc::operator=(const ImgProc& v)
{
   if( this == &v ){ return *this; }
   if( Nx != v.Nx || Ny != v.Ny || Nc != v.Nc )
   {
      clear();
      Nx = v.Nx;
      Ny = v.Ny;
      Nc = v.Nc;
      Nsize = v.Nsize;
   }
   img_data = new float[Nsize];
#pragma omp parallel for
   for( long i=0;i<Nsize;i++){ img_data[i] = v.img_data[i]; }
   return *this;
}


void ImgProc::operator*=(float v)
{
   if( img_data == nullptr ){ return; }
#pragma omp parallel for
   for( long i=0;i<Nsize;i++ ){ img_data[i] *= v; }
}


void ImgProc::operator/=(float v)
{
   if( img_data == nullptr ){ return; }
#pragma omp parallel for
   for( long i=0;i<Nsize;i++ ){ img_data[i] /= v; }
}


void ImgProc::operator+=(float v)
{
   if( img_data == nullptr ){ return; }
#pragma omp parallel for
   for( long i=0;i<Nsize;i++ ){ img_data[i] += v; }
}


void ImgProc::operator-=(float v)
{
   if( img_data == nullptr ){ return; }
#pragma omp parallel for
   for( long i=0;i<Nsize;i++ ){ img_data[i] -= v; }
}


void ImgProc::compliment()
{
   if( img_data == nullptr ){ return; }

   for( long i=0;i<Nsize;i++ ){ img_data[i] = 1.0 - img_data[i]; }
}


void ImgProc::increaseBrightness(float brightnessNumber)
{
   if( img_data == nullptr ){ return; }

   //float brightnessNumber = 1.5;

   for( long i=0;i<Nsize;i++ ){ img_data[i] *= brightnessNumber; }
}


void ImgProc::decreaseBrightness(float brightnessNumber)
{
   if( img_data == nullptr ){ return; }

   //float brightnessNumber = 0.75;

   for( long i=0;i<Nsize;i++ ){ img_data[i] *= brightnessNumber; }
}


void ImgProc::increaseBias()
{
   if( img_data == nullptr ){ return; }

   float bias = 0.075; 

   for( long i=0;i<Nsize;i++ ){ img_data[i] += bias; }
}


void ImgProc::decreaseBias(float bias)
{
   if( img_data == nullptr ){ return; }

   //float bias = 0.075; 

   for( long i=0;i<Nsize;i++ ){ img_data[i] -= bias; }
}


void ImgProc::increaseGamma()
{
   if( img_data == nullptr ){ return; }

   float gamma = 1.5; 

   for( long i=0;i<Nsize;i++ ){ img_data[i] = pow(img_data[i], gamma); }
}


void ImgProc::decreaseGamma()
{
   if( img_data == nullptr ){ return; }

   float gamma = 0.75; 

   for( long i=0;i<Nsize;i++ ){ img_data[i] = pow(img_data[i], gamma); }
}


//The definition of grayscale dictates there be 3 channels: R, G, B.
void ImgProc::createGrayscale()
{
   if( img_data == nullptr ){ return; }

   float grayscale = 0;

   for( long i=0;i<Nsize;i++ )
   { 
      if (i % 3 == 0){grayscale = (img_data[i] * 0.2126);} //red
      
      else if (i % 3 == 1){grayscale += (img_data[i] * 0.7152);} //green

      else 
      {
         grayscale += (img_data[i] * 0.0722); //blue
         
         img_data[i] = grayscale;
         img_data[i - 1] = grayscale;
         img_data[i - 2] = grayscale;

         grayscale = 0;
      }
   }
}


void ImgProc::quantize()
{
   if( img_data == nullptr ){ return; }

   float numberOfColors = 5.0; 

   int intPart = 0;

   float intensityOfChannel = 0;

   for( long i=0;i<Nsize;i++ )
   { 
      intensityOfChannel = img_data[i];

      intPart = (int)(intensityOfChannel * numberOfColors);

      img_data[i] = ((float)(intPart))/numberOfColors;	
   }
}


void ImgProc::rmsContrast()
{
   if( img_data == nullptr ){ return; }

   double channelTotal [Nc] = {0};
   double channelMean [Nc] = {0};
   double channelStdDevSquared [Nc] = {0};
   double channelStdDev [Nc] = {0};

   double numOfPixels = (double)(Nx * Ny);
   
   long cntr = 0; //this counter will be used to move back to the first element of an array
                  //each time the channels are cycled through 
   
   for( long i=0;i<Nsize;i++ )
   { 
      channelTotal [cntr] += img_data[i];

      if (cntr == Nc - 1) {cntr = 0;}
      else {cntr ++;}
   }

   for( long i=0;i<Nc;i++ )
   {
      channelMean [i] = channelTotal [i] / numOfPixels; 
   }

   cntr = 0;
   
   //I am using this STD DEV formula that requires an extra loop because Dr. Tessendorf said it is safer
   for( long i=0;i<Nsize;i++ ) 
   { 
      channelStdDevSquared [cntr] += pow(img_data[i] - channelMean [cntr], 2);

      if (cntr == Nc - 1) {cntr = 0;}
      else {cntr ++;}
   }

   for( long i=0;i<Nc;i++ ) 
   {
      channelStdDevSquared [i] /= numOfPixels;

      channelStdDev [i] = sqrt(channelStdDevSquared[i]);
   }

   cntr = 0;

   for( long i=0;i<Nsize;i++ )
   { 
      img_data[i] = (img_data[i] - channelMean [cntr]) / channelStdDev [cntr];

      if (cntr == Nc - 1) {cntr = 0;}
      else {cntr ++;}
   }
}


void ImgProc::createImageFile()
{

   const char *filename = "newFile.exr";
   const int xres = Nx;
   const int yres = Ny;
   const int channels = Nc;

   ImageOutput *out = (ImageOutput::create (filename)).get();
   if (!out) {return;}

   ImageSpec spec (xres, yres, channels, TypeDesc::FLOAT);
   out -> open (filename, spec);
   out -> write_image (TypeDesc::FLOAT, img_data);
   out -> close();
   ImageOutput::destroy (out);
}


void ImgProc::realSpaceCoherentEstimation(int firstChannel, int secondChannel, float firstWeight, float secondWeight)
{
   if ( img_data == nullptr ){ return; }

   double firstChannelAverage = 0;
   double secondChannelAverage = 0;
   double firstChannelSigma = 0;
   double secondChannelSigma = 0;
   double crossChannelSigma = 0;
   float coherentEstimation = 0;


   for (int j = 0; j < Ny; j++)
   {
      for(int i = 0; i < Nx; i++)
      {
         std::vector<float> pixel;

         value(i, j, pixel);

         firstChannelAverage += (pixel[firstChannel]);
         secondChannelAverage += (pixel[secondChannel]);
      }
   }
   
   firstChannelAverage /= (Nx * Ny);
   secondChannelAverage /= (Nx * Ny);

   for (int j = 0; j < Ny; j++)
   {
      for(int i = 0; i < Nx; i++)
      {
         std::vector<float> pixel;

         value(i, j, pixel);

         firstChannelSigma += (pixel[firstChannel] - firstChannelAverage) * (pixel[firstChannel] - firstChannelAverage);
         secondChannelSigma += (pixel[secondChannel] - secondChannelAverage) * (pixel[secondChannel] - secondChannelAverage);
         crossChannelSigma += (pixel[firstChannel] - firstChannelAverage) * (pixel[secondChannel] - secondChannelAverage);

      }
   }
 
   firstChannelSigma /= (Nx * Ny);
   secondChannelSigma /= (Nx * Ny);
   crossChannelSigma /= (Nx * Ny);

   float determinant = firstChannelSigma * secondChannelSigma - crossChannelSigma * crossChannelSigma;

   for (int j = 0; j < Ny; j++)
   {
#pragma omp parallel for
      for(int i = 0; i < Nx; i++)
      { 
         std::vector<float> pixel;

         value(i, j, pixel);

         coherentEstimation = firstWeight * (pixel[firstChannel] - firstChannelAverage) * secondChannelSigma
                            - firstWeight * (pixel[secondChannel] - secondChannelAverage) * crossChannelSigma
                            - secondWeight * (pixel[firstChannel] - firstChannelAverage) * crossChannelSigma
                            + secondWeight * (pixel[secondChannel] - secondChannelAverage) * firstChannelSigma;

         float val = coherentEstimation / determinant;

         //Fill each channel in pixel with the new value. Hard coded to work inside the #pragma
         pixel[0] = val;
         pixel[1] = val;
         pixel[2] = val;

         set_value(i, j, pixel);
      }
   }
   decreaseBias(3.75);
}


long ImgProc::index(int i, int j, int c) const
{
   return (long) c + (long) Nc * index(i,j); // interleaved channels

   // return index(i,j) + (long)Nx * (long)Ny * (long)c; // sequential channels
}


long ImgProc::index(int i, int j) const
{
   return (long) i + (long)Nx * (long)j;
}


void img::swap(ImgProc& u, ImgProc& v)
{
   float* temp = v.img_data;
   int Nx = v.Nx;
   int Ny = v.Ny;
   int Nc = v.Nc;
   long Nsize = v.Nsize;

   v.Nx = u.Nx;
   v.Ny = u.Ny;
   v.Nc = u.Nc;
   v.Nsize = u.Nsize;
   v.img_data = u.img_data;

   u.Nx = Nx;
   u.Ny = Ny;
   u.Nc = Nc;
   u.Nsize = Nsize;
   u.img_data = temp;
}
