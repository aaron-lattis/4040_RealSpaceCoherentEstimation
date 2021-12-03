

#ifndef IMGPROC_H
#define IMGPROC_H

#include <string>
#include <vector>

namespace img
{


class ImgProc
{

  public:

    //! Construct with no content
    ImgProc();
   ~ImgProc();
   
    //! delete existing content and leave in a blank state
    void clear();
    //! delete existing content and re-initialize to the input dimensions with value 0.0
    void clear(int nX, int nY, int nC);

    //! Load an image from a file.  Deletes exising content.
    bool load( const std::string& filename );

    //! Retrieve the width
    int nx() const { return Nx; }
    //! Retrieve the height
    int ny() const { return Ny; }
    //! Retrieve the number of channels
    int depth() const { return Nc; }

    //! Retrieve the (multichannel) value at a pixel.  Copies the value into parameter 'pixel'.
    void value( int i, int j, std::vector<float>& pixel) const;
    //! Set the (multichannel) value at a pixel.
    void set_value( int i, int j, const std::vector<float>& pixel);

    //! Copy constructor. Clears existing content.
    ImgProc(const ImgProc& v);
    //! Copy assignment. Clears existing content.
    ImgProc& operator=(const ImgProc& v);

    friend void swap(ImgProc& u, ImgProc& v);

    //! multiplies all pixels and channels by a value
    void operator*=(float v);
    //! divides all pixels and channels by a value
    void operator/=(float v);
    //! adds a value to all pixels and channels
    void operator+=(float v);
    //! subtracts a value from all pixels and channels
    void operator-=(float v);

    //! converts image to its compliment in-place
    void compliment();

    //! increases the brightness of the image in-place
    void increaseBrightness(float brightnessNumber);

    //! decreases the brightness of the image in-place
    void decreaseBrightness(float brightnessNumber);

    //! increases the bias of the image in-place
    void increaseBias();

    //! decreases the bias of the image in-place
    void decreaseBias(float bias);

    //! increases the gamma of the image in-place
    void increaseGamma();

    //! decreases the gamma of the image in-place
    void decreaseGamma();

    //! applies gray scale to the image in-place
    void createGrayscale();

    //! quantizes the image in-place
    void quantize();

    //! performs an RMS contrast on the image in-place
    void rmsContrast();

    //! writes the image to a new .EXR file
    void createImageFile();

    //! indexing to a particular pixel and channel
    long index(int i, int j, int c) const;
   
    //! indexing to a particular pixel
    long index(int i, int j) const;

    //! returns raw pointer to data (dangerous)
    float* raw(){ return img_data; }

    //! Performs real space coherent estimation on the image
    void realSpaceCoherentEstimation(int firstChanel, int secondchannel, float firstWeight, float secondWeight); 



  private:

    int Nx, Ny, Nc;
    long Nsize;
    float * img_data;
};


//! swaps content of two images
void swap(ImgProc& u, ImgProc& v);


}
#endif
