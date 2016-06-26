#ifndef MASKED_IMAGE_H_
#define MASKED_IMAGE_H_
#include "Point.h"
#include "CImg.h"

using namespace cimg_library;


template<typename T> class MaskedImage:public CImg<T>{
 private:
  CImg<bool>mask;


 public:
 MaskedImage<T>(const CImg<T>&img, const CImg<bool>& img_mask):CImg<T>(img), mask(img_mask){}

  MaskedImage<T>():CImg<T>(){mask = CImg<bool>();}

  MaskedImage<T>(const MaskedImage<T>& copy):CImg<T>(copy), mask(copy.getMask()){}
  //int width() const {return image.width();}
  //int height() const {return image.height();}
  void halfSize();
  void doubleSize();
  //CImg<int> getImage() const;
  CImg<bool> getMask() const;
  bool isMasked(Point pt) const;
  //int operator()(int x, int y, int band) const;
  //void setValue(int x, int y, int band, int value);
  void setMask(int x, int y, bool masked);
  void setMaskImage(const CImg<bool>&img);
  MaskedImage<T>& operator=(const MaskedImage<T>& img);
};


#endif
