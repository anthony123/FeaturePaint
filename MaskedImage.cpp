#include "MaskedImage.h"
#include "CImg.h"
#include "Point.h"
#include <iostream>

using namespace cimg_library;

template<typename T> void MaskedImage<T>::halfSize()
{
    CImg<T>::resize_halfXY();
  int newH = mask.height()/2;
  int newW = mask.width()/2;

  CImg<bool>new_mask(newW,newH, 1,1,true);
  cimg_forXY(mask, x,y){
    if(mask(x,y,0) == false){
      int newX = x/2;
      int newY = y/2;
      new_mask(newX, newY, 0) = false;
    }
  }

  mask = new_mask;
}

template<typename T> void MaskedImage<T>::doubleSize()
{
  CImg<T>::resize_doubleXY();
  int newH = mask.height()*2;
  int newW = mask.width()*2;
  CImg<bool>new_mask(newW, newH, 1,1,true);
  cimg_forXY(new_mask, x, y)
  {
    int oldX = x/2;
    int oldY = y/2;

    if(mask(oldX,oldY, 0) == false){
      new_mask(x,y,0) = false;
    }
  }

  mask = new_mask;
}

template<typename T> CImg<bool> MaskedImage<T>::getMask() const
{
  return mask;
}


template<typename T> bool MaskedImage<T>::isMasked(Point pt) const
{
  return mask(pt.x, pt.y, 0);
}


template<typename T> void MaskedImage<T>::setMask(int x, int y, bool masked)
{
  mask(x,y,0) = masked;
}
template<typename T> void MaskedImage<T>::setMaskImage(const CImg<bool>&img){
  mask = img;
}
template<typename T> MaskedImage<T>& MaskedImage<T>::operator=(const MaskedImage& img)
{

  CImg<T>::operator=(img);
  mask  = img.getMask();

  return *this;
}
