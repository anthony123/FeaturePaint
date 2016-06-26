#include "Color.h"

int Color::operator()(int index)const{
  switch(index){
  case 0: return red; break;
  case 1: return green; break;
  case 2: return blue; break;
  default: return -1; break;
  }
}

void Color::setValues(int r, int g, int b){
  red = r;
  green = g;
  blue = b;
}

bool Color::operator==(const Color& color){
  return (red == color.red)&&(green == color.green)&&(blue == color.blue);
}

bool Color::operator!=(const Color&color){
  return (red != color.red)||(green != color.green) ||(blue != color.blue);
    
}
