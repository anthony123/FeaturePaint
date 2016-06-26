#ifndef COLOR_H_
#define COLOR_H_
class Color{
  int  red, blue,green;
 public:
 Color(int r = 0, int g = 0, int b = 0):red(r),blue(b), green(g){}

  int operator()(int index) const;
  void setValues(int red, int green, int blue);

  bool operator== (const Color& color);
  bool operator!= (const Color& color);
  
};

#endif
