#ifndef RANGE_H_
#undef max
#undef min
struct Range{
  int min;
  int max;

 Range(int _min = -1, int _max = -2):min(_min), max(_max){}
};


#define RANGE_H_
#endif
