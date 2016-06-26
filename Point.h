#ifndef POINT_H_
#define POINT_H_
struct Point{
  float x;
  float y;

  //constructors
Point(float xx = -1, float yy = -1):x(xx),y(yy){}

Point(const Point &point):x(point.x),y(point.y){}

  //move vertically, make sure return a valid point
  Point verticallyMove(Point point, int pace);

  //overload =
  Point& operator=(const Point& point);

  //overload !=
  bool operator!=(const Point& point)const;

  //overload ==
  bool operator== (const Point& point) const;
    
    //overload +
    Point operator+(const Point & point)const ;
    
    //overload -
    Point operator-(const Point& point)const;
    
    //overload +=
    Point & operator+= (const Point &point);
    
    //overload *
    Point operator*(float v)const;
    
    //overload /
    Point operator/(float v)const;

  //overload ()
  float operator()(int dimension);

  //overload <
  bool operator<(const Point& point)const;

  //for set
  bool operator<(Point& point);

};



#endif















