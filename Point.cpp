#include "Point.h"

//move vertically
/*Point Point::verticallyMove(Point point, int pace){
  Point res;
  int i,j;
  int diff_x,diff_y;
  diff_x = point.x-x;
  diff_y = point.y-y;
    
  
  //search from neighbor 8 points
  for(i = -1; i < 2; i++){
    for(j = -1; j < 2; j++){
      //vertical
      if((diff_x*i + diff_y+j) == 0)
	break;
    }
  }

  int ver_x = pace*i;
  int ver_y = pace*j;

  return Point(x+ ver_x, y + ver_y);
  }*/

Point& Point::operator=(const Point&point){
  x = point.x;
  y = point.y;

  return *this;
}

bool Point::operator!=(const Point& point)const{
  return (x != point.x || y != point.y);
}

bool Point::operator==(const Point&point)const{
  return (x == point.x && y == point.y );
}

Point Point::operator+(const Point& point)const{
  return Point(x+point.x, y+point.y);
}

Point Point::operator-(const Point& point)const{
  return Point(x-point.x, y-point.y);
}

Point Point::operator*(float v)const{
  return Point(x*v, y*v);
}

Point Point::operator/(float v)const{
  return Point(x/v, y/v);
}

Point& Point::operator+=(const Point& point){
    x += point.x;
    y += point.y;
    
    return *this;
}

bool Point::operator<(const Point& point)const{
  return (x < point.x)||(x == point.x && y < point.y);
}

float Point::operator()(int dimension){
  switch(dimension){
  case 0: return x;
  case 1: return y;
  default: return -1.0f;
  }


}

bool Point::operator<(Point& point){
  return x < point.x;
}
















