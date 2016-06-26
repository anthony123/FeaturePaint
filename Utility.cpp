#include <vector>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <armadillo>
#include "Point.h"
#include "Utility.h"
#include "typeRigid.h"
#include "FeaturePaint.h"
#include "Range.h"

using namespace arma;

//find four edge points centered in a point

bool Utility::compx(const Point& p1, const Point& p2)
{
  return p1.x < p2.x;
}

bool Utility::compy(const Point& p1, const Point& p2)
{
  return p1.y < p2.y;
}

float Utility::computeDistance(const Point& p1, const Point& p2)
{
  return sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y));
}

Point Utility::findNormal(int index, const std::vector<Point>&pts)
{
  int sz = pts.size();
  Point normal, p0, p1, p2;
  float dx, dy, mag;

  if(index == 0){
    p0 = pts[0];
    p1 = pts[1];

    dx = p1.y - p0.y;
    dy = p0.x - p1.x;

    mag = sqrt(dx*dx + dy*dy);

    dx /= mag;
    dy /= mag;
  }else if(index == sz-1){
    p0 = pts[sz-2];
    p1 = pts[sz-1];

    dx = p1.y - p0.y;
    dy = p0.x - p1.x;

    mag = sqrt(dx*dx + dy*dy);
    dx /= mag;
    dy /= mag;
  }else{
    p0 = pts[index-1];
    p1 = pts[index];
    p2 = pts[index+1];

    dx = p2.y - p0.y;
    dy = p0.x - p2.x;

    mag = sqrt(dx*dx + dy*dy);

    dx /= mag;
    dy /= mag;
  }

  return Point(dx, dy);
}

Point Utility::findTangent(int index, const std::vector<Point>&pts)
{
  Point normal = findNormal(index, pts);
  return Point(normal.y, -normal.x);
}

std::vector<Point>Utility::findEdgePoints(int index, const std::vector<Point>&pts)
{
  float height_radius = PATCH_HEIGHT/2.0f;
  float width_radius = PATCH_WIDTH/2.0f;
  //float radius = 1.5f;
  Point v0,v1, p0,p1,p2, pt1, pt2, pt3, pt4, tangent, normal;

  // if(index == 0){
  //   p0 = pts[0];
  //   p1 = pts[1];

  //   dx = p1.y - p0.y;
  //   dy = p0.x - p1.x;

  //   mag = sqrt(dx*dx + dy*dy);

  //   dx /= mag;
  //   dy /= mag;

  //   v0 = Point(p0.x + height_radius * dx, p0.y + height_radius * dy);
  //   v1 = Point(p0.x - height_radius * dx, p0.y - height_radius * dy);


  // }else if(index == sz-1){
  //   p0 = pts[sz-2];
  //   p1 = pts[sz-1];

  //   dx = p1.y - p0.y;
  //   dy = p0.x - p1.x;

  //   mag = sqrt(dx*dx + dy*dy);
  //   dx /= mag;
  //   dy /= mag;

  //   v0 = Point(p1.x + height_radius * dx, p1.y + height_radius * dy);
  //   v1 = Point(p1.x - height_radius * dx, p1.y - height_radius * dy);



  // }else{
  //   p0 = pts[index-1];
  //   p1 = pts[index];
  //   p2 = pts[index+1];

  //   dx = p2.y - p0.y;
  //   dy = p0.x - p2.x;

  //   mag = sqrt(dx*dx + dy*dy);

  //   dx /= mag;
  //   dy /= mag;

  //   v0 = Point(p1.x + height_radius * dx, p1.y + height_radius * dy);
  //   v1 = Point(p1.x - height_radius * dx, p1.y - height_radius * dy);
  // }

  //tangent = Point(dy, -dx);
  //normal = Point(dx, dy);
  tangent = findTangent(index, pts);
  normal = findNormal(index, pts);
  Point p = pts[index];
  v0 = p + normal*height_radius;
  v1 = p - normal*height_radius;


  //find four edge points

  pt1 = v1 - tangent*width_radius;
  pt2 = v1 + tangent*width_radius;
  pt3 = v0 - tangent*width_radius;
  pt4 = v0 + tangent*width_radius;

  std::vector<Point>res;
  res.push_back(pt1);
  res.push_back(pt2);
  if(computeDistance(pt2,pt3) - computeDistance(pt2,pt4) >= 0){
    res.push_back(pt4);
    res.push_back(pt3);
  }
  else{
    res.push_back(pt3);
    res.push_back(pt4);
  }

  //must arrange the sequence
  return res;
}

//determine whether the point is in the rectangle centered in pts[index];
bool Utility::pointInRectangle(const Point& p, int index, const std::vector<Point>&pts)
{
  Point  pt1, pt2, pt3, pt4;
  std::vector<Point>edges = findEdgePoints(index, pts);
  pt1 = edges[0];
  pt2 = edges[1];
  pt3 = edges[2];
  pt4 = edges[3];

  //cout<<"the target point is ["<<p.x<<", "<<p.y<<"]\n";
  //std::cout<<"the center point is ["<<pts[index].x<<", "<<pts[index].y<<"]\n";
  //std::cout<<"edge points are ["<<pt1.x<<","<<pt1.y<<"], ["<<pt2.x<<","<<pt2.y<<"], ["<<pt3.x<<","<<pt3.y<<"], ["<<pt4.x<<","<<pt4.y<<"]\n";
  //to determine if the point in the rectangle, just do to-left test four times
  bool p12 = toLeft(pt1, pt2, p);
  bool p23 = toLeft(pt2, pt3, p);
  bool p34 = toLeft(pt3, pt4, p);
  bool p41 = toLeft(pt4, pt1, p);

  return (p12 == p23)&&(p23 == p34)&&(p34==p41);

}

inline bool Utility::toLeft(const Point& p1, const Point& p2, const Point& q)
{
  return Utility::area2(p1,p2,q) > 0;
}

inline int Utility::area2(const Point& p1,  const Point& p2,  const Point& p3)
{
  //return p1.x*p2.y + p2.x*q.y + p1.y*q.x - p2.y*q.x - p1.y*p2.x - p1.x*q.y;
  return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
}

inline  bool Utility::isBetween(const Point& p1, const Point& p2, const Point& q){
  return Utility::area2(p1, p2, q) == 0 && (q.x-p1.x)*(p2.x-p1.x) + (q.y - p1.y)*(p2.y - p1.y) < (p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y);
}

bool Utility::isExist(const vector<Point>&pts, Point pt, const vector<Range>&range)
{
  int min_x = pts.begin()->x;
  int max_x = pts.rbegin()->x;
  int x = pt.x;
  int y = pt.y;
  int delta = min_x;
  if(x == 126){
    cout<<"miny = "<<range[x-delta].min<<", maxy = "<<range[x-delta].max<<endl;
  }
  if( x >= min_x && x <= max_x && y >= range[x-delta].min && y <= range[x-delta].max){
    return true;
  }else{
    return false;
  }

}

int Utility::findPoint(const vector<Point>&pts, Point pt, const vector<Range>&range)
{
  int x = pt.x;
  int y = pt.y;
  int delta = pts.begin()->x;
  int index = 0;
  int offset = x - delta;
  if(isExist(pts, pt, range)){
    for(int i = 0; i < offset; i++){
      index = index +  (range[i].max-range[i].min+1);
      //cout<<"x = "<<i+delta<<", minX = "<<range[i].min<<", maxX = "<<range[i].max<<endl;
      //cout<<"the current point is ("<<pts[index].x<<", "<<pts[index].y<<")"<<endl;
      assert(range[i].max == pts[index-1].y);
    }

    assert(pts[index].x == x && pts[index-1].x == (x-1));
    //cout<<"point in points: "<<"("<<pts[index].x<<", "<<pts[index].y<<")"<<", the target point: ("<<x<<", "<<y<<")\n";
    while(index < pts.size()){
      if(pts[index].y == y)
	break;
      else
	++index;

    }
    assert(pts[index].x == x && pts[index].y == y);
    return index;
  }else{
    return -1;
  }
}

/*unordered_set<Point, Utility::hash_func> Utility:: getEffectivePoints(unordered_set<Point, Utility::hash_func>&points, int size){
  unordered_set<Point, Utility::hash_func>res;
  for(const auto &elem:points){
  int x = elem.x;
  int y = elem.y;

  int exist = true;
  for(int i = 1; i < size; i++){
  for(int j = 1; j < size; j++){
  Point candidate(x+i, y+j);
  if(points.find(candidate) == points.end()){
  exist = false;
  }
  }
  }

  if(exist)
  res.insert(elem);
  }

  return res;
  }*/

double * Utility::generateD2S()
{
  double base[] = {1.0, 0.99, 0.96, 0.83, 0.38, 0.11, 0.02, 0.005, 0.0006, 0.0001, 0 };
  //double base[] = {1.0, 0.995, 0.99, 0.975, 0.96, 0.895,0.83, 0.605,0.38, 0.245,0.11,0.065, 0.02,0.0125, 0.005,0.0028, 0.0006, 0.00035,0.0001,0.00005, 0 };
  //double base[] = {1.0, 0.9975, 0.995, 0.9925, 0.99,0.9825, 0.975, 0.9675,0.96,0.9275, 0.895, 0.8625,0.83,0.7175, 0.605, 0.4925, 0.38,0.3125, 0.245, 0.1775,0.11, 0.0875,0.065, 0.0425,0.02, 0.01625, 0.0125,0.00875, 0.005, 0.0039, 0.0028, 0.0017, 0.0006, 0.000475, 0.00035,0.000255,0.0001,0.000075,0.00005,0.000025, 0 };
  int base_length = sizeof(base)/sizeof(*base);

  int length = DSCALE+1;
  double * similarity = (double*)calloc(length, sizeof(double));
  for ( int i=0 ; i<length ; ++i) {
    double t = (double)i/length;
    int j = (int)(100*t);
    int k=j+1;
    double vj = (j<base_length)?base[j]:0;
    double vk = (k<base_length)?base[k]:0;
    similarity[i] = vj + (100*t-j)*(vk-vj);
    //cout<<"similarity["<<i<<"] = "<<similarity[i]<<endl;
  }

  return similarity;
}

double Utility::distanceToSimilarity(double* similarity, int distance)
{
  //cout<<"int the Utility function, distance = "<<distance<<endl;
  if(distance > DSCALE)
    distance = DSCALE;
  double ret = similarity[distance];
  //cout<<"ret = "<<ret<<endl;
  //free(similarity);
  return ret;
}

mat Utility::precomputeWeights(mat p, mat v, double d)
{
  mat w;
  w.zeros(p.n_cols, v.n_cols);
  mat p_resize;
  mat norms;
  //norms.zeros(2, v.n_cols);
  mat norms_a;
  mat p_v;

  //Iterate through the control points
  for(int i = 0 ; i<p.n_cols; i++)
    {
      //compute the norms
      p_resize = repmat(p.col(i),1, v.n_cols);
      p_v = p_resize -v;
      p_v = pow(p_v,2);
      norms = p_v.row(0) + p_v.row(1);
      norms_a = pow(norms,d);
        
      //compute the weights
      w.row(i) = 1.0 / norms_a;
    }
  return w;
}

mat Utility::precomputeWCentroids(mat p, mat w)
{
  mat Pstar;
  mat mult;
  mat resize;
  mat sum = zeros(1, w.n_cols);
  mult= p*w;

  for(int i = 0; i< w.n_rows; i++)
    sum += w.row(i);
  resize = repmat(sum, p.n_rows, 1);

  Pstar = mult / resize;
  return Pstar;
}

vector <_typeA> Utility::precomputeA(mat Pstar, vector <mat> Phat, mat v, mat w)
{
  vector <_typeA> A;

  //fixed part
  mat R1 = v - Pstar;
  mat R2;
  //vconcat(R1.row(1),-R1.row(0),R2);
  //R2 = join_cols(R1.row(1), -1*R1.row(0));
    R2 = -1*join_cols(R1.row(1), -1*R1.row(0));
  //R2 = join_rows(R1.row(1), -R1.row(0));

  for(int i = 0 ; i< Phat.size() ; i++)
    {
      //precompute
      typeA temp;
      mat L1 = Phat.at(i);
      mat L2;
      //vconcat(L1.row(1),(L1.row(0)).mul(-1),L2);
      //L2 = join_cols(L1.row(1), -1*(L1.row(0)));
        L2 = -1*join_cols(L1.row(1), -1*(L1.row(0)));
      //L2 = join_rows(L1.row(1), -1*(L1.row(0)));

      // cout<<"L1: "<<"("<<L1.n_cols<<","<<L1.n_rows<<")"<<endl;
      // cout<<"R1: "<<"("<<R1.n_cols<<","<<R1.n_rows<<")"<<endl;
      // cout<<"L2: "<<"("<<L2.n_cols<<","<<L2.n_rows<<")"<<endl;
      // cout<<"R2: "<<"("<<R2.n_cols<<","<<R2.n_rows<<")"<<endl;
      
      //mat L1R1 = L1*R1;
      
      mat L1R1 = L1%(R1);
      mat sumL1R1 = zeros(1,L1R1.n_cols);

      //mat L1R2 = L1*R2;
      mat L1R2 = L1%(R2);
      mat sumL1R2 = zeros(1, L1R2.n_cols);

      //mat L2R1 = L2*R1;
      mat L2R1 = L2%(R1);
      mat sumL2R1 = zeros(1, L2R1.n_cols);

      //mat L2R2 = L2*R2;
      mat L2R2 = L2%(R2);
      mat sumL2R2 = zeros(1, L2R2.n_cols);

      for(int j = 0; j<L1R1.n_rows; j++)
	sumL1R1 += L1R1.row(j);

      for(int j = 0; j<L1R2.n_rows; j++)
	sumL1R2 += L1R2.row(j);

      for(int j = 0; j<L2R1.n_rows; j++)
	sumL2R1 += L2R1.row(j);

      for(int j = 0; j<L2R2.n_rows; j++)
	sumL2R2 += L2R2.row(j);

      // cout<<"sumL1R1: "<<"("<<sumL1R1.n_cols<<","<<sumL1R1.n_rows<<")"<<endl;
      

      temp.a = (w.row(i))%(sumL1R1);
      temp.b = (w.row(i))%(sumL1R2);
      temp.c = (w.row(i))%(sumL2R1);
      temp.d = (w.row(i))%(sumL2R2);

      A.push_back(temp);
    }

  
  return A;
}

typeRigid Utility::precomputeRigid(mat p, mat v, mat w)
{

  typeRigid data;
  mat Pstar = precomputeWCentroids(p, w);
  cout<<"pre compute W centroids done"<<endl;
  vector <mat> Phat;
  for(int i =0 ; i<p.n_cols; i++)
    {
      mat t = repmat(p.col(i) ,1, Pstar.n_cols ) - Pstar;
      Phat.push_back(t);
    }

  vector <_typeA> A = precomputeA(Pstar, Phat, v, w);
  cout<<"pre compute A done"<<endl;
  mat v_Pstar = v-Pstar;
  mat vpower;
  vpower = pow(v_Pstar, 2);
  mat sum = zeros(1,vpower.n_cols);
  for(int i =0; i<vpower.n_rows; i++)
    sum += vpower.row(i);

  data.normof_v_Pstar = arma::sqrt(sum);
  data.A = A;
  return data;
}

mat Utility::PointsTransformRigid(mat w,typeRigid mlsd, mat q)
{
  mat Qstar = precomputeWCentroids(q, w);
  mat Qhat;
  mat fv2 = zeros(Qstar.n_rows, Qstar.n_cols);
  mat prod1, prod2;
  mat con1, con2;
  mat update;
  mat repeatmat;
  mat npower;
  mat normof_fv2;
  mat fv = zeros(Qstar.n_rows, Qstar.n_cols);
  for( int i = 0 ; i< q.n_cols; i++)
    {
      Qhat = repmat(q.col(i) ,1, Qstar.n_cols ) - Qstar;

      //vconcat((mlsd.A.at(i)).a, (mlsd.A.at(i)).c, con1);
      con1 = join_cols((mlsd.A.at(i)).a, (mlsd.A.at(i)).c);
      prod1 = Qhat%(con1);
      mat sum1 = zeros(1, prod1.n_cols);
      for(int j = 0 ; j<prod1.n_rows; j++)
	sum1+=prod1.row(j);

      //vconcat((mlsd.A.at(i)).b, (mlsd.A.at(i)).d, con2);
      con2 = join_cols((mlsd.A.at(i)).b, (mlsd.A.at(i)).d);
      prod2 = Qhat%(con2);
      mat sum2 = zeros(1, prod2.n_cols);
      for(int j = 0 ; j<prod2.n_rows; j++)
	sum2+=prod2.row(j);

      //vconcat(sum1, sum2, update);
      update = join_cols(sum1, sum2);
      fv2 = fv2 + update;
    }
  npower = fv2%(fv2);

  mat sumfv2 = zeros(1,npower.n_cols);
  for(int i =0; i<npower.n_rows; i++)
    sumfv2 += npower.row(i);


  normof_fv2 = 	sqrt(sumfv2);

  mat norm_fact = (mlsd.normof_v_Pstar)%(1/normof_fv2);

  repeatmat = repmat(norm_fact, fv2.n_rows, 1);
  fv = fv2%repeatmat + Qstar;


  return fv;
}
