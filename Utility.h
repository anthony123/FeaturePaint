#ifndef UTILITY_H_
#define UTILITY_H_
#include <cstdlib>
#include <vector>
#include <armadillo>
#include <cmath>
#include "Point.h"
#include "typeRigid.h"
#include "Range.h"

using namespace arma;


namespace Utility {
  bool compx(const Point& p1, const Point& p2);
  bool compy(const Point& p1, const Point& p2);
  float computeDistance(const Point& p1, const Point& p2);
  bool pointInRectangle(const Point& p, int index, const std::vector<Point>&pts);
  bool toLeft(const Point &p1, const Point &p2, const Point &q);
  std::vector<Point>findEdgePoints(int index, const std::vector<Point>&pts);
  Point findTangent(int index, const std::vector<Point>&pts);
  Point findNormal(int index, const std::vector<Point>&pts);
  int area2(const Point& p1, const Point& p2, const Point& p3);

  //determine if a point is between two other points on a line segment
  bool isBetween(const Point& p1, const Point& p2, const Point& q);

  //find the index of a specific point in a sequence of points
  int findPoint(const std::vector<Point>&pts, Point pt, const std::vector<Range>&range);

  //determine if the point is exist in the points set
  bool isExist(const std::vector<Point>&pts, Point pt, const std::vector<Range>&range);

  double* generateD2S();

  /* struct hash_func{ */
  /*   int operator()(const Point& p) const */
  /*   { */
  /*       return pow(2, p.x) + pow(3, p.y); */
  /*   } */
  /* }; */

  /*struct cmp_func{
    bool operator()( const Point& p1,  const Point& p2)const
    {
        return p1.x == p2.x && p1.y == p2.y;
    }
  };

  struct hash_func
  {
    int operator()(Point const & p) const
    {
      return (
	      (51 + std::hash<int>()(p.x)) * 51
	      + std::hash<int>()(p.y)
	      );
    }
  };

  std::unordered_set<Point, Utility::hash_func> getEffectivePoints(std::unordered_set<Point, Utility::hash_func>&points, int size);*/

  double distanceToSimilarity(double* similarity, int distance);

  //related to MLS
  mat precomputeWeights(mat p, mat v, double d);
  mat precomputeWCentroids(mat p, mat w);
    vector <_typeA> precomputeA(mat Pstar, vector <mat> Phat, mat v, mat w);
  typeRigid precomputeRigid(mat p, mat v, mat w);
  mat PointsTransformRigid(mat w,typeRigid mlsd, mat q);
}

#endif
