#ifndef FEATURE_PAINT_H
#define FEATURE_PAINT_H
#include <vector>
#include <stack>
#include <set>
#include "CImg.h"
#include "Color.h"
#include "Range.h"
#include "Point.h"
#include "Utility.h"
#include "NNF.h"
#include "MaskedImage.h"

using namespace cimg_library;
using namespace std;

const Color RED(255,0,0);
const Color GREEN(0,255,0);
const Color BLUE(0, 0, 255);
const Color BLACK(0,0,0);
const Color WHITE(255,255,255);
const int DSCALE = 65535;



//stroke width
const int PATCH_WIDTH = 7;
//const int PATCH_WIDTH = 10;

//sample pixel width
const int PATCH_HEIGHT = 15;
//const int PATCH_HEIGHT = 10;

//threshold of gradient
const int THRESHOLD = 30;

//the range of finding the corresponding point on the contour
const int RANGE = PATCH_HEIGHT/2;

//width and height of fill patch
const int FILL_PATCH = 7;

//const int MIN_PIXELS = PATCH_WIDTH*PATCH_WIDTH;
//const int MIN_PIXELS = PATCH_WIDTH;
const int MIN_PIXELS = PATCH_WIDTH;

//number of iterators for patch match
//const int ITERS = 15;

//parameter r
//const float PARAM_R = 1.3;

class FeaturePaint{
  //the source image as feature palette
  CImg<int> source_image;

  //the stroke we select
  CImg<int> stroke_image;

  //the area we select
  CImg<int>area_image;

  //the contour image
  CImg<int> contour_image;

  //the output image
  CImg<int>target_image;


  //stroke points
  vector< CImgList<int> >strokes;

  //source stoke points
  vector<Point>source_stroke_points;

  //target stroke points
  vector<Point>target_stroke_points;

  double *simi;

  //number of stroke
  int nStroke;

  //number of area
  int nArea;

  int TARGET_IMAGE_WIDTH;
  int TARGET_IMAGE_HEIGHT;

  /****************private functions*************************/

  /***************************brush******************************/

  //find how many neighbors in the same color
  int numberOfNeighbors(CImg<> img, Point point, const Color& color);

  //find next point in the same color
  Point nextPoint(CImg<> img, CImgList<int>points, Point curr, Color color);

  //get user's stroke
  CImgList<int>getStroke(CImg<>img, const Color& color, bool isClosed);

  //warping the stroke to a  straight line: operate on the input_image but can not alter the input_image  and the stroke


  //sample the stroke
  CImgList<int>sampleStroke(CImgList<int>stroke);

  //adjust user's input-stroke to contour
  CImgList<int> adhereToContour(CImgList<int>stroke);



  //get the stroke of contour_image


  //rendering
  //void render(Point src, Point tar);
  void render(const vector<int>&path,CImgList<int>sample_stroke,CImgList<int>sample_contour_stroke);

  //distance of two points
  //int dist(const unordered_set<Point, Utility::hash_func>& target_pts, int sourceX, int sourceY,const unordered_set<Point, Utility::hash_func>& source_pts, int targetX, int targetY, CImg<int>&source_image, CImg<int>&target_image);


  /*****************************FILL******************************/
  //unordered_set<Point, Utility::hash_func>readArea(CImg<int>&img, const Color& color, bool isFilledArea = false);
  vector<Point>readArea(CImg<int>&img, const Color& color, bool isFilledArea = false);


 public:
   //colors of stroke/area
  vector<Color>stroke_colors;
  vector<Color>area_colors;

  //target color point
  Point target_color_point;

  //constructor
  FeaturePaint(char* source, char* stroke, char* area, char* contour);

  //set number of stroke
  void setNumberOfStroke(int n);

  //set number of area
  void setNumberOfArea(int n);

  //process one line: aggragate all steps into one function
  void processALine(const Color& color, bool isColosed);

  //Moving Least Square Warping
  void MLSWarping(const vector<Point>& source_control_points, const vector<Point>&target_control_points);

  void afterProcess();
  void multiScaleFill(const Color& target_color);

  //find the patch of a specified point
  int findPatch(Point p, const vector<Point>&pts);

  //find all stroke points
  void findAllStrokePoints(const vector<Point>&pts, const CImg<int>& image);

  //find corresponding pixel of source patch for the specified pixel in the target patch
  Point findCorrespondingPixel(Point p, int index_of_target_patch, int index_of_source_patch, const vector<Point>&tpts, const vector<Point>&spts);

  //mapping corresponding pixel value
  void mapPixel(Point target, Point source, CImg<int>&target_image, CImg<int>&source_image);
  void mapPatchPixel(Point target, Point source, CImg<int>&target_image, CImg<int>&source_image);

  CImg<float> gradientImage(CImg<int>image);

  //find target color position
  Point findTargetColorPoint(Color color, bool isFilledArea);

  bool testColorInTwoNeighbors(CImg<int>contour_image, Point point, Color color);

  MaskedImage<int>* ExpectationMaximization(vector<MaskedImage<int> >sources, vector<MaskedImage<int> >targets, NNF& S2Tnnf,NNF& T2Snnf,int level);
    void ExpectationStep(const NNF& nnf, int sourceToTarget, CImg<double>&vote,  const MaskedImage<int>& source, bool upscaled);

    void weightedCopy(const MaskedImage<int>& src, int xs, int ys, CImg<double>& vote, int xd, int yd, double w);

    void MaximizationStep(MaskedImage<int>& target, const CImg<double>& vote);


    //destructor
    ~FeaturePaint();


};

//help function
bool findPoint(CImgList<>points, Point point, int num);

#endif
