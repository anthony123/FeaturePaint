#include <vector>
#include <set>
#include <iostream>
#include <cassert>
#include <ctime>
#include <stack>
#include <algorithm>
#include <armadillo>
#include "MaskedImage.h"
#include "NNF.h"
#include "Range.h"
#include "CImg.h"
#include "Point.h"
#include "FeaturePaint.h"
#include "Utility.h"

using namespace cimg_library;
using namespace std;
using namespace arma;

template<typename T> void MaskedImage<T>::setMask(int x, int y, bool masked)
{
  mask(x,y,0) = masked;
}

template<typename T> void MaskedImage<T>::setMaskImage(const CImg<bool>&img){
  mask = img;
}

template<typename T> CImg<bool> MaskedImage<T>::getMask() const
{
  return mask;
}

template<typename T> bool MaskedImage<T>::isMasked(Point pt) const
{
  return mask(pt.x, pt.y, 0);
}

template<typename T> void MaskedImage<T>::halfSize()
{
  CImg<T>::resize_halfXY();
  int newH = mask.height()/2;
  int newW = mask.width()/2;

  CImg<bool>new_mask(newW, newH, 1,1,true);
  cimg_forXY(mask, x,y){
    if(mask(x,y,0) == false){
      int newX = x/2;
      int newY = y/2;
      //std::cout<<"newX = "<<newX<<", newY = "<<newY<<std::endl;
      //std::cout<<"newH = "<<newH<<", newW = "<<newW<<std::endl;
      new_mask(newX, newY, 0)  = false;
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

template<typename T> MaskedImage<T>& MaskedImage<T>::operator=(const MaskedImage& img)
{

  CImg<T>::operator=(img);
  mask  = img.getMask();

  return *this;
}


//constructor
FeaturePaint::FeaturePaint(char* source, char* stroke, char* area, char* contour){
  source_image.assign(source);
  stroke_image.assign(stroke);
  area_image.assign(area);
  contour_image.assign(contour);
  simi = Utility::generateD2S();
  target_image.assign(contour_image.width(),contour_image.height(),1,3).fill(255);

  TARGET_IMAGE_WIDTH = contour_image.width();
  TARGET_IMAGE_HEIGHT = contour_image.height();
}

//set number of stroke
void FeaturePaint::setNumberOfStroke(int n){
  nStroke = n;
}

//destructor
FeaturePaint::~FeaturePaint()
{
  free(simi);
}

//set number of area
void FeaturePaint::setNumberOfArea(int n){
  nArea = n;
}

CImgList<int>FeaturePaint::getStroke(CImg<>img, const Color& color, bool isClosed){
  CImgList<int>points;
  const int r = color(0);
  const int g = color(1);
  const int b = color(2);

  cimg_forXY(img, x, y){
    if(img(x,y,0) == r  && img(x,y,1) == g && img(x,y,2) == b){
      int nei = numberOfNeighbors(img, Point(x,y), color);

      //first node
      if(isClosed || nei == 1){
	Point curr(x,y);
	Point first(x,y);
	//insert the first node
	points.insert(CImg<int>::vector(curr.x,curr.y));
	while(1){
	  int index = points.size()-1;
	  curr = Point(points(index,0), points(index,1));
	  Point next =  nextPoint(img, points, curr, color);
	  if(next == Point(-1,-1) || next == first){
	    //the last point

	    //refine points
	    int i = 1;
	    while(i < points.size()-1){
	      if((points(i,0) == points(i-1,0) && points(i,1) == points(i+1,1))||(points(i,1) == points(i-1,1) && points(i,0) == points(i+1,0))){
		points.remove(i);
	      }else{
		++i;
	      }
	    }
	    //for debug
	    /*cout<<"----------------------------------------"<<endl;
	      cimglist_for(points, p){
	      int x = points(p, 0);
	      int y = points(p, 1);
	      cout <<"["<<x<<", "<<y<<"]"<<endl;
	      }
	      cout<<"----------------------------------------"<<endl;*/
	    return points;
	  }

	  points.insert(CImg<int>::vector(next.x,next.y));

	}

      }

    }
  }
  strokes.push_back(points);
  return points;
}


//find out how many neighbors in the same color
int FeaturePaint::numberOfNeighbors(CImg<>img, Point point, const Color& color){
  int res = 0;
  int x = point.x;
  int y = point.y;
  int xx, yy;
  int
    red = color(0),
    green = color(1),
    blue = color(2);
  int width = img.width();
  int height = img.height();

  //upper left
  xx = x-1;
  yy = y-1;

  //for debug

  if(xx >= 0 && yy >= 0 && img(xx,yy,0) == red && img(xx,yy,1) == green && img(xx,yy,2) == blue)
    ++res;

  //upper
  xx = x;
  yy = y-1;

  if(yy >= 0 &&  img(xx,yy,0) == red && img(xx,yy,1) == green && img(xx,yy,2) == blue)
    ++res;

  //upper right
  xx = x+1;
  yy = y-1;

  if(xx<width && yy>=0 &&  img(xx,yy,0) == red && img(xx,yy,1) == green && img(xx,yy,2) == blue)
    ++res;

  //left
  xx = x-1;
  yy = y;

  if(xx>=0 && img(xx,yy,0) == red && img(xx,yy,1) == green && img(xx,yy,2) == blue)
    ++res;

  //right
  xx = x+1;
  yy = y;

  if(xx < width && img(xx,yy,0) == red && img(xx,yy,1) == green && img(xx,yy,2) == blue)
    ++res;

  //lower left
  xx = x-1;
  yy = y+1;

  if(xx >= 0 && yy <height &&img(xx,yy,0) == red && img(xx,yy,1) == green && img(xx,yy,2) == blue)
    ++res;

  //lower
  xx = x;
  yy = y+1;

  if(xx < width && yy <height && img(xx,yy,0) == red && img(xx,yy,1) == green && img(xx,yy,2) == blue)
    ++res;

  //lower right
  xx = x+1;
  yy = y+1;
  if(xx < width && yy < height && img(xx,yy,0) == red && img(xx,yy,1) == green && img(xx,yy,2) == blue)
    ++res;


  return res;
}

//find next point in the same color
Point FeaturePaint::nextPoint(CImg<>img, CImgList<int> points, Point curr, Color color){
  Point res;
  int x = curr.x;
  int y = curr.y;
  int xx,yy;
  int width = img.width();
  int height = img.height();
  int
    red = color(0),
    green = color(1),
    blue = color(2);

  int size = points.size();  //>= 1

  //left
  xx = x-1;
  yy = y;

  if(xx>=0 && img(xx,yy,0) == red && img(xx,yy,1) == green && img(xx,yy,2) == blue){
    int num = (size>10)?10:size;
    bool found = findPoint(points, Point(xx,yy), num);
    if(!found)
      return Point(xx,yy);
  }

  //right
  xx = x+1;
  yy = y;

  if(xx < width && img(xx,yy,0) == red && img(xx,yy,1) == green && img(xx,yy,2) == blue){
    int num = (size>10)?10:size;
    bool found = findPoint(points, Point(xx,yy), num);
    if(!found)
      return Point(xx,yy);
  }

  //upper
  xx = x;
  yy = y-1;
  if(xx<width && yy>=0 &&  img(xx,yy,0) == red && img(xx,yy,1) == green && img(xx,yy,2) == blue){
    int num = (size>10)?10:size;
    bool found = findPoint(points, Point(xx,yy), num);
    if(!found)
      return Point(xx,yy);
  }

  //lower
  xx = x;
  yy = y+1;
  if(xx < width && yy <height &&img(xx,yy,0) == red && img(xx,yy,1) == green && img(xx,yy,2) == blue){
    int num = (size>10)?10:size;
    bool found = findPoint(points, Point(xx,yy), num);
    if(!found)
      return Point(xx,yy);
  }

  //upper left
  xx = x-1;
  yy = y-1;

  if(xx >= 0 && yy >= 0 && img(xx,yy,0) == red && img(xx,yy,1) == green && img(xx,yy,2) == blue){
    int num = (size>10)?10:size;
    bool found = findPoint(points, Point(xx,yy), num);
    if(!found)
      return Point(xx,yy);
  }


  //upper right
  xx = x+1;
  yy = y-1;

  if(xx<width && yy>=0 &&  img(xx,yy,0) == red && img(xx,yy,1) == green && img(xx,yy,2) == blue){
    int num = (size>10)?10:size;
    bool found = findPoint(points, Point(xx,yy), num);
    if(!found)
      return Point(xx,yy);
  }

  //lower left
  xx = x-1;
  yy = y+1;

  if(xx >= 0 && yy <height &&img(xx,yy,0) == red && img(xx,yy,1) == green && img(xx,yy,2) == blue){
    int num = (size>10)?10:size;
    bool found = findPoint(points, Point(xx,yy), num);
    if(!found)
      return Point(xx,yy);
  }

  //lower right
  xx = x+1;
  yy = y+1;
  if(xx < width && yy < height && img(xx,yy,0) == red && img(xx,yy,1) == green && img(xx,yy,2) == blue){
    int num = (size>10)?10:size;
    bool found = findPoint(points, Point(xx,yy), num);
    if(!found)
      return Point(xx,yy);
  }


  return Point(-1,-1);
}


//sample the stroke
CImgList<int>FeaturePaint::sampleStroke(CImgList<int>stroke){
  CImgList<int>sample;
  //cout<<"-------------------------------"<<endl;
  int SAMPLE_WIDTH = PATCH_WIDTH/2;
  //SAMPLE_WIDTH = (SAMPLE_WIDTH>=2)?SAMPLE_WIDTH:2;
  cimglist_for(stroke, index){
    if(!(index%SAMPLE_WIDTH)){
      int x = stroke(index,0);
      int y = stroke(index, 1);
      //for debug

      // cout<<"[ "<<x<<", "<<y<<"]";
      sample.insert(CImg<>::vector(x,y));
    }
  }

  return sample;
}

//process one line: selection and synthesis
void FeaturePaint::processALine(const Color& color, bool isClosed){

  //get user's stroke
  CImgList<int> original_stroke = getStroke(stroke_image,color, false);
  cout<<"get stroke done"<<endl;

  //adjust user's input-stroke to contour
  CImgList<int> stroke = adhereToContour(original_stroke);
  //CImgList<int> stroke = original_stroke;

  //sample the stroke
  CImgList<int> sample_stroke = sampleStroke(stroke);

  //get the stroke of contour_image
  CImgList<int>contour_stroke = getStroke(contour_image, color, isClosed);
  cout<<"get contour stroke done"<<endl;

  //sample the contour stroke
  CImgList<int>sample_contour_stroke = sampleStroke(contour_stroke);

  cout<<"the original target points is "<<contour_stroke.size()<<", the sampled target points is "<<sample_contour_stroke.size()<<endl;
  

  //cout<<"the number of  source patch is "<<sample_stroke.size()<<endl;
  int target_length = sample_contour_stroke.size();
  int source_length = sample_stroke.size();
  //cout<<"source length is "<<source_length<<endl;
  int curr_size = 0;
  vector<int>target_path;

  srand(time(NULL));
  while(1){
    int start = rand()%source_length;
    //left direction or right direction
    int direction = rand()%2;
    if(direction == 0){
      for(int i = start; i < source_length; i++){
	target_path.push_back(i);
	++curr_size;
	if(curr_size >= target_length)
	  break;
      }
    }else{
      for(int i = start; i >= 0; --i){
	target_path.push_back(i);
	++curr_size;
	if(curr_size >= target_length)
	  break;
      }
    }

    if(curr_size >= target_length)
      break;
  }

  cout<<"path length = "<<target_path.size()<<", contour length = "<<sample_contour_stroke.size()<<endl;
  assert(target_path.size() == sample_contour_stroke.size());

  render(target_path, sample_stroke, sample_contour_stroke);
  //cout<<"mapping process done"<<endl;

}

void FeaturePaint::render(const vector<int> &path, CImgList<int>sample_stroke, CImgList<int>sample_contour_stroke)
{
  vector<Point>tpts;
  vector<Point>spts;

  int sz = path.size();

  int source_sz = sample_stroke.size();
  for(int i = 0; i < source_sz; i++){
    int x1 = sample_stroke(i,0);
    int y1 = sample_stroke(i,1);
    spts.push_back(Point(x1,y1));
  }

  for(int i = 0; i < sz; i++){

    int x2 = sample_contour_stroke(i,0);
    int y2 = sample_contour_stroke(i,1);
    tpts.push_back(Point(x2,y2));
  }

  //added on 2016/4/3
  //warp the target points into a straight line
  std::vector<Point> straight_target_points;
  auto temp_target = CImg<int>(1024, 1024,1,3,255);
  int gap = PATCH_WIDTH/2;
  int target_size = tpts.size();

  for(int i = 0; i < target_size; i++){
    straight_target_points.push_back(Point(64+i*gap, 64));
  }

  for(int i = 0; i < straight_target_points.size();i++){
    int x = straight_target_points[i].x;
    int y = straight_target_points[i].y;

    temp_target(x,y,0) = 0;
    temp_target(x,y,1) = 0;
    temp_target(x,y,2) = 0;
    
  }
  
  // for(Point point: tpts){
  //   //straight_target_points.push_back(Point(point.x, tpts[0].y));
  //   straight_target_points.push_back(Point(point.x, 64));
    
  // }

  //find all target points in the stroke
  findAllStrokePoints(straight_target_points, temp_target);

  std::vector<Point> corresponding_source_points;
  for(const auto point : source_stroke_points){
    int index = findPatch(point, straight_target_points);
    assert(index != -1);

    int index_of_source_patch = path[index];
    //int index_of_source_patch = index_of_target_patch;
    Point corresponding = findCorrespondingPixel(point, index, index_of_source_patch, straight_target_points,spts);
    corresponding_source_points.push_back(corresponding);
  }
  for(int i = 0; i < source_stroke_points.size(); i++){
    mapPixel(source_stroke_points[i], corresponding_source_points[i], temp_target, source_image);
  }
  
   temp_target.save_png("temp_target.png");


  //using MLS to warp the straight line
  MLSWarping(straight_target_points,tpts);
  cout<<"MLSWarping done"<<endl;
  //map the pixels
  assert(source_stroke_points.size() == target_stroke_points.size());
  
  for(int i = 0; i < source_stroke_points.size(); i++){
    cout<<"the "<<i+1<<" iterations"<<endl;
    Point point = target_stroke_points[i];
    if(point.x >= 0 && point.x < target_image.width() && point.y >= 0 && point.y < target_image.height()){
      mapPixel( target_stroke_points[i],corresponding_source_points[i], target_image, source_image);
    }else{
      //use four neighbors pixels to fill
      
    }
   
  }

  // //sort by x axis
  // vector<Point>tempx = tpts;
  // sort(tempx.begin(), tempx.end(), Utility::compx);
  //
  // //sort by y axis
  // vector<Point>tempy = tpts;
  // sort(tempy.begin(), tempy.end(), Utility::compy);
  //
  // //most left x and toppest y
  // int most_left_x = (tempx[0].x-PATCH_WIDTH >= 0)?tempx[0].x-PATCH_WIDTH:0;
  // int most_top_y = (tempy[0].y-PATCH_HEIGHT >= 0)?tempy[0].y-PATCH_HEIGHT:0;
  //
  //
  // //most right x and lowest y
  // int most_right_x = (tempx[sz-1].x+PATCH_WIDTH < TARGET_IMAGE_WIDTH)?(tempx[sz-1].x+PATCH_WIDTH):(TARGET_IMAGE_WIDTH-1);
  // int most_low_y = (tempy[sz-1].y+PATCH_HEIGHT < TARGET_IMAGE_HEIGHT)?(tempy[sz-1].y+PATCH_HEIGHT):(TARGET_IMAGE_HEIGHT-1);
  // //cout<<"width = "<<TARGET_IMAGE_WIDTH<<", height = "<<TARGET_IMAGE_HEIGHT<<endl;
  // //cout<<"Most left x is "<<most_left_x<<", most right x is "<<most_right_x<<", most top y is "<<most_top_y<<", most low y is "<<most_low_y<<endl;
  // //find points in range Point([most_left_x, most_right_x], [most_top_y, most_low_y]) belong to which patch
  // for(int xx = most_left_x; xx <= most_right_x; xx++){
  //   for(int yy = most_top_y; yy <= most_low_y; yy++){
  //     int index_of_target_patch = findPatch(Point(xx, yy), tpts);
  //     //for debug
  //     //cout<<"the index of target patch is "<<index_of_target_patch<<endl;
  //     if(index_of_target_patch != -1){
	// int index_of_source_patch = path[index_of_target_patch];
	// //int index_of_source_patch = index_of_target_patch;
	// Point corresponding = findCorrespondingPixel(Point(xx,yy), index_of_target_patch, index_of_source_patch, tpts,spts);
	// mapPixel(Point(xx,yy), corresponding, target_image, source_image);
  //     }
  //   }
  // }

  
  //  //print out straight stroke tripe
  //  CImg<>straight_image(contour_image.width(), contour_image.height(), 1, 3);
  //  straight_image.fill(255);
  //  Point origin = Point(64,64);
  //  int SAMPLE_WIDTH = PATCH_WIDTH/4;
  //  SAMPLE_WIDTH = (SAMPLE_WIDTH>=2)?SAMPLE_WIDTH:2;
  //  int sour_sz = spts.size();
  //  vector<Point>straight_points;
  //  for(int i = 0; i < sour_sz; i++){
  //    straight_points.push_back(Point(origin.x+i*SAMPLE_WIDTH, origin.y));
  //  }
  
  //  cimg_forXY(straight_image,xx,yy){
  //    int index = findPatch(Point(xx,yy), straight_points);
  //    if(index != -1){
  //      Point corresponding = findCorrespondingPixel(Point(xx,yy), index, index, straight_points, spts);
  //      mapPixel(Point(xx,yy), corresponding, target_image, source_image);
  //    }
  //  }
  // straight_image.save_png("straight.png");
}

void FeaturePaint::MLSWarping(const vector<Point>& source_control_points, const vector<Point>&target_control_points)
{
    //1. convert vectors to matrix
    int source_size = source_control_points.size();
    mat p(2, source_size);

    for(int i = 0; i < source_size; i++){
      p(0, i)  = source_control_points[i].x;
      p(1, i)  = source_control_points[i].y;
    }

  int target_size = target_control_points.size();
  
  mat q(2, target_size);

  for(int i = 0; i < target_size; i++){
    q(0, i) = target_control_points[i].x;
    q(1, i) = target_control_points[i].y;
  }

  int all_points_size = source_stroke_points.size();
  mat v(2, all_points_size);

  for (int i = 0; i < all_points_size; i++) {
    v(0,i) = source_stroke_points[i].x;
    v(1,i) = source_stroke_points[i].y;
  }

  cout<<"convert to matrix done!"<<endl;

  //print out source control point and target control point
  CImg<int>test(1024,1024,1,3,255);
  for(int i = 0; i < target_size; i++){
    int x = p(0,i);
    int y = p(1,i);
    test(x,y,0) = 0;
    
  }

  for(int i = 0; i < source_size; i++){
    int x = q(0,i);
    int y = q(1,i);

    test(x,y,1) = 0;
  }
  test.save_png("test.png");
  assert(source_size == target_size);

  //2. precompute W
  mat w = Utility::precomputeWeights(p,v,2.0);
  cout<<"compute w done"<<endl;

  
    
  //3. precomputeRigid
  typeRigid mlsd = Utility::precomputeRigid(p, v, w);
  cout<<"pre-computed rigid done"<<endl;

  //4. PointsTransformRigid
  mat fv = Utility::PointsTransformRigid(w,mlsd, q);
  cout << "Rigid transform done"<<endl;

  cout<<"fv:("<<fv.n_cols<<","<<fv.n_rows<<")"<<endl;
  int number_of_columns = fv.n_cols;
  for (int i = 0; i < number_of_columns; i++) {
    int x = fv(0,i);
    int y = fv(1,i);
    target_stroke_points.push_back(Point(x,y));
  }

}

CImgList<int> FeaturePaint::adhereToContour(CImgList<int>stroke)
{
  int sz = stroke.size();
  vector<Point>pts;

  for(int i = 0; i < sz; i++){
    int xx = stroke(i,0);
    int yy = stroke(i, 1);

    pts.push_back(Point(xx,yy));
  }

  CImgList<int>result;

  //gradient image
  CImg<float> grad = gradientImage(source_image);
  //grad.save_png("grad.png");

  for(int i = 0; i < sz; i++){
    Point normal = Utility::findNormal(i, pts);
    int xx,yy;
    for(int j = 1; j <= RANGE; j++){
      Point p1 = pts[i] + normal*j;
      Point p2 = pts[i] - normal*j;
      xx = static_cast<int>(p1.x);
      yy = static_cast<int>(p1.y);
      //cout<<"grad is "<<grad(xx, yy)<<endl;
      if(xx >= 0 && xx < TARGET_IMAGE_WIDTH && yy >= 0 && yy < TARGET_IMAGE_HEIGHT && static_cast<int>(grad(xx, yy)) > THRESHOLD){
	result.insert(CImg<int>::vector(xx,yy));
	break;
      }

      xx = static_cast<int>(p2.x);
      yy = static_cast<int>(p2.y);
      //cout<<"grad is "<<grad(xx, yy)<<endl;
      if(xx >= 0 && xx < TARGET_IMAGE_WIDTH && yy >= 0 && yy < TARGET_IMAGE_HEIGHT && static_cast<int>(grad(xx, yy)) > THRESHOLD){
	result.insert(CImg<int>::vector(xx,yy));
	break;
      }

    }
  }
  //cout<<"the original size is "<<stroke.size()<<", the result size is "<<result.size()<<endl;
  return result;

}

CImg<float>FeaturePaint::gradientImage(CImg<int>img)
{
  CImgList<float> grad = img.get_gradient("XY");

  int width = grad[0].width();
  int height = grad[0].height();
  CImg<float>result(width, height, 3);
  for(int i = 0; i < width; i++){
    for(int j = 0; j < height; j++){
      result(i,j, 0) = sqrt(grad(0, i,j, 0)*grad(0, i,j, 0)+grad(1,i,j, 0)*grad(1,i,j, 0));
      result(i,j, 1) = sqrt(grad(0, i,j, 1)*grad(0, i,j, 1)+grad(1,i,j, 1)*grad(1,i,j, 1));
      result(i,j, 2) = sqrt(grad(0, i,j, 2)*grad(0, i,j, 2)+grad(1,i,j, 2)*grad(1,i,j, 2));

    }
  }

  CImg<>res = result.normalize(0,255);

  return res;
}

void FeaturePaint::findAllStrokePoints(const vector<Point>&spts, const CImg<int>& image)
{

  cimg_forXY(image, x, y){
    if(findPatch(Point(x,y), spts) != -1)
      source_stroke_points.push_back(Point(x,y));
  }

  cout<<"in the finAllStrokePoints function, the source points size is "<<spts.size()<<", the returned points is "<<source_stroke_points.size()<<endl;

  
}



void FeaturePaint::mapPixel(Point target, Point source, CImg<int>&target_image, CImg<int>&source_image)
{

  // cout<<"Mapping pixels\n";
  cout<<"source: ("<<source.x<<", "<<source.y<<")";
  cout<<"target: ("<<target.x<<", "<<target.y<<")"<<endl;
  int r = source_image(source.x, source.y, 0);
  int g = source_image(source.x, source.y, 1);
  int b = source_image(source.x, source.y, 2);
  // cout<<"("<<r<<","<<g<<","<<b<<")\n";
  target_image(target.x, target.y, 0) = r;
  target_image(target.x, target.y, 1) = g;
  target_image(target.x, target.y, 2) = b;
}

void FeaturePaint::mapPatchPixel(Point target, Point source, CImg<int>&target_image, CImg<int>&source_image)
{
  for(int i = 0; i < FILL_PATCH; i++){
    for(int j = 0; j < FILL_PATCH; j++){
      int r = source_image(source.x+i, source.y+j, 0);
      int g = source_image(source.x+i, source.y+j, 1);
      int b = source_image(source.x+i, source.y+j, 2);

      target_image(target.x+i, target.y+j, 0) = r;
      target_image(target.x+i, target.y+j, 1) = g;
      target_image(target.x+i, target.y+j, 2) = b;
    }
  }
}

Point FeaturePaint::findCorrespondingPixel(Point p, int index_of_target_patch, int index_of_source_patch, const vector<Point>&tpts, const vector<Point>&spts){

  Point source_center_point = spts[index_of_source_patch];
  Point target_center_point = tpts[index_of_target_patch];
  Point source_normal = Utility::findNormal(index_of_source_patch, spts);
  Point source_tangent = Utility::findTangent(index_of_source_patch, spts);
  Point target_normal = Utility::findNormal(index_of_target_patch, tpts);
  Point target_tangent = Utility::findTangent(index_of_target_patch, tpts);

  Point d_target = p - target_center_point;
  mat A(2,2);
  A(0, 0) = target_normal.x;
  A(1, 0) = target_normal.y;
  A(0, 1) = target_tangent.x;
  A(1, 1) = target_tangent.y;

  vec b = {d_target.x, d_target.y};
  vec x = solve(A,b);

  Point res = source_center_point + source_normal*x(0) + source_tangent*x(1);
  return Point(static_cast<int>(res.x), static_cast<int>(res.y));

}

//find which patch the point p is belonging to
int FeaturePaint::findPatch(Point p, const vector<Point>& pts)

{
  //find the point which is nearest to the p
  int sz = pts.size();
  float min_distance = 2*(TARGET_IMAGE_WIDTH+TARGET_IMAGE_HEIGHT);
  int min_index = -1;
  for(int i = 0; i < sz; i++){
    float dis = Utility::computeDistance(p, pts[i]);
    if(dis < min_distance){
      min_distance = dis;
      min_index = i;
    }
  }

  //test if p in the rectangle centered in pts[min_index]
  if(Utility::pointInRectangle(p, min_index, pts) == true)
    return min_index;
  else
    return -1;
}

void FeaturePaint::afterProcess(){
  target_image.save_png("target.png");
}

//find the source filling color, return the sorted points
//set<Point>FeaturePaint::readArea(const Color& color, bool isFilledArea)
//unordered_set<Point, Utility::hash_func>FeaturePaint::readArea(CImg<int>& img, const Color& color, bool isFilledArea)
vector<Point>FeaturePaint::readArea(CImg<int>& img, const Color& color, bool isFilledArea)
{

  //set<Point>res;
  vector<Point>res;
  int
    r = color(0),
    g = color(1),
    b = color(2);

  Color white(255,255,255);
  if(!isFilledArea){
    cimg_forXY(img, x, y){
      if(img(x,y,0) == r && img(x,y,1) == g && img(x,y,2) == b){
	res.push_back(Point(x,y));
      }
    }

  }else{

    //1. find the position of target color point
    //Maybe a bug here, what if there two area to be filled
    //const Point& target_color_point = findTargetColorPoint(color, isFilledArea);
    //cout<<"target color point is ["<<target_color_point.x<<", "<<target_color_point.y<<"]\n";
    //2. for all points,find all points are in the same region with the target point

    stack<Point>pts;

    //push target point to pts
    pts.push(target_color_point);
    res.push_back(target_color_point);

    img(target_color_point.x, target_color_point.y, 0) = 250;
    img(target_color_point.x, target_color_point.y, 1) = 0;
    img(target_color_point.x, target_color_point.y, 2) = 0;

    int tr,tg,tb;
    while(!pts.empty()){
      Point temp = pts.top();
      int x = temp.x;
      int y = temp.y;

      //cout<<"read ["<<x<<", "<<y<<"]\n";
      //cout<<"the stack size is "<<pts.size()<<endl;
      pts.pop();
      if((x+1) < TARGET_IMAGE_WIDTH){
	tr = img(x+1, y, 0);
	tg = img(x+1, y, 1);
	tb = img(x+1, y, 2);
	Color curr_color(tr,tg,tb);
	//cout<<"["<<(x+1)<<", "<<(y)<<"]: ("<<tr<<", "<<tg<<", "<<tb<<")\n";
	if((curr_color == color)||(curr_color == white)){
	  pts.push(Point(x+1, y));
	  res.push_back(Point(x+1, y));
	  img(x+1, y, 0) = 250;
	  img(x+1, y, 1) = 0;
	  img(x+1, y, 2) = 0;

	}

      }

      if((x-1) >= 0){
	tr = img(x-1, y, 0);
	tg = img(x-1, y, 1);
	tb = img(x-1, y, 2);
	Color curr_color(tr,tg,tb);
	//cout<<"["<<(x-1)<<", "<<y<<"]: ("<<tr<<", "<<tg<<", "<<tb<<")\n";
	if((curr_color == color)||(curr_color == white)){
	  pts.push(Point(x-1, y));
	  res.push_back(Point(x-1,y));
	  img(x-1, y, 0) = 250;
	  img(x-1, y, 1) = 0;
	  img(x-1, y, 2) = 0;
	}

      }

      if((y+1) < TARGET_IMAGE_HEIGHT){
	tr = img(x, y+1, 0);
	tg = img(x, y+1, 1);
	tb = img(x, y+1, 2);
	Color curr_color(tr,tg,tb);
	//cout<<"["<<(x)<<", "<<(y+1)<<"]: ("<<tr<<", "<<tg<<", "<<tb<<")\n";
	if((curr_color == color)||(curr_color == white)){
	  pts.push(Point(x, y+1));
	  res.push_back(Point(x, y+1));
	  img(x, y+1, 0) = 250;
	  img(x, y+1, 1) = 0;
	  img(x, y+1, 2) = 0;
	}

      }

      if((y-1) > 0){
	tr = img(x, y-1, 0);
	tg = img(x, y-1, 1);
	tb = img(x, y-1, 2);
	Color curr_color(tr,tg,tb);
	//cout<<"["<<(x)<<", "<<(y-1)<<"]: ("<<tr<<", "<<tg<<", "<<tb<<")\n";
	if((curr_color == color)||(curr_color == white)){
	  pts.push(Point(x, y-1));
	  res.push_back(Point(x,y-1));
	  img(x, y-1, 0) = 250;
	  img(x, y-1, 1) = 0;
	  img(x, y-1, 2) = 0;
	}

      }
    }

    //debug
    // int res_sz = res.size();
    // for(int i = 0; i < res_sz; i++){
    //   int xx = res[i].x;
    //   int yy = res[i].y;

    //   target_image(xx,yy,0) = 250;
    //   target_image(xx,yy,1) = 0;
    //   target_image(xx,yy,2) = 0;
    // }
    for(auto& elem: res){
      int x = elem.x;
      int y = elem.y;
      if(img(x,y,0) == 250 && img(x,y,1) == 0 && img(x,y,2) == 0){
	img(x,y,0) = 255;
	img(x,y,1) = 255;
	img(x,y,2) = 255;
      }
    }
  }


  return res;
}

//find target color position
Point FeaturePaint::findTargetColorPoint(Color color, bool isFilledArea)
{
  int
    r = color(0),
    g = color(1),
    b = color(2);

  Point target;
  if(isFilledArea){
    cimg_forXY(contour_image,x,y){
      if(contour_image(x,y,0) == r && contour_image(x,y,1) == g && contour_image(x,y,2) == b){
	target = Point(x,y);
	break;
      }
    }
  }else{
    cimg_forXY(area_image,x,y){
      // if(area_image(x,y,0) == 255 && area_image(x,y,1)==255&&area_image(x,y,2)==0){
      // 	//cout<<"haha, i found the target color\n";
      // }
      if(area_image(x,y,0) == r && area_image(x,y,1) == g && area_image(x,y,2) == b){
	target = Point(x,y);
	break;
      }
    }
  }

  return target;
}

void FeaturePaint::multiScaleFill(const Color& target_color){



  MaskedImage<int> new_target;

  /** 1. generate source and target masked images **/
  vector<Point> filling_pts, filled_pts;
  target_color_point = findTargetColorPoint(target_color, true);

  //read the filling area
  filling_pts = readArea(area_image, target_color,false);


  //read the filled area
  filled_pts = readArea(contour_image, target_color, true);
  int source_width = source_image.width();
  int source_height = source_image.height();

  int target_width = target_image.width();
  int target_height = target_image.height();

  CImg<bool>source_mask(source_width, source_height, 1,1,true);
  CImg<bool>target_mask(target_width, target_height, 1, 1, true);

  for (auto pt:filling_pts){
    int x = pt.x;
    int y = pt.y;

    source_mask(x,y,0) = false;
  }

  for(auto pt:filled_pts){
    int x = pt.x;
    int y = pt.y;

    target_mask(x,y,0) = false;
  }

  MaskedImage<int> source(source_image, source_mask);
  MaskedImage<int> target(target_image, target_mask);

  /** 2.generate levels of MaskedImage **/
  int max_level = 0;
  vector<MaskedImage<int>>sources;
  vector<MaskedImage<int>>targets;

  while (true) {
    int min_source_size = min(source.width(), source.height());
    int min_target_size = min(target.width(), target.height());
    //cout<<"min_source = "<<min_source_size<<", min_target = "<<min_target_size<<endl;
    if (min_target_size >= FILL_PATCH/2 && min_source_size >= FILL_PATCH/2) {
      ++max_level;

      sources.push_back(source);
      targets.push_back(target);
      //cout<<"source.width = "<<source.width()<<", source.height = "<<source.height()<<endl;
      //cout<<"target.width = "<<target.width()<<", target.height = "<<target.height()<<endl;


      //half size source and target
      source.halfSize();
      target.halfSize();

    }else{

      break;
    }
  }

  cout<<"the level of images is "<<max_level<<endl;



  /** for each level,  init+prepagation&random search **/
  NNF S2Tnnf, T2Snnf, new_S2Tnnf, new_T2Snnf;

  for (int level = max_level - 1; level >= 0; --level) {
    cout<<"********level progress "<<level<<": "<<max_level<<"************"<<endl;

    //1. init
    if (level == max_level-1) {
      S2Tnnf = NNF(sources[level], targets[level]);
      // S2Tnnf = NNF(sources[level], targets[level]);
      T2Snnf = NNF(targets[level], sources[level]);

      //just random init
      //cout<<"****************random init and asign pixel values*************\n";
      S2Tnnf.randomInit(true);
      T2Snnf.randomInit(false);
    }else{
      //cout<<"new_target.width = "<<new_target.width()<<", new_target.height = "<<new_target.height()<<endl;
      //cout<<"level = "<<level<<endl;
      //cout<<"source.height = "<<source.height()<<", source.width = "<<source.width()<<endl;
      //new target is updated Is it reasonable?

      //new_S2Tnnf = NNF(sources[level], new_target);
      new_S2Tnnf = NNF(sources[level], targets[level]);
      //new_S2Tnnf.randomInit(true);
      new_S2Tnnf.initFromOtherNNF(S2Tnnf,true);
      S2Tnnf = new_S2Tnnf;

      //new_T2Snnf = NNF(new_target, sources[level]);
      new_T2Snnf = NNF(targets[level], sources[level]);
      //new_T2Snnf.randomInit(false);
      new_T2Snnf.initFromOtherNNF(T2Snnf,false);
      T2Snnf = new_T2Snnf;
      //S2Tnnf.setOutput(new_target,true);
      //T2Snnf.setInput(new_target,false);
      //cout<<"*****************init from other nnf*******\n";

      //init from smaller-sized nnf



      //cout<<"this nnf size:("<<new_S2Tnnf.getInput().width()<<","<<new_S2Tnnf.getInput().height()<<")"<<endl;
      // cout<<"previous nnf size:("<<S2Tnnf.getInput().width()<<", "<<S2Tnnf.getInput().height()<<")"<<endl;



      //cout<<"*******************init from other nnf********\n";

    }

    //update targets in targets
    //targets[level] = new_target;
    // cout<<"******************Just After initialize Step**************\n";
    //   for(int x = 0; x < S2Tnnf.getInput().width(); x++){
    //   for(int y = 0; y < S2Tnnf.getInput().height(); y++){
    //     if(!S2Tnnf.getInput().isMasked(Point(x,y)))
    //       cout<<"("<<x<<","<<y<<"):"<<S2Tnnf.getXField(x,y)<<","<<S2Tnnf.getYField(x,y)<<endl;
    //     }
    //   }
    cout<<"init done"<<endl;
    //Propagation & random search return best target image  of each level further update S2Tnnf and T2Snnf
    MaskedImage<int> *imgptr =  ExpectationMaximization(sources, targets,S2Tnnf,T2Snnf,level);
    cout<<"**********EM process done**************\n";
    new_target = *imgptr;
    /*for(int x = 0; x < new_target.width(); x++){
      for(int y = 0; y < new_target.height(); y++){
      if (new_target.isMasked(Point(x,y)) == false)
      cout<<"("<<x<<","<<y<<"): "<<"("<<new_target(x,y,0)<<","<<new_target(x,y,1)<<","<<new_target(x,y,2)<<")\n";
      }
      }*/
    delete imgptr;
  }



  //target_image = CImg<int>(new_target);
  //target_image = CImg<int>(S2Tnnf.getOutput());

  for(int y = 0; y < new_target.height(); y++){
    for(int x = 0; x < new_target.width(); x++){
        if(!new_target.isMasked(Point(x,y))){
          target_image(x,y,0) = new_target(x,y,0);
          target_image(x,y,1) = new_target(x,y,1);
          target_image(x,y,2) = new_target(x,y,2);
          //cout<<"("<<x<<","<<y<<"): "<<"("<<new_target(x,y,0)<<","<<new_target(x,y,1)<<","<<new_target(x,y,2)<<")\n";
          }
        }
     }
     cout<<"S2Tnnf.getOutput().height = "<<S2Tnnf.getOutput().height()<<"S2Tnnf.getOutput().width = "<<S2Tnnf.getOutput().width()<<endl;

  cout<<"multi-scale and patch match done\n";
  cout<<"**************Quality Measurement*****************"<<endl;
  long dist = 0;
  int count = 0;
  for(int y = 0; y < new_target.height(); y++){
    for(int x = 0; x < new_target.width(); x++){
        if(!new_target.isMasked(Point(x,y))){
          ++count;
          int xs = T2Snnf.getXField(x,y);
          int ys = T2Snnf.getYField(x,y);
          dist += T2Snnf.distance(x,y,xs,ys);
          //cout<<"("<<x<<","<<y<<"): "<<"("<<new_target(x,y,0)<<","<<new_target(x,y,1)<<","<<new_target(x,y,2)<<")\n";
          }
        }
     }
     cout<<"the average distance is "<<dist/count<<endl;
}


MaskedImage<int>* FeaturePaint::ExpectationMaximization(vector<MaskedImage<int>>sources, vector<MaskedImage<int>>targets, NNF& S2Tnnf,NNF& T2Snnf,int level)
{
  int iterEM = (1 + 2*level);
  int iterNNF = int(min(7, 1 + level));

  // bool upscaled;

  //MaskedImage<int> new_source;
  //MaskedImage<int> source = S2Tnnf.getInput();
  //MaskedImage<int> target = S2Tnnf.getOutput();
  MaskedImage<int> source = sources[level];
  MaskedImage<int> target = targets[level];
  //cout<<"in the EM, level = "<<level<<endl;
  //cout<<"S2Tnnf.getInput().height() = source.height(): "<<S2Tnnf.getInput().height()<<" == "<<source.height()<<endl;
  assert(S2Tnnf.getInput().height() == source.height());
  //cout<<"S2Tnnf.getOutput().height() = target.height(): "<<S2Tnnf.getOutput().height() <<"==" <<target.height()<<endl;
  assert(S2Tnnf.getOutput().height() == target.height());
  //cout<<"(T2Snnf.getInput().height() == target.height()): "<<T2Snnf.getInput().height()<<"="<<target.height()<<endl;
  assert(T2Snnf.getInput().height() == target.height());
  //cout<<"(T2Snnf.getOutput().height() == source.height()): "<< T2Snnf.getOutput().height() <<"="<<source.height()<<endl;
  assert(T2Snnf.getOutput().height() == source.height());

  MaskedImage<int> new_target;
  bool target_set = false;

  for(int emloop = 1; emloop <= iterEM; emloop++){
    cout<<"processing "<<emloop<<"/"<<iterEM<<endl;
    if(target_set == true){
      target = new_target;
      target.setMaskImage(targets[level].getMask());
      //S2Tnnf = NNF(source, target);
      //T2Snnf = NNF(target, source);

      //assign values or not?
      //cout<<"about to set target"<<endl;
      S2Tnnf.setOutput(target, true);
      T2Snnf.setInput(target, false);
      //cout<<"set target done"<<endl;
      //targets[level] = target;

    }

    //cout<<"in the EM, target.height() = "<<target.height()<<", target.width = "<<target.width()<<endl;
    //cout<<"in the EM, source.height() = "<<source.height()<<", source.width = "<<source.width()<<endl;

    //int H = source.height();
    //int W = source.width();

    //update (S2Tnnf S2Tnnfd T2Fnnf T2Snnfd)'s field
    S2Tnnf.minimizeNNF(iterNNF);
    T2Snnf.minimizeNNF(iterNNF);
    // cout<<"******************Just After minimizeNNF Step**************\n";
    //   for(int x = 0; x < S2Tnnf.getInput().width(); x++){
    //   for(int y = 0; y < S2Tnnf.getInput().height(); y++){
    //     if(!S2Tnnf.getInput().isMasked(Point(x,y)))
    //       cout<<"("<<x<<","<<y<<"):"<<S2Tnnf.getXField(x,y)<<","<<S2Tnnf.getYField(x,y)<<endl;
    //     }
    //   }


    //cout<<"input.height = "<<S2Tnnf.getInput().height()<<", input.width = "<<S2Tnnf.getInput().width()<<endl;
    //cout<<"output.height = "<<S2Tnnf.getOutput().height()<<", output.height = "<<S2Tnnf.getOutput().width()<<endl;

    //now upscale the original image
    bool upscaled = false;

    //cout<<"in the EM, target.width = "<<target.width()<<", target.height = "<<target.height()<<endl;

    if(level >= 1 && emloop == iterEM){
      //upscale the image
      source = sources[level-1];

      //新生成的target没用？？ 直接用上一层的target 代替？？
       //new_target = targets[level-1];
      new_target = target;
      new_target.doubleSize();
      new_target.setMaskImage(targets[level-1].getMask());
      upscaled = true;
    }else{
      source = sources[level];
      new_target = target;
      new_target.setMaskImage(targets[level].getMask());
      upscaled = false;
    }

    int newH = new_target.height();
    int newW = new_target.width();

    //cout<<"target.width = "<<target.width()<<", target.height = "<<target.height()<<endl;
    //cout<<"source.width = "<<source.width()<<", source.height = "<<source.height()<<endl;

    CImg<double>vote(newW, newH, 1, 4, 0.0);
    //cout<<"newW = "<<newW<<", newH = "<<newH<<endl;

    //update votes 如果upscaled为true，则说明已经更新了下一次的field，应该
    //更新field
    // cout<<"******************before entering Expectation Step**************\n";
    //   for(int x = 0; x < S2Tnnf.getInput().width(); x++){
    //   for(int y = 0; y < S2Tnnf.getInput().height(); y++){
    //     if(!S2Tnnf.getInput().isMasked(Point(x,y)))
    //       cout<<"("<<x<<","<<y<<"):"<<S2Tnnf.getXField(x,y)<<","<<S2Tnnf.getYField(x,y)<<endl;
    //     }
    //   }

    cout<<"Entering Expectation step from source to target"<<endl;
    ExpectationStep(S2Tnnf, 1, vote, source, upscaled);
    cout<<"Expection Step from Source to Target done\n";
    cout<<"Entering Expectation Step from target to source"<<endl;
    ExpectationStep(T2Snnf, 0, vote, source, upscaled);
    //cout<<"**********************Source pixels*************************\n";
    /*for(int x = 0; x < source.width(); x++)
      for(int y = 0; y < source.height(); y++)
      if(source.isMasked(Point(x,y)) == false)
      cout<<"("<<x<<","<<y<<"): "<<"("<<source(x,y,0)<<","<<source(x,y,1)<<", "<<source(x,y,2)<<")\n";*/

    cout<<" Expection Step from Target to Source done\n";

    //update pixel values
    cout<<"Entering Maximization Step\n";

    MaximizationStep(new_target, vote);

    cout<<"Maximization Step done\n";
    target_set = true;

    //cout<<"*****************target values**************\n";
    //int new_target_width = new_target.width();
    //int new_target_height = new_target.height();
    /*for(int x = 0; x < new_target_width; x++){
      for(int y = 0; y < new_target_height; y++){
      cout<<"("<<x<<","<<y<<"): "<<"("<<new_target(x,y,0)<<","<<new_target(x,y,1)<<", "<<new_target(x,y,2)<<")\n";
      }
      }*/

  }
  MaskedImage<int>*imgptr = new MaskedImage<int>(new_target);

  S2Tnnf.setOutput(new_target,true);
  T2Snnf.setInput(new_target,false);
  return imgptr;
}

void FeaturePaint::ExpectationStep(const NNF& nnf, int sourceToTarget, CImg<double>&vote,  const MaskedImage<int>& source, bool upscaled){

  //double *similarity = Utility::generateD2S();
  int H = nnf.getInput().height();
  int W = nnf.getInput().width();

  int Ho = nnf.getOutput().height();
  int Wo = nnf.getOutput().width();
  int source_H,source_W,target_H,target_W;

  //cout<<"in the Expectation Step, H = "<<H<<", W = "<<W<<", Ho = "<<Ho<<", Wo = "<<Wo<<endl;
  //cout<<"source.height = "<<source.height()<<", source.width = "<<source.width()<<endl;
  int R = FILL_PATCH/2;

  for(int x = 0; x < W; x++){
    for(int y = 0; y < H; y++){
      //auto pt = nnf.field(Point(x,y));
      //int dp = T2Snnf(Point(x,y));
      //cout<<"x = "<<x<<", y = "<<y<<endl;
      //if(!nnf.getInput().isMasked(Point(x,y))){
        int xp = nnf.getXField(x,y);
        int yp = nnf.getYField(x,y);
        int dp = nnf.getDField(x,y);
        //cout<<"xp = "<<xp<<",yp = "<<yp<<endl;

        double w = Utility::distanceToSimilarity(simi,dp);

	//vote for each pixel inside the input patch
	int xs, ys, xt, yt;
	for(int dy = -R; dy <= R; ++dy){
	  for(int dx = -R; dx <= R; ++dx){
	    if(sourceToTarget){
	      xs = x+dx;
	      ys = y+dy;
	      xt = xp+dx;
	      yt = yp+dy;
        source_H = H;
        source_W = W;
        target_H = Ho;
        target_W = Wo;
	    }else{
	      xs=xp+dx;
	      ys=yp+dy;
	      xt=x+dx;
	      yt=y+dy;
        target_H = H;
        target_W = W;
        source_H = Ho;
        source_W = Wo;
	    }

	    if (xs<0 || xs>=source_W) continue;
	    if (ys<0 || ys>=source_H) continue;
	    if (xt<0 || xt>=target_W) continue;
	    if (yt<0 || yt>=target_H) continue;

	    // add vote for the value
	    if (upscaled) {
	      weightedCopy(source, 2*xs,   2*ys,   vote, 2*xt,   2*yt,   w);
	      weightedCopy(source, 2*xs+1, 2*ys,   vote, 2*xt+1, 2*yt,   w);
	      weightedCopy(source, 2*xs,   2*ys+1, vote, 2*xt,   2*yt+1, w);
	      weightedCopy(source, 2*xs+1, 2*ys+1, vote, 2*xt+1, 2*yt+1, w);
	    } else {
	      weightedCopy(source, xs, ys, vote, xt, yt, w);
	    }
	  }
	}
      //}
    }
  }

}

void FeaturePaint::weightedCopy(const MaskedImage<int>& src, int xs, int ys, CImg<double>& vote, int xd, int yd, double w)
{

  if (src.isMasked(Point(xs, ys))) {
    return;
  }

  //cout<<"xd = "<<xd<<", yd = "<<yd<<", width = "<<vote.width()<<", height = "<<vote.height()<<endl;

  vote(xd,yd,0) += w*src(xs, ys, 0);
  vote(xd,yd,1) += w*src(xs, ys, 1);
  vote(xd,yd,2) += w*src(xs, ys, 2);
  vote(xd, yd,3) += w;
}

void FeaturePaint::MaximizationStep(MaskedImage<int>& target, const CImg<double>& vote)
{
  int y, x, H, W, r, g, b;
  H = target.height();
  W = target.width();
  for( x=0 ; x<W ; ++x){
    for( y=0 ; y<H ; ++y){
      if(x == 198 && y == 364)
        cout<<"vote for (198, 364) = "<<vote(x,y,3)<<endl;
      //cout<<"vote ( "<<x<<","<<y<<") = "<<vote(x,y,3)<<endl;
      if (vote(x,y,3)>0 && !target.isMasked(Point(x,y))) {
	//cout<<"change target values\n";
	r = (int) (vote(x,y,0)/vote(x,y,3));
	g = (int) (vote(x,y,1)/vote(x,y,3));
	b = (int) (vote(x,y,2)/vote(x,y,3));

	//target.setValue(x, y, 0, r );
	//target.setValue(x, y, 1, g );
	//target.setValue(x, y, 2, b );
	//target.setMask( x, y, false);
	target(x,y,0) = r;
	target(x,y,1) = g;
	target(x,y,2) = b;
	//cout<<"r = "<<r<<", g = "<<g<<", b = "<<b<<endl;
	//target.setMask(x,y,false);
      }


    }
  }
}

//4956  38082  23514

bool findPoint(CImgList<>points, Point point, int num){
  int pSize = points.size();
  bool res = false;
  for(int i = 0; i < num; i++){
    int index = pSize-1-i;
    int x = points(index,0);
    int y = points(index,1);
    if(Point(x,y) == point){
      res = true;
      break;
    }
  }

  return res;
}
