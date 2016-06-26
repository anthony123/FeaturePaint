#include "FeaturePaint.h"
#include "Color.h"
#include <iostream>
#include <ctime>

using namespace std;

#define DIR "./dataTest/4"

int main(void){
  clock_t tic, toc;
  tic = clock();
  //to select features
  char source[] = DIR"/source.png";

  //to select stroke
  char stroke[] = DIR"/stroke.png" ;

  //to select area
  char area[] =   DIR"/area.png";

  //target contour
  char contour[] = DIR"/contour.png";

  //number of stroke we select
  int nStroke = 3;

  //number of area we select
  int nArea = 3;

  FeaturePaint fp(source, stroke, area, contour);
  fp.setNumberOfStroke(nStroke);
  fp.setNumberOfArea(nArea);

  //process stroke
  const Color s1(0, 255,79);
  fp.stroke_colors.push_back(s1);
  const Color s2(255,0,23);
  fp.stroke_colors.push_back(s2);
  const Color s3(0,255,254);
  fp.stroke_colors.push_back(s3);
  //  Color s4(255,12,176);
  //fp.stroke_colors.push_back(a4);

  const Color f1(64,64,64);
  //fp.area_colors.push_back(f1);
  const Color f2(255,128,128);
  //fp.area_colors.push_back(f2);
  const Color f3(0,255,0);
  //fp.area_colors.push_back(f3);


  cout<<"------------------process first stroke-------------------"<<endl;
  fp.processALine(s1, false);

  cout<<"------------------process second stroke-------------------"<<endl;
  //fp.processALine(s2,false);

  cout<<"------------------process third stroke-------------------"<<endl;
  //fp.processALine(s3, true);

   cout<<"------------------process fourth stroke-------------------"<<endl;
   //fp.processALine(s4);


  cout<<"------------------process first area-------------------"<<endl;
  //fp.multiScaleFill(f1);


  cout<<"------------------process second area-------------------"<<endl;
  //fp.multiScaleFill(f2);

  cout<<"------------------process third area-------------------"<<endl;
  //
  //fp.multiScaleFill(f3);

  //test efla
  //fp.inTheSameArea(Point(391,501), Point(197, 239));

  fp.afterProcess();
  toc = clock();
  cout<<"The process time is "<<(toc-tic)/CLOCKS_PER_SEC<<" seconds"<<endl;

  return 0;
}
