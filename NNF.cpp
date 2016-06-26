#include <cstdlib>
#include <iostream>
#include "NNF.h"
#include "MaskedImage.h"
#include "FeaturePaint.h"
#include "Point.h"
#include "CImg.h"

using namespace cimg_library;

template<typename T> void MaskedImage<T>::setMask(int x, int y, bool masked)
{
  mask(x,y,0) = masked;
}

template<typename T> CImg<bool> MaskedImage<T>::getMask() const
{
  return mask;
}

template<typename T> bool MaskedImage<T>::isMasked(Point pt) const
{
  return mask(pt.x, pt.y, 0);
}


NNF::NNF(const MaskedImage<int> &source_img, const MaskedImage<int>&target_img)
{
  input = source_img;
  output = target_img;

  int W = input.width();
  int H = input.height();

  field = CImg<int>(W,H,1,3);
  cimg_forXY(field, x,y){
    field(x,y,2) = DSCALE;
  }
}

void NNF::setInput(const MaskedImage<int>& input_img, bool sourceToTarget)
{
  input = input_img;
  //initialize(sourceToTarget);
}

void NNF::setOutput(const MaskedImage<int>& output_img, bool sourceToTarget)
{
  output = output_img;
  //initialize(sourceToTarget);
}

void NNF::initFromOtherNNF(const NNF& otherNNF, bool sourceToTarget)
{
  //typically from a half-sized nnf
  int H = input.height();
  int W = input.width();
  int Ho = output.height();
  int Wo = output.width();
  //field = CImg<int>(W, H, 1, 3);
  //cout<<"in the initFromOtherNNF H = "<<H<<", W = "<<W<<", Hother = "<<otherNNF.getInput().height()<<", Whother = "<<otherNNF.getInput().width()<<endl;
  for (int x=0;x<W;++x) {
    for (int y=0;y<H;++y) {
      if(!input.isMasked(Point(x,y))){
	//int xlow = int(min(x/2, W-1));
	//int ylow = int(min(y/2, H-1));
          int xlow = x/2;
          int ylow = y/2;

	 field(x,y,0) = max(min(otherNNF.getXField(xlow, ylow)*2, Wo),0);
	  field(x,y,1) = max(min(otherNNF.getYField(xlow, ylow)*2, Ho),0);
	 field(x,y,2) = distance(x,y, field(x,y,0), field(x,y,1));
  // int xt = max(min(otherNNF.getXField(xlow, ylow)*2, Wo),0);
  // int yt = max(min(otherNNF.getYField(xlow, ylow)*2, Ho),0);
  // int dt = distance(x,y, xt, yt);
  // if(field(x,y,2) > dt){
  //   field(x,y,0) = xt;
  //   field(x,y,1) = yt;
  //   field(x,y,2) = dt;
  // }
  //field(x,y,2) = DSCALE;
  //cout<<"xt  = "<<field(x,y,0)<<", yt = "<<field(x,y,1)<<endl;
      }
    }
  }
  //cout<<"in the initFromOtherNNF, about to enter into initialize"<<endl;
  initialize(sourceToTarget);
}

NNF::NNF(const NNF& nnf)
{
  input = nnf.input;
  output = nnf.output;
  field = nnf.field;
}

int NNF::distance1(int sourceX, int sourceY, int targetX, int targetY) const
{


  int dis = 0;
  int r1,g1,b1,r2,g2,b2, dx,dy;
  int  S = FILL_PATCH/2;
  srand(time(NULL));
  for ( dy=-S ; dy<=S ; ++dy ) {
    for ( dx=-S ; dx<=S ; ++dx ) {
      int xs = sourceX+dx;
      int ys = sourceY+dy;

      int xt = sourceX+dx;
      int yt = sourceY+dy;

      if(xs < 0 || xs >= input.width()|| ys < 0|| ys >= input.height()|| input.isMasked(Point(xs,ys))){
	r1 = rand()%255;
	g1 = rand()%255;
	b1 = rand()%255;
      }else{
	r1 = input(xs, ys, 0);
	g1 = input(xs, ys, 1);
	b1 = input(xs, ys, 2);
      }

      if(xt < 0 || xt >= output.width() ||yt<0||yt >= output.height()|| output.isMasked(Point(xt,yt))){
	r2 = rand()%255;
	g2 = rand()%255;
	b2 = rand()%255;
      }else{
	r2 = output(xt,yt,0);
	g2 = output(xt,yt,1);
	b2 = output(xt,yt,2);
      }

      dis += ((r1-r2)*(r1-r2) + (g1-g2)*(g1-g2) + (b1 - b2)*(b1 - b2));
    }
  }
  //cout<<"dis = "<<dis<<endl;

  return dis;
}

int  NNF::distance(int sourceX, int sourceY, int targetX, int targetY) const
{
  long double distance=0;
  long double wsum=0, ssdmax = 9*255*255;
  int dy, dx, band;
  int xks, yks;
  int xkt, ykt;
  long double ssd;
  long res;
  int s_value, t_value, s_gx, t_gx, s_gy, t_gy;
  int S = FILL_PATCH/2;

  // for each pixel in the source patch
  for ( dy=-S ; dy<=S ; ++dy ) {
    for ( dx=-S ; dx<=S ; ++dx ) {

      //source position
      xks = sourceX+dx;
      yks = sourceY+dy;

      //target position
      xkt=targetX+dx;
      ykt=targetY+dy;

      //number of pixels to add up
      wsum++;

      //if not in the img range
      if ( xks<1 || xks>=input.width()-1 ) {distance++; continue;}
      if ( yks<1 || yks>=input.height()-1 ) {distance++; continue;}

      // cannot use masked pixels as a valid source of information
      if (input.isMasked(Point(xks, yks))) {distance++; continue;}

      // corresponding pixel in the target patch
      if (xkt<1 || xkt>=output.width()-1) {distance++; continue;}
      if (ykt<1 || ykt>=output.height()-1) {distance++; continue;}

      // cannot use masked pixels as a valid source of information
      if (output.isMasked(Point(xkt, ykt))) {distance++; continue;}

      ssd=0;
      for (band=0; band<3; ++band) {
	// pixel values
	s_value =  input(xks, yks, band);
	t_value = output(xkt, ykt, band);

	// pixel horizontal gradients (Gx)
	s_gx = 128 + ( input(xks+1, yks, band) -  input( xks-1, yks, band))/2;
	t_gx = 128 + (output(xkt+1, ykt, band) - output( xkt-1, ykt, band))/2;

	// pixel vertical gradients (Gy)
	//s_gy = 128+(getSampleMaskedImage(source, xks, yks+1, band) - getSampleMaskedImage(source, xks, yks-1, band))/2;
	//t_gy = 128+(getSampleMaskedImage(target, xkt, ykt+1, band) - getSampleMaskedImage(target, xkt, ykt-1, band))/2;
	s_gy = 128 + ( input(xks, yks+1, band) -  input(xks, yks-1, band))/2;
	t_gy = 128 + (output(xkt, ykt+1, band) - output(xkt, ykt-1, band))/2;

	ssd += pow((long double)s_value-t_value , 2); // distance between values in [0,255^2]

	ssd += pow((long double)s_gx-t_gx , 2); // distance between Gx in [0,255^2]
	ssd += pow((long double)s_gy-t_gy , 2); // distance between Gy in [0,255^2]
      }

      // add pixel distance to global patch distance
      distance += ssd/ssdmax;
      //cout<<"ssd = "<<ssd<<", ssdmax = "<<ssdmax<<", distance = "<<distance<<endl;
    }
  }
  //cout<<"in the distance function, distance = "<<distance<<endl;
  //cout<<"in the distance function, wsum = "<<wsum<<endl;
  res = (int)(DSCALE*distance/wsum);
  //cout<<"in the distance function, distance = "<<distance<<endl;
  //cout<<"in the distance function, wsum = "<<wsum<<endl;
  //cout<<"in the distance function, res = "<<res<<endl;
  if (res < 0 || res > DSCALE) return DSCALE;
  return res;

}

void NNF::randomInit(bool sourceToTarget)
{
  srand(time(NULL));
  int W = input.width();
  int H = input.height();
  int Ho = output.height();
  int Wo = output.width();

  for(int i = 0; i < W; i++){
    for(int j = 0; j < H; j++){
      if(input.isMasked(Point(i,j)) == false){
	int xt,yt;
	do{xt = rand()%Wo;
	  yt = rand()%Ho;
	  field(i, j, 0) = xt;
	  field(i, j, 1) = yt;
	  field(i, j, 2) = distance(i,j,xt,yt);
	}while(output.isMasked(Point(xt,yt)));
      }
    }
  }

  initialize(sourceToTarget);
}


void NNF::initialize(bool sourceToTarget)
{
  int H = input.height();
  int W = input.width();
  int ho = output.height();
  int wo = output.width();

  //std::cout<<"H = "<<H<<", W = "<<W<<", ho = "<<ho<<", wo = "<<wo<<endl;
  int max_iter = max(output.width(), output.height());
  srand(time(NULL));


  //only deal with unmaksed pixels
  for(int x = 0; x < W; x++){
    for(int y = 0; y < H; y++){
      if(!input.isMasked(Point(x,y))) {
	field(x,y,2) = distance(x,y,field(x,y,0), field(x,y,1));
	int iter = 0;

	while(field(x,y,2) == DSCALE && iter < max_iter){
	  ++iter;

	  field(x,y,0) = rand()%wo;
	  field(x,y,1) = rand()%ho;
	  field(x,y,2) = distance(x,y,field(x,y,0), field(x,y,1));
	  //cout<<"in the while loop xo = "<<field(x,y,0)<<", yo = "<<field(x,y,1)<<endl;
	}
  //cout<<"x = "<<x<<", w = "<<W<<", y = "<<y<<", H = "<<H<<", x0 = "<<field(x,y,0)<<", w0 = "<<wo<<", y = "<<field(x,y,1)<<", h0 = "<<ho<<endl;
    if(sourceToTarget){
      output(field(x,y,0), field(x,y,1),0) = input(x,y,0);
      output(field(x,y,0), field(x,y,1),1) = input(x,y,1);
      output(field(x,y,0), field(x,y,1),2) = input(x,y,2);
    }else{
      input(x,y,0) = output(field(x,y,0), field(x,y,1),0);
      input(x,y,1) = output(field(x,y,0), field(x,y,1),1);
      input(x,y,2) = output(field(x,y,0), field(x,y,1),2);
    }
      }
      //cout<<"field = "<<field(x,y,2)<<endl;
    }
  }
}

NNF& NNF::operator=(const NNF& nnf){
  input = nnf.input;
  output = nnf.output;
  field = nnf.field;
  return *this;
}

MaskedImage<int> NNF::getInput() const
{
  return input;

}

MaskedImage<int> NNF::getOutput() const
{
  return output;
}

void NNF::minimizeNNF(int pass)
{
  int i, y, x;
  int min_x=0, min_y=0, max_y=input.height()-1, max_x=input.width()-1;

  // multi-pass minimization
  for (i=0;i<pass;i++) {
    // scanline order
    for (x=min_x;x<= max_x;++x)
      for (y=min_y;y<=max_y;++y)
	if (!input.isMasked(Point(x,y)))
	  minimizeLinkNNF(x,y,+1);

    // reverse scanline order
    for (x=max_x;x>=min_x;x--)
      for (y=max_y;y>=min_y;y--)
	if (!input.isMasked(Point(x,y)))
	  minimizeLinkNNF(x,y,-1);
  }
}

void NNF::minimizeLinkNNF(int x, int y, int dir)
{
  int xp,yp,dp,wi, xpi, ypi;

  //Propagation Up/Down
  if (x-dir>0 && x-dir<input.width()) {
    xp = field(x-dir,y,0)+dir;
    yp = field(x-dir,y,1);

    //must ensure that (xp, yp) is in the image
    dp = distance(x,y, xp,yp);
    if (dp<field(x,y,2) && xp < output.width() && yp < output.height()&& xp >=0 && yp >= 0&&!output.isMasked(Point(xp,yp))) {
      field(x,y,0) = xp;
      field(x,y,1) = yp;
      field(x,y,2) = dp;
    }
  }
  //Propagation Left/Right
  if (y-dir>0 && y-dir<input.height()) {
    xp = field(x,y-dir,0);
    yp = field(x,y-dir,1) + dir;
    dp = distance(x,y, xp,yp);

    if (dp<field(x,y,2) && xp < output.width() && xp >=0 && yp >= 0 && yp < output.height()&&!output.isMasked(Point(xp,yp))) {
      field(x,y,0) = xp;
      field(x,y,1) = yp;
      field(x,y,2) = dp;
    }
  }

  //Random search
  wi=output.width();

  xpi = field(x,y,0);
  ypi = field(x,y,1);
  int r=0;
  while (wi>0) {
    r=(rand() % (2*wi)) -wi;
    xp = xpi + r;
    r=(rand() % (2*wi)) -wi;
    yp = ypi + r;
    xp = int(max(0, min(output.width()-1, xp )));
    yp = int(max(0, min(output.height()-1, yp )));

    dp = distance(x,y, xp,yp);
    if (dp<field(x,y,2) && !output.isMasked(Point(xp,yp))) {
      field(x,y,0) = xp;
      field(x,y,1) = yp;
      field(x,y,2) = dp;
    }
    wi/=2;
  }
  if(field(x,y,2) == DSCALE){
    cout<<"distance of ("<<x<<","<<y<<") is "<<DSCALE<<endl;
  }
}

int NNF::getXField(int x, int y) const
{
  return field(x,y,0);
}

int NNF::getYField(int x, int y) const
{
  return field(x,y,1);
}

int NNF::getDField(int x, int y) const
{
  return field(x, y, 2);
}
