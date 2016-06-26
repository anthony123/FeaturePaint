#ifndef NNF_H_
#define NNF_H_

#include "MaskedImage.h"

using namespace cimg_library;

class NNF{
 private:
  MaskedImage<int> input;
  MaskedImage<int> output;
  CImg<int>field;

 public:
  NNF(){input = MaskedImage<int>(); output = MaskedImage<int>(); field = CImg<int>();}
  NNF(const MaskedImage<int>& source, const MaskedImage<int>& target);
  NNF(const NNF& nnf);
  void randomInit(bool sourceToTarget);
  int distance(int sourceX, int sourceY, int targetX, int targetY) const;
  int distance1(int sourceX, int sourceY, int targetX, int targetY) const;
  void initialize(bool sourceToTarget);
  NNF& operator=(const NNF& nnf);
  MaskedImage<int> getInput() const ;
  MaskedImage<int> getOutput() const ;
  void minimizeNNF(int pass);
  void minimizeLinkNNF(int x, int y, int dir);
  int getXField(int x, int y) const;
  int getYField(int x, int y) const;
  int getDField(int x, int y) const;
  void initFromOtherNNF(const NNF& otherNNF, bool sourceToTarget);
  void setInput(const MaskedImage<int>& img, bool sourceToTarget);
  void setOutput(const MaskedImage<int>& img, bool sourceToTarget);
};
#endif
