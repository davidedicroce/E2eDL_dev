#include "E2eDL/FrameProducers/interface/FrameStriding.h"

// Striding input frames (vDetFrame by rowstrides and colstrides accordingly)
e2e::Frame1D frameStriding(std::vector<float>& vDetFrame, int rows, int columns, int rowstrides, int colstrides){
  e2e::Frame1D vStridedFrame ((rows*rowstrides)*(columns*colstrides),0);
  for (int rowidx=0; rowidx<rows; rowidx++){
    for (int colidx=0; colidx<columns; colidx++){
      for (int kernelrow=0; kernelrow<rowstrides; kernelrow++){
        for (int kernelcol=0; kernelcol<colstrides; kernelcol++){
          vStridedFrame[(rowstrides*rowidx+kernelrow)*columns+colstrides*colidx+kernelcol] = vDetFrame[rowidx*columns+colidx]/(rowstrides*colstrides);
          //if(rowidx<5 && colidx<5) std::cout<<"("<<rowstrides*rowidx+kernelrow<<","<<colstrides*colidx+kernelcol<<"): "<<vStridedFrame[rowstrides*rowidx+kernelrow][colstrides*colidx+kernelcol]<<" "<<vDetFrame[rowidx*columns+colidx]/(rowstrides*colstrides);
        }
      }
    }
  }
  return vStridedFrame;
}
