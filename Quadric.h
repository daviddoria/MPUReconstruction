#ifndef QUADRIC_H
#define QUADRIC_H 

#include <stdio.h>
#include "PointSet.h"
#include "LocalApprox.h"

class ImplicitOctCell;

class Quadric: public LocalApprox{

public:
  float _c0, _cx, _cy, _cz, _cxx, _cyy, _czz, _cxy, _cyz, _czx;
  
  Quadric(PointSet* ps, float R_error, float R_laf, ImplicitOctCell *cell, int* index_list, int listN, float p_ave[3]);
  void computePolySVD(PointSet* ps, float R, ImplicitOctCell *cell, int* indexList, int ListN, float p_ave[3]);
  void computePolySVD2(PointSet* ps, float R, ImplicitOctCell *cell, int* indexList, int ListN, float p_ave[3]);
  float computeMaxError(PointSet* ps, float R, ImplicitOctCell *cell, int* indexList, int ListN);
  
  float value(float x, float y, float z){
    return _c0 + 
      (_cx + _cxx*x + _cxy*y)*x +
        (_cy + _cyy*y + _cyz*z)*y +
          (_cz + _czz*z + _czx*x)*z;
  }
  
  void gradient(float g[], float x, float y, float z){
    g[0] = _cx + 2.0f*_cxx*x + _cxy*y + _czx*z;
    g[1] = _cy + 2.0f*_cyy*y + _cyz*z + _cxy*x;
    g[2] = _cz + 2.0f*_czz*z + _czx*x + _cyz*y;
  }
  
  float valueAndGradient(float g[3], float x, float y, float z){
    g[0] = _cx + 2.0f*_cxx*x + _cxy*y + _czx*z;
    g[1] = _cy + 2.0f*_cyy*y + _cyz*z + _cxy*x;
    g[2] = _cz + 2.0f*_czz*z + _czx*x + _cyz*y;
    
    return _c0 + 
      (_cx + _cxx*x + _cxy*y)*x +
        (_cy + _cyy*y + _cyz*z)*y +
          (_cz + _czz*z + _czx*x)*z;
  }
};

#endif