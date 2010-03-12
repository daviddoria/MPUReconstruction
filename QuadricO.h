#ifndef QUADRICO_H
#define QUADRICO_H 

#include <stdio.h>
#include "PointSet.h"
#include "LocalApprox.h"

class ImplicitOctCell;

class QuadricO: public LocalApprox{

public:
  float _c0, _cx, _cy, _cz, _cxx, _cyy, _czz, _cxy, _cyz, _czx;
  
  QuadricO(PointSet* ps, float R_error, float R_laf, ImplicitOctCell *cell, 
           int* index_list, int listN, float p_ave[3], float n_ave[3]);
  
  float value(float x, float y, float z){
    return _c0 + 
      (_cx + _cxx*x + _cxy*y)*x +
        (_cy + _cyy*y + _cyz*z)*y +
          (_cz + _czz*z + _czx*x)*z;	   
  }
  
  void gradient(float g[], float x, float y, float z){
    g[0] = _cx + 2.0f*_cxx*x + _cxy*y + _czx*z;
    g[1] = _cy + 2.0f*_cyy*y + _cxy*x + _cyz*z;
    g[2] = _cz + 2.0f*_czz*z + _cyz*y + _czx*x;
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