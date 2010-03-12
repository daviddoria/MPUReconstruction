#ifndef QUADRICOEDGE_H
#define QUADRICOEDGE_H 

#include <stdio.h>
#include "PointSet.h"
#include "LocalApprox.h"
#include <math.h>

class ImplicitOctCell;

class QuadricOEdge: public LocalApprox{

public:
  float _c10, _c1x, _c1y, _c1z, _c1xx, _c1yy, _c1zz, _c1xy, _c1yz, _c1zx;
  float _c20, _c2x, _c2y, _c2z, _c2xx, _c2yy, _c2zz, _c2xy, _c2yz, _c2zx;
  bool _is_convex;
  
  QuadricOEdge(PointSet* ps, float R_error, float R_laf, ImplicitOctCell *cell,
               int* index_list, int listN, int index1, int index2);
  
  void computeCoefficient(PointSet* ps, float R, ImplicitOctCell *cell,
                          int* index_list, int start, int end, int f_index,
                          float o[3], float n[3]);
  
  float value(float x, float y, float z);
  
  void gradient(float g[], float x, float y, float z){
    float f1 = _c10 + _c1x*x + _c1y*y + _c1z*z + _c1xx*x*x + _c1yy*y*y + _c1zz*z*z +_c1xy*x*y + _c1yz*y*z + _c1zx*z*x;
    float f2 = _c20 + _c2x*x + _c2y*y + _c2z*z + _c2xx*x*x + _c2yy*y*y + _c2zz*z*z +_c2xy*x*y + _c2yz*y*z + _c2zx*z*x;
    
    if(_is_convex){
      if(f1 > f2){
        g[0] = _c1x + 2.0f*_c1xx*x + _c1xy*y + _c1zx*z;
        g[1] = _c1y + 2.0f*_c1yy*y + _c1xy*x + _c1yz*z;
        g[2] = _c1z + 2.0f*_c1zz*z + _c1yz*y + _c1zx*x;
      }
      else{
        g[0] = _c2x + 2.0f*_c2xx*x + _c2xy*y + _c2zx*z;
        g[1] = _c2y + 2.0f*_c2yy*y + _c2xy*x + _c2yz*z;
        g[2] = _c2z + 2.0f*_c2zz*z + _c2yz*y + _c2zx*x;
      }
    }
    else{
      if(f1 < f2){
        g[0] = _c1x + 2.0f*_c1xx*x + _c1xy*y + _c1zx*z;
        g[1] = _c1y + 2.0f*_c1yy*y + _c1xy*x + _c1yz*z;
        g[2] = _c1z + 2.0f*_c1zz*z + _c1yz*y + _c1zx*x;
      }
      else{
        g[0] = _c2x + 2.0f*_c2xx*x + _c2xy*y + _c2zx*z;
        g[1] = _c2y + 2.0f*_c2yy*y + _c2xy*x + _c2yz*z;
        g[2] = _c2z + 2.0f*_c2zz*z + _c2yz*y + _c2zx*x;
      }
    }
  }
  
  float valueAndGradient(float g[3], float x, float y, float z){
    float f1 = _c10 + 
      (_c1x + _c1xx*x + _c1xy*y)*x +
        (_c1y + _c1yy*y + _c1yz*z)*y +
          (_c1z + _c1zz*z + _c1zx*x)*z;
    float f2 = _c20 + 
      (_c2x + _c2xx*x + _c2xy*y)*x +
        (_c2y + _c2yy*y + _c2yz*z)*y +
          (_c2z + _c2zz*z + _c2zx*x)*z;
    
    if(_is_convex){
      if(f1 > f2){
        g[0] = _c1x + 2.0f*_c1xx*x + _c1xy*y + _c1zx*z;
        g[1] = _c1y + 2.0f*_c1yy*y + _c1xy*x + _c1yz*z;
        g[2] = _c1z + 2.0f*_c1zz*z + _c1yz*y + _c1zx*x;
        
        return f1;
      }
      else{
        g[0] = _c2x + 2.0f*_c2xx*x + _c2xy*y + _c2zx*z;
        g[1] = _c2y + 2.0f*_c2yy*y + _c2xy*x + _c2yz*z;
        g[2] = _c2z + 2.0f*_c2zz*z + _c2yz*y + _c2zx*x;
        
        return f2;
      }
    }
    else{
      if(f1 < f2){
        g[0] = _c1x + 2.0f*_c1xx*x + _c1xy*y + _c1zx*z;
        g[1] = _c1y + 2.0f*_c1yy*y + _c1xy*x + _c1yz*z;
        g[2] = _c1z + 2.0f*_c1zz*z + _c1yz*y + _c1zx*x;
        
        return f1;
      }
      else{
        g[0] = _c2x + 2.0f*_c2xx*x + _c2xy*y + _c2zx*z;
        g[1] = _c2y + 2.0f*_c2yy*y + _c2xy*x + _c2yz*z;
        g[2] = _c2z + 2.0f*_c2zz*z + _c2yz*y + _c2zx*x;
        
        return f2;
      }
    }
  }
};

#endif