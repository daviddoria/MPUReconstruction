#ifndef QUADRICOCORNER_H
#define QUADRICOCORNER_H 

#include <stdio.h>
#include "PointSet.h"
#include "LocalApprox.h"
#include <math.h>

class ImplicitOctCell;

class QuadricOCorner: public LocalApprox{

public:
  int _deg;
  float *_c0, *_cx, *_cy, *_cz, *_cxx, *_cyy, *_czz, *_cxy, *_cyz, *_czx;
  bool *_op;
  
  QuadricOCorner(PointSet* ps, float R_error, float R_laf, ImplicitOctCell *cell,
                 int* index_list, int listN, int index1, int index2, int index3, bool &OK, float edgeT);
  
  ~QuadricOCorner(){
    if(_deg != 0){
      delete[] _c0;
      delete[] _cx; delete[] _cy; delete[] _cz; 
      delete[] _cxx; delete[] _cyy; delete[] _czz;
      delete[] _cxy; delete[] _cyz; delete[] _czx;
      delete[] _op;
    }
  }
  
  void swapCoefficient(int i, int j){
    float tmp;
    tmp = _c0[i];  _c0[i] = _c0[j];    _c0[j] = tmp;
    
    tmp = _cx[i];  _cx[i] = _cx[j];    _cx[j] = tmp;
    tmp = _cy[i];  _cy[i] = _cy[j];    _cy[j] = tmp;
    tmp = _cz[i];  _cz[i] = _cz[j];    _cz[j] = tmp;
    
    tmp = _cxx[i]; _cxx[i] = _cxx[j];  _cxx[j] = tmp;
    tmp = _cyy[i]; _cyy[i] = _cyy[j];  _cyy[j] = tmp;
    tmp = _czz[i]; _czz[i] = _czz[j];  _czz[j] = tmp;
    
    tmp = _cxy[i]; _cxy[i] = _cxy[j];  _cxy[j] = tmp;
    tmp = _cyz[i]; _cyz[i] = _cyz[j];  _cyz[j] = tmp;
    tmp = _czx[i]; _czx[i] = _czx[j];  _czx[j] = tmp;
  }
  
  void computeCoefficient(PointSet* ps, float R, ImplicitOctCell *cell,
                          int* index_list, int start, int end, int f_index,
                          float o[3], float n[3]);
  
  float value(float x, float y, float z){
    float f =  _c0[0] + 
      (_cx[0] + _cxx[0]*x + _cxy[0]*y)*x +
        (_cy[0] + _cyy[0]*y + _cyz[0]*z)*y +
          (_cz[0] + _czz[0]*z + _czx[0]*x)*z;
    for(int i=1; i<_deg; i++){
      float fi = _c0[i] + 
        (_cx[i] + _cxx[i]*x + _cxy[i]*y)*x +
          (_cy[i] + _cyy[i]*y + _cyz[i]*z)*y +
            (_cz[i] + _czz[i]*z + _czx[i]*x)*z;
      
      if(_op[i-1]){
        if(f < fi)
          f = fi;
      }
      else{
        if(f > fi)
          f = fi;
      }
    }
    return f;
  }
  
  void gradient(float g[], float x, float y, float z){
    float f = _c0[0] + 
      (_cx[0] + _cxx[0]*x + _cxy[0]*y)*x +
        (_cy[0] + _cyy[0]*y + _cyz[0]*z)*y +
          (_cz[0] + _czz[0]*z + _czx[0]*x)*z;
    int index = 0;
    for(int i=1; i<_deg; i++){
      float fi = _c0[i] + 
        (_cx[i] + _cxx[i]*x + _cxy[i]*y)*x +
          (_cy[i] + _cyy[i]*y + _cyz[i]*z)*y +
            (_cz[i] + _czz[i]*z + _czx[i]*x)*z;
      
      if(_op[i-1]){
        if(f < fi)
          index = i;
      }
      else{
        if(f > fi)
          index = i;
      }
    }
    g[0] = _cx[index] + 2.0f*_cxx[index]*x + _cxy[index]*y + _czx[index]*z;
    g[1] = _cy[index] + 2.0f*_cyy[index]*y + _cxy[index]*x + _cyz[index]*z;
    g[2] = _cz[index] + 2.0f*_czz[index]*z + _cyz[index]*y + _czx[index]*x;
  }
  
  float valueAndGradient(float g[3], float x, float y, float z){
    float f = _c0[0] + 
      (_cx[0] + _cxx[0]*x + _cxy[0]*y)*x +
        (_cy[0] + _cyy[0]*y + _cyz[0]*z)*y +
          (_cz[0] + _czz[0]*z + _czx[0]*x)*z;
    int index = 0;
    for(int i=1; i<_deg; i++){
      float fi = _c0[i] + 
        (_cx[i] + _cxx[i]*x + _cxy[i]*y)*x +
          (_cy[i] + _cyy[i]*y + _cyz[i]*z)*y +
            (_cz[i] + _czz[i]*z + _czx[i]*x)*z;
      
      if(_op[i-1]){
        if(f < fi){
          f = fi;
          index = i;
        }
      }
      else{
        if(f > fi){
          f = fi;
          index = i;
        }
      }
    }
    g[0] = _cx[index] + 2.0f*_cxx[index]*x + _cxy[index]*y + _czx[index]*z;
    g[1] = _cy[index] + 2.0f*_cyy[index]*y + _cxy[index]*x + _cyz[index]*z;
    g[2] = _cz[index] + 2.0f*_czz[index]*z + _cyz[index]*y + _czx[index]*x;
    
    return f;
  }
};

#endif