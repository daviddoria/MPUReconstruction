#ifndef IMPLICITOCT_H 
#define IMPLICITOCT_H 

#include <stdio.h>
#include <math.h>
#include "PointSet.h"
#include "ImplicitFunction.h"

#include "LocalApprox.h"
#include "Quadric.h"
#include "QuadricO.h"
#include "QuadricOEdge.h"
#include "QuadricOCorner.h"
#include "AxisKdTree.h"

#include <math.h>

#define TWOSQRT3 3.46410161513775f

class ImplicitOctCell{
public:
  ImplicitOctTree* _tree;
  
  LocalApprox* _laf;
  bool _checked;
  float _error;
  
  ImplicitOctCell** _sub_cell;
  bool _leaf;
  unsigned char _level; 
  
  float _cx, _cy, _cz;
  float _size;
  
  ImplicitOctCell(ImplicitOctTree* tree, unsigned char level, float cx, float cy, float cz, float size);
  ~ImplicitOctCell();
  
  void evaluateFunction(float x, float y, float z, 
                        int until, float e, double &wf, double &w);
  void evaluateGradient(float x, float y, float z, int until, float e, 
                        double &wf, double &w, double gwf[], double gw[]);
  
  float supportSize();
  void computeLA();
  void split();
  
  //Is not used
  void buildFunction(int until, float e);
  
  //Is not used
  void computeConfidence(int target_level, float T);
  
  //Is not used
  void smoothNormal(int level, float (*temp_normal)[3]);
  
  inline void cellCenter(float c[3]){
    c[0] = _cx;
    c[1] = _cy;
    c[2] = _cz;
  }
  
  inline double weight(double d, float R){
    if(d > R)
      return 0;
    
    //B-spline (degree = 2)
    else{
      d = 1.5*(d/R);
      if(d < 0.5)
        return (-d*d + 0.75);
      else
        return (0.5*(1.5-d)*(1.5-d));
    }
  }
  
  inline void weightG(double g[3], float x, float y, float z, double d, float R){
    if(d > R){
      g[0] = g[1] = g[2] = 0;	
    }
    //B-spline (degree = 2)
    else{
      d = 1.5*(d/R);
      double c;
      if(d < 0.5)
        c = -4.5/(R*R);
      else
        c = -2.25*(1.5-d)/(d*R*R);
      g[0] = c*x;
      g[1] = c*y;
      g[2] = c*z;
    }
  }
  
  double theWeight(float x, float y, float z){
    double d = sqrt((x-_cx)*(x-_cx) + (y-_cy)*(y-_cy) + (z-_cz)*(z-_cz));
    float R = supportSize();
    
    return weight(d, R);
  }
  
  inline bool isIn(float x, float y, float z){
    if(x < _cx-_size || y < _cy-_size || z < _cz-_size || 
       x > _cx+_size || y > _cy+_size || z > _cz+_size)
      return false;
    else
      return true;
  }
  
  void freeChind(){
    if(!_leaf){
      for(int i=0; i<8; i++)
        delete _sub_cell[i];
      delete[] _sub_cell;
    }
    _leaf = true;
    _sub_cell = NULL;
  }
  
  void evaluateLevel(float x, float y, float z, 
                     int until, float e, double &wl, double &w){
    double d = sqrt((x-_cx)*(x-_cx) + (y-_cy)*(y-_cy) + (z-_cz)*(z-_cz));
    float R = supportSize();
    if(d > R)
      return;
    
    if(!_checked){
      computeLA();
      
      if(_error > e && _laf != NULL && _level < until){
        delete _laf;
        _laf = NULL;
      }
    }
    
    if(_error > e && _level < until){
      if(_leaf)
        split();
      
      for(int i=0; i<8; i++)
        _sub_cell[i]->evaluateLevel(x, y, z, until, e, wl, w);
    }	
    else{
      double w1 = weight(d, R);
      wl += w1* _level;
      w += w1;
    }
  }
};

class ImplicitOctTree{
public:
  PointSet* _ps;
  AxisKdTree* _tree;
  ImplicitOctCell* _root;
  
  int _visit_counter;
  int _cell_count;
  
  float _support;
  int _Nmin;
  float _edgeT;
  float _cornerT;
  float _lambda;
  bool _sharp;
  
  ImplicitOctTree(PointSet* ps, float min[], float max[]){
    _ps = ps;
    
    int N = ps->_pointN;
    float size = max[0] - min[0];
    if(size < max[1] - min[1])
      size = max[1] - min[1];
    if(size < max[2] - min[2])
      size = max[2] - min[2];
    _ps->_tree = this;
    _root = new ImplicitOctCell(this, 0, 0.5f*(max[0]+min[0]),
                                0.5f*(max[1]+min[1]),
                                0.5f*(max[2]+min[2]),
                                0.5f*size);
    
    _cell_count = 1;
    
    _tree = new AxisKdTree(ps);
    
    //recommended values
    _support = 0.75;
    _Nmin = 15;
    _edgeT = 0.9f;
    _cornerT = 0.7f;
    _lambda = 0.1f;
    _sharp = true;
  }
  
  ~ImplicitOctTree(){
    delete _root;
    delete _tree;
  }
  
  float value(float x, float y, float z, int max_level, float tol){
    if(!_root->isIn(x,y,z))
      return -100000;
    
    double wf = 0;
    double w = 0;
    
    _root->evaluateFunction(x, y, z, max_level, tol, wf, w);
    
    return (float)(wf/w);
  }
  
  float level(float x, float y, float z, int max_level, float tol){
    if(!_root->isIn(x,y,z))
      return 0;
    
    double wl = 0;
    double w = 0;
    
    _root->evaluateLevel(x, y, z, max_level, tol, wl, w);
    return (float)(wl/w);
  }
  
  void build(int max_level, float tol){
    _root->buildFunction(max_level, tol);
  }
  
  void gradient(float g[3], float x, float y, float z, int max_level, float tol){
    if(!_root->isIn(x,y,z)){
      g[0] = g[1] = g[2] = 0;
      return;
    }
    
    double wf = 0;
    double w = 0;
    double gwf[] = {0,0,0};
    double gw[] = {0,0,0};
    
    _root->evaluateGradient(x, y, z, max_level, tol, wf, w, gwf, gw);
    
    g[0] = (float)(gwf[0]/w - wf*gw[0]/(w*w));
    g[1] = (float)(gwf[1]/w - wf*gw[1]/(w*w));
    g[2] = (float)(gwf[2]/w - wf*gw[2]/(w*w));
  }
  
  float valueAndGradient(float g[3], float x, float y, float z, int max_level, float tol){
    if(!_root->isIn(x,y,z)){
      g[0] = g[1] = g[2] = 0;
      return -100000;
    }
    
    double wf = 0;
    double w = 0;
    double gwf[] = {0,0,0};
    double gw[] = {0,0,0};
    
    _root->evaluateGradient(x, y, z, max_level, tol, wf, w, gwf, gw);
    
    g[0] = (float)((w*gwf[0] - wf*gw[0])/(w*w));
    g[1] = (float)((w*gwf[1] - wf*gw[1])/(w*w));
    g[2] = (float)((w*gwf[2] - wf*gw[2])/(w*w));
    
    return (float)(wf/w);
  }
  
  int collectPointsInSphere(int* &index_list, float c[], float R){
    int listN;
    _tree->collectPointIndexInSphere(index_list, listN, c, R);
    return listN;
  }
  
  void smoothNormal(int level){
    int N = _ps->_pointN;
    float (*temp_normal)[3] = new float[N][3];
    
    for(int i=0; i<N; i++)
      temp_normal[i][0] = temp_normal[i][1] = temp_normal[i][2] = 0;
    
    _root->smoothNormal(level, temp_normal);
    
    float (*normal)[3] = _ps->_normal;
    for(int i=0; i<N; i++){
      float* m = temp_normal[i];
      double len = sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
      if((float)len != 0){
        normal[i][0] = (float)(m[0]/len);
        normal[i][1] = (float)(m[1]/len);
        normal[i][2] = (float)(m[2]/len);
      }
    }
    delete[] temp_normal;
  }
  
  void computeConfidence(int level, float T){
    int N = _ps->_pointN;
    bool* bound = _ps->_bound;
    for(int i=0; i<N; i++)
      bound[i] = false;
    
    _root->computeConfidence(level, T);
  }
};

class ImplicitPOU : public ImplicitFunction{
public:
  ImplicitOctTree* imp;
  int _max_level;
  float _tol;
  float _iso;
  
  int count;
  
  ~ImplicitPOU(){
    delete imp;
  }
  
  //Not used.
  void setMaxCellN(int max){
    imp->_visit_counter = 0;
    //imp->_max_cell_N = max;
  }
  
  float value(float x, float y, float z){
    count++;
    return imp->value(x, y, z, _max_level, _tol) + _iso;
  }
  
  void build(){
    imp->build(_max_level, _tol);
  }
  
  float valueAndGradient(float g[3], float x, float y, float z){
    count++;
    return imp->valueAndGradient(g, x, y, z, _max_level, _tol);
  }
  
  void gradient(float g[3], float x, float y, float z){
    count++;
    imp->gradient(g, x, y, z, _max_level, _tol);
  }
  
  int getCellCount(){
    return imp->_cell_count;
  }
  
  bool isIn(float x, float y, float z){
    return imp->_root->isIn(x,y,z);
  }
  
  float level(float x, float y, float z){
    count++;
    return imp->level(x, y, z, _max_level, _tol);
  }
  
  //Should be |d|=1
  bool searchZeroOnLine(float p[], float n[], float o[], float d[]){
    //compute intersection of line and bounding sphere 
    float c[3], R;
    imp->_root->cellCenter(c);
    R = imp->_root->_size;
    float a[] = {o[0]-c[0], o[1]-c[1], o[2]-c[2]};
    float dot = a[0]*d[0] + a[1]*d[1] + a[2]*d[2];
    float D = dot*dot - (a[0]*a[0] + a[1]*a[1] + a[2]*a[2] - R*R);
    if(D <= 0)
      return false;
    
    //two end points
    float ts = -dot - (float)sqrt(D);
    float te = -dot + (float)sqrt(D);
    
    if(ts < 0){
      if(te < 0)
        return false;
      else
        ts = 0;
    }
    
    p[0] = o[0] + ts*d[0];
    p[1] = o[1] + ts*d[1];
    p[2] = o[2] + ts*d[2];
    float f1 = value(p[0], p[1], p[2]);
    float dt = (float)fabs(f1);
    for(float t = ts + dt; t < te; t += dt){
      p[0] = o[0] + t*d[0];
      p[1] = o[1] + t*d[1];
      p[2] = o[2] + t*d[2];
      
      float f2 = value(p[0], p[1], p[2]);
      //rare case
      if(f1*f2 < 0){
        float t1 = t - dt;
        float t2 = t;
        float t3 = (float)((fabs(f2)*t1 + fabs(f1)*t2)/(fabs(f1)+fabs(f2)));
        for(int i=0; i<20; i++){
          p[0] = o[0] + t3*d[0];
          p[1] = o[1] + t3*d[1];
          p[2] = o[2] + t3*d[2];
          float f = valueAndGradient(n, p[0], p[1], p[2]);
          
          if(fabs(f) < 0.001)
            break;
          
          if(f1*f < 0){
            f2 = f;
            t2 = t3;
          }
          else{
            f1 = f;
            t1 = t3;
          }
          
          float df = n[0]*d[0] + n[1]*d[1] + n[2]*d[2];
          t = t3 - f/df;
          if(t1 > t || t > t2 || df == 0)
            t3 = (float)((fabs(f2)*t1 + fabs(f1)*t2)/(fabs(f1)+fabs(f2)));
          else
            t3 = t;
        }
        float len = (float)sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
        if(len == 0)
          return false;
        n[0] /= len;
        n[1] /= len;
        n[2] /= len;
        return true;
      }
      else{
        f1 = f2;
        dt = (float)fabs(f1);
        if(dt < 0.001){
          gradient(n, p[0], p[1], p[2]);
          float len = (float)sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
          if(len == 0)
            return false;
          n[0] /= len;
          n[1] /= len;
          n[2] /= len;
          return true;
        }
      }
    }
    return false;
  }
};

#endif