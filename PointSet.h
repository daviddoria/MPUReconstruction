#ifndef POINTSET_H 
#define POINTSET_H 

#include <cstdio>
#include <cmath>
#include <limits>
#include <vnl/vnl_math.h>

// Dummy
class ImplicitOctTree;

class PointSet
{

public:
  float (*_point)[3];
  float (*_normal)[3];
  bool *_bound;
  int _pointN;
  ImplicitOctTree* _tree;
  
public:
  PointSet()
  {
    _pointN = 0;
    _point = NULL;
    _normal = NULL;
    _bound = NULL;
  }
  
  ~PointSet()
  {
    if(_point != NULL)
      delete[] _point;
    if(_normal != NULL)
      delete[] _normal;
    if(_bound != NULL)
      delete[] _bound;
  };
  
  void setPointSize(int N)
  {
    _pointN = N;
    _point = new float[N][3];
    _normal = new float[N][3];
    //_bound = new bool[N];
    //for(int i=0; i<N; i++)
      //_bound[i] = false;
  }
  
  inline void swapIndex(int i, int j)
  {
    float tmp = _point[i][0];
    _point[i][0] = _point[j][0];
    _point[j][0] = tmp;
    
    tmp = _point[i][1];
    _point[i][1] = _point[j][1];
    _point[j][1] = tmp;
    
    tmp = _point[i][2];
    _point[i][2] = _point[j][2];
    _point[j][2] = tmp;
    
    tmp = _normal[i][0];
    _normal[i][0] = _normal[j][0];
    _normal[j][0] = tmp;
    
    tmp = _normal[i][1];
    _normal[i][1] = _normal[j][1];
    _normal[j][1] = tmp;
    
    tmp = _normal[i][2];
    _normal[i][2] = _normal[j][2];
    _normal[j][2] = tmp;
    
    //bool tmpB = _bound[i];
    //_bound[i] = _bound[j];
    //_bound[j] = tmpB;
  }
  
  void setPoint(int i, float x, float y, float z)
  {
    _point[i][0] = x;
    _point[i][1] = y;
    _point[i][2] = z;
  }
  
  void setNormal(int i, float x, float y, float z)
  {
    _normal[i][0] = x;
    _normal[i][1] = y;
    _normal[i][2] = z;
  }
  
  //Index "end" is not taken int account
  void centroid(float c[3], int start, int end)
  {
    c[0] = c[1] = c[2] = 0;
    for(int i=start; i<end; i++){
      c[0] += _point[i][0];
      c[1] += _point[i][1];
      c[2] += _point[i][2];
    }
    c[0] /= (end-start);
    c[1] /= (end-start);
    c[2] /= (end-start);
  }
  
  //Index "end" is not taken int account
  void averagedNormal(float n[3], int start, int end)
  {
    n[0] = n[1] = n[2] = 0;
    for(int i=start; i<end; i++){
      n[0] += _normal[i][0];
      n[1] += _normal[i][1];
      n[2] += _normal[i][2];
    }
    double len = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    if((float)len != 0)
    {
      n[0] = (float)(n[0]/len);
      n[1] = (float)(n[1]/len);
      n[2] = (float)(n[2]/len);
    }
  }
  
  void bound(float min[], float max[])
  {
    min[0] = min[1] = min[2] = std::numeric_limits<float>::infinity();
    max[0] = max[1] = max[2] = -std::numeric_limits<float>::infinity();
    for(int i=0; i<_pointN; i++)
    {
      float *p = _point[i];
      if(p[0] < min[0])
        min[0] = p[0];
      if(p[0] > max[0])
        max[0] = p[0];
      
      if(p[1] < min[1])
        min[1] = p[1];
      if(p[1] > max[1])
        max[1] = p[1];
      
      if(p[2] < min[2])
        min[2] = p[2];
      if(p[2] > max[2])
        max[2] = p[2];
    }
  }
  
  void bound(float min[], float max[], int start, int end)
  {
    min[0] = min[1] = min[2] = std::numeric_limits<float>::infinity();
    max[0] = max[1] = max[2] = -std::numeric_limits<float>::infinity();
    for(int i=start; i<end; i++){
      float *p = _point[i];
      if(p[0] < min[0])
        min[0] = p[0];
      if(p[0] > max[0])
        max[0] = p[0];
      
      if(p[1] < min[1])
        min[1] = p[1];
      if(p[1] > max[1])
        max[1] = p[1];
      
      if(p[2] < min[2])
        min[2] = p[2];
      if(p[2] > max[2])
        max[2] = p[2];
    }
  }
  
  void shift(float sx, float sy, float sz)
  {
    for(int i=0; i<_pointN; i++)
    {
      _point[i][0] += sx;
      _point[i][1] += sy;
      _point[i][2] += sz;
    }
  }
  
  void scale(float s)
  {
    for(int i=0; i<_pointN; i++)
    {
      _point[i][0] *= s;
      _point[i][1] *= s;
      _point[i][2] *= s;
    }
  }
  
  void rotate(float rx, float ry, float rz)
  {
    rx = (float)(rx*vnl_math::pi/180.0);
    ry = (float)(ry*vnl_math::pi/180.0);
    rz = (float)(rz*vnl_math::pi/180.0);
    if(rx != 0)
    {
      double s = sin(rx);
      double c = cos(rx);
      for(int i=0; i<_pointN; i++)
      {
        float y = _point[i][1];
        float z = _point[i][2];
        
        _point[i][1] = (float)(c*y + s*z);
        _point[i][2] = (float)(-s*y + c*z);
        
        y = _normal[i][1];
        z = _normal[i][2];
        
        _normal[i][1] = (float)(c*y + s*z);
        _normal[i][2] = (float)(-s*y + c*z);
      }
    }
    if(ry != 0)
    {
      double s = sin(ry);
      double c = cos(ry);
      for(int i=0; i<_pointN; i++)
      {
        float z = _point[i][2];
        float x = _point[i][0];
        
        _point[i][2] = (float)(c*z + s*x);
        _point[i][0] = (float)(-s*z + c*x);
        
        z = _normal[i][2];
        x = _normal[i][0];
        
        _normal[i][2] = (float)(c*z + s*x);
        _normal[i][0] = (float)(-s*z + c*x);
      }
    }
    if(rz != 0)
    {
      double s = sin(rz);
      double c = cos(rz);
      for(int i=0; i<_pointN; i++)
      {
        float x = _point[i][0];
        float y = _point[i][1];
        
        _point[i][0] = (float)(c*x + s*y);
        _point[i][1] = (float)(-s*x + c*y);
        
        x = _normal[i][0];
        y = _normal[i][1];
        
        _normal[i][0] = (float)(c*x + s*y);
        _normal[i][1] = (float)(-s*x + c*y);
      }
    }
  }
  
  PointSet* copyPoints()
  {
    PointSet* c = new PointSet();
    c->setPointSize(_pointN);
    for(int i=0; i<_pointN; i++)
    {
      c->setPoint(i, _point[i][0], _point[i][1], _point[i][2]);
      c->setNormal(i, _normal[i][0], _normal[i][1], _normal[i][2]);
    }
    return c;
  }
};

#endif