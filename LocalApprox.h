#ifndef LOCALAPPROX_H 
#define LOCALAPPROX_H 

#include <stdio.h>
#include <math.h>

class LocalApprox{
public:  
  virtual float value(float x, float y, float z)=0;
  
  virtual void gradient(float g[3], float x, float y, float z)=0;
  
  virtual float valueAndGradient(float g[3], float x, float y, float z)=0;
  
};

#endif