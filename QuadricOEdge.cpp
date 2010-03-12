#include "QuadricOEdge.h"
#include "ImplicitOctTree.h"

#include <vnl/algo/vnl_svd.h>
#include <vnl/vnl_math.h>

#define SVD_T_E 0.001f

float QuadricOEdge::value(float x, float y, float z){
  float f1 = _c10 + 
    (_c1x + _c1xx*x + _c1xy*y)*x +
      (_c1y + _c1yy*y + _c1yz*z)*y +
        (_c1z + _c1zz*z + _c1zx*x)*z;
  
  float f2 = _c20 + 
    (_c2x + _c2xx*x + _c2xy*y)*x +
      (_c2y + _c2yy*y + _c2yz*z)*y +
        (_c2z + _c2zz*z + _c2zx*x)*z;
  
  if(_is_convex){
    if(f1 > f2)
      return f1;
    else
      return f2;
  }
  else{
    if(f1 < f2)
      return f1;
    else
      return f2;
  }
}

QuadricOEdge::QuadricOEdge(PointSet* ps, float R_error, float R_laf, ImplicitOctCell *cell,
                           int* index_list, int listN, int index1, int index2){
  float c[3];
  cell->cellCenter(c);
  
  float (*point)[3] = ps->_point;
  float (*normal)[3] = ps->_normal;
  
  //separate two point set
  float *n1 = normal[index1];
  float *n2 = normal[index2];
  int i = 0;
  int j = listN-1;
  while(i <= j){
    float*m = normal[index_list[i]];
    if(n1[0]*m[0] + n1[1]*m[1] + n1[2]*m[2] > n2[0]*m[0] + n2[1]*m[1] + n2[2]*m[2])
      i++;
    else{
      int tmp = index_list[i];
      index_list[i] = index_list[j];
      index_list[j] = tmp;
      j--;
    }  
  }
  int split = i;
  
  //compute functions
  float c1[3], c2[3];
  float m1[3], m2[3];
  computeCoefficient(ps, R_laf, cell, index_list, 0, split, 0, c1, m1);
  computeCoefficient(ps, R_laf, cell, index_list, split, listN, 1, c2, m2);
  
  //check convex or concave
  float vx = c1[0] - c2[0];
  float vy = c1[1] - c2[1];
  float vz = c1[2] - c2[2];
  double len = sqrt(vx*vx + vy*vy + vz*vz);
  
  n1 = m1;
  n2 = m2;
  
  double dot = (vx*n1[0] + vy*n1[1] + vz*n1[2])/len;
  if(dot > 1.0)
    dot = 1.0;
  else if(dot < -1.0)
    dot = -1.0;
  double angle1 = acos(dot);
  
  dot = -(vx*n2[0] + vy*n2[1] + vz*n2[2])/len;
  if(dot > 1.0)
    dot = 1.0;
  else if(dot < -1.0)
    dot = -1.0;
  double angle2 = acos(dot);
  
  _is_convex = (angle1 + angle2 < vnl_math::pi );
  
  /*
  if(R_error == 0){
    cell->_error = 0;
    return;
  }*/
  
  //max error
  float error = 0;
  
  float R2 = R_error*R_error;
//  bool* bound = ps->_bound;
  for(i=0; i<listN; i++){
    j = index_list[i];
    //if(bound[j])
      //continue;
    
    float* p = point[j];
    float vx = p[0] - c[0];
    float vy = p[1] - c[1];
    float vz = p[2] - c[2];
    
    if(R2 < vx*vx + vy*vy + vz*vz)
      continue;
    
    float g[3];
    float f = (float)fabs(valueAndGradient(g, p[0], p[1], p[2]));
    float l = (float)sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);
    if(l != 0)
      f /= l;
    
    if(f > error)
      error = f;
  }
  cell->_error = error;
}

void QuadricOEdge::computeCoefficient(PointSet* ps, float R, ImplicitOctCell *cell,
                                      int* index_list, int start, int end, int f_index, float o[3], float n[3]){
  float (*point)[3] = ps->_point;
  float (*normal)[3] = ps->_normal;
  
  float c[3];
  cell->cellCenter(c);
  
  float t1[3], t2[3]; 
  
  n[0] = n[1] = n[2] = 0;
  o[0] = o[1] = o[2] = 0;
  double totalW = 0;
  for(int i=start; i<end; i++){
    int j = index_list[i];
    float* p = point[j];
    float* m = normal[j];
    
    float vx = p[0] - c[0];
    float vy = p[1] - c[1];
    float vz = p[2] - c[2];
    
    double w = cell->weight(sqrt(vx*vx+vy*vy+vz*vz), R);
    totalW += w;
    
    o[0] += (float)(w*p[0]);
    o[1] += (float)(w*p[1]);
    o[2] += (float)(w*p[2]);
    
    n[0] += (float)(w*m[0]);
    n[1] += (float)(w*m[1]);
    n[2] += (float)(w*m[2]);
  }
  double len = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  if((float)len != 0){
    n[0] = (float)(n[0]/len);
    n[1] = (float)(n[1]/len);
    n[2] = (float)(n[2]/len);
  }
  o[0] = (float)(o[0]/totalW);
  o[1] = (float)(o[1]/totalW);
  o[2] = (float)(o[2]/totalW);
  
  if(end - start < 6 || (float)len == 0){
    if(f_index == 0){
      _c1xx = _c1yy = _c1zz = 0;
      _c1xy = _c1yz = _c1zx = 0;
      _c1x = n[0];
      _c1y = n[1];
      _c1z = n[2];
      _c10 = -(_c1x*o[0] + _c1y*o[1] + _c1z*o[2]);
    }
    else{
      _c2xx = _c2yy = _c2zz = 0;
      _c2xy = _c2yz = _c2zx = 0;
      _c2x = n[0];
      _c2y = n[1];
      _c2z = n[2];
      _c20 = -(_c2x*o[0] + _c2y*o[1] + _c2z*o[2]);
    }
    return;
  }
  
  //compute base
  if(fabs(n[0]) < fabs(n[1])){
    double l = sqrt(n[1]*n[1] + n[2]*n[2]);
    t1[0] = 0;
    t1[1] = -(float)(n[2]/l);
    t1[2] = (float)(n[1]/l);
  }
  else{
    double l = sqrt(n[0]*n[0] + n[2]*n[2]);
    t1[0] = (float)(n[2]/l);
    t1[1] = 0;
    t1[2] = -(float)(n[0]/l);
  }
  t2[0] = n[1]*t1[2] - n[2]*t1[1];
  t2[1] = n[2]*t1[0] - n[0]*t1[2];
  t2[2] = n[0]*t1[1] - n[1]*t1[0];
  
  //least square fitting
  vnl_matrix< float > A( 6, 6, 0. );
  vnl_vector< float > b( 6, 0. );
  
  for(int i=start; i<end; i++){
    int in = index_list[i];
    float* p = point[in];
    
    float vx = p[0] - c[0];
    float vy = p[1] - c[1];
    float vz = p[2] - c[2];
    
    float w = (float)cell->weight(sqrt(vx*vx+vy*vy+vz*vz), R);
    
    vx = p[0] - o[0];
    vy = p[1] - o[1];
    vz = p[2] - o[2];
    
    float u = t1[0]*vx + t1[1]*vy + t1[2]*vz;
    float v = t2[0]*vx + t2[1]*vy + t2[2]*vz;
    float g = n[0]*vx + n[1]*vy + n[2]*vz;
    
    float tmp[6];
    tmp[0] = w;
    tmp[1] = w*u;
    tmp[2] = w*v;
    tmp[3] = w*u*u;
    tmp[4] = w*u*v;
    tmp[5] = w*v*v;
    
    for(int j=0; j<6; j++){
      for(int k=0; k<6; k++)
        A[j][k] += tmp[j]*tmp[k];
      b[j] += w*tmp[j]*g;
    }
  }
  /*
  for(int j=0; j<6; j++){
    for(int k=j; k<6; k++)
      A[j+1][k+1] /= (float)totalW;
    b[j+1] /= (float)totalW;
  }*/
  
  for(int i=1; i<6; i++)
    for(int j=0; j<i; j++)
      A[i][j] = A[j][i];
  
  vnl_svd< float > svd( A );

  vnl_vector< float > w = svd.W().diagonal();
  
  float wmax=0.0f;
  for (int k=0;k<6;k++)
    if (fabs(w[k]) > wmax) wmax=(float)fabs(w[k]);
  
  if(wmax < 0.00000000001f){
    if(f_index == 0){
      _c1xx = _c1yy = _c1zz = 0;
      _c1xy = _c1yz = _c1zx = 0;
      _c1x = n[0];
      _c1y = n[1];
      _c1z = n[2];
      _c10 = -(_c1x*o[0] + _c1y*o[1] + _c1z*o[2]);
    }
    else{
      _c2xx = _c2yy = _c2zz = 0;
      _c2xy = _c2yz = _c2zx = 0;
      _c2x = n[0];
      _c2y = n[1];
      _c2z = n[2];
      _c20 = -(_c2x*o[0] + _c2y*o[1] + _c2z*o[2]);
    }
    return;
  }

  svd.zero_out_relative( SVD_T_E );
  vnl_vector< float > x = svd.solve( b );
  
  float c0 = x[0];
  float cu = x[1];
  float cv = x[2];
  float cuu = x[3];
  float cuv = x[4];
  float cvv = x[5];
  
  //convert into world coordinates (u, x, w) -> (x, y, z)
  float cx = n[0] - cu*t1[0] - cv*t2[0];
  float cy = n[1] - cu*t1[1] - cv*t2[1];
  float cz = n[2] - cu*t1[2] - cv*t2[2];
  
  if(f_index == 0){
    _c1xx = -(cuu*t1[0]*t1[0] + cvv*t2[0]*t2[0] + cuv*t1[0]*t2[0]);
    _c1yy = -(cuu*t1[1]*t1[1] + cvv*t2[1]*t2[1] + cuv*t1[1]*t2[1]);
    _c1zz = -(cuu*t1[2]*t1[2] + cvv*t2[2]*t2[2] + cuv*t1[2]*t2[2]);
    
    _c1xy = -(2.0f*(cuu*t1[0]*t1[1] + cvv*t2[0]*t2[1]) + cuv*(t1[0]*t2[1] + t1[1]*t2[0]));
    _c1yz = -(2.0f*(cuu*t1[1]*t1[2] + cvv*t2[1]*t2[2]) + cuv*(t1[1]*t2[2] + t1[2]*t2[1]));
    _c1zx = -(2.0f*(cuu*t1[2]*t1[0] + cvv*t2[2]*t2[0]) + cuv*(t1[2]*t2[0] + t1[0]*t2[2]));
    
    _c1x = cx - _c1xy*o[1] - _c1zx*o[2] - 2.0f*_c1xx*o[0];
    _c1y = cy - _c1yz*o[2] - _c1xy*o[0] - 2.0f*_c1yy*o[1];
    _c1z = cz - _c1zx*o[0] - _c1yz*o[1] - 2.0f*_c1zz*o[2];
    
    _c10 = -c0 - cx*o[0] - cy*o[1] - cz*o[2] 
      + _c1xy*o[0]*o[1] + _c1yz*o[1]*o[2] + _c1zx*o[2]*o[0]
        + _c1xx*o[0]*o[0] + _c1yy*o[1]*o[1] + _c1zz*o[2]*o[2];
  }
  else{
    _c2xx = -(cuu*t1[0]*t1[0] + cvv*t2[0]*t2[0] + cuv*t1[0]*t2[0]);
    _c2yy = -(cuu*t1[1]*t1[1] + cvv*t2[1]*t2[1] + cuv*t1[1]*t2[1]);
    _c2zz = -(cuu*t1[2]*t1[2] + cvv*t2[2]*t2[2] + cuv*t1[2]*t2[2]);
    
    _c2xy = -(2.0f*(cuu*t1[0]*t1[1] + cvv*t2[0]*t2[1]) + cuv*(t1[0]*t2[1] + t1[1]*t2[0]));
    _c2yz = -(2.0f*(cuu*t1[1]*t1[2] + cvv*t2[1]*t2[2]) + cuv*(t1[1]*t2[2] + t1[2]*t2[1]));
    _c2zx = -(2.0f*(cuu*t1[2]*t1[0] + cvv*t2[2]*t2[0]) + cuv*(t1[2]*t2[0] + t1[0]*t2[2]));
    
    _c2x = cx - _c2xy*o[1] - _c2zx*o[2] - 2.0f*_c2xx*o[0];
    _c2y = cy - _c2yz*o[2] - _c2xy*o[0] - 2.0f*_c2yy*o[1];
    _c2z = cz - _c2zx*o[0] - _c2yz*o[1] - 2.0f*_c2zz*o[2];
    
    _c20 = -c0 - cx*o[0] - cy*o[1] - cz*o[2] 
      + _c2xy*o[0]*o[1] + _c2yz*o[1]*o[2] + _c2zx*o[2]*o[0]
        + _c2xx*o[0]*o[0] + _c2yy*o[1]*o[1] + _c2zz*o[2]*o[2];
  }
}