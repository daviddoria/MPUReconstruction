#include "QuadricO.h"
#include "ImplicitOctTree.h"

#include <vnl/algo/vnl_svd.h>

#define SVD_T_O 0.00001f

QuadricO::QuadricO(PointSet* ps, float R_error, float R_laf, ImplicitOctCell *cell, 
                   int* index_list, int listN, float p_ave[3], float n_ave[3]){
  float c[3];
  cell->cellCenter(c);
  
  float (*point)[3] = ps->_point;
  float (*normal)[3] = ps->_normal;
  
  float o[3], n[3], t1[3], t2[3];
  o[0] = p_ave[0];
  o[1] = p_ave[1];
  o[2] = p_ave[2];
  
  n[0] = n_ave[0];
  n[1] = n_ave[1];
  n[2] = n_ave[2];
  
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

    double totalW = 0;
    for(int i=0; i<listN; i++){
      int in = index_list[i];
      float* p = point[in];
      
      float vx = p[0] - c[0];
      float vy = p[1] - c[1];
      float vz = p[2] - c[2];
      
      double d2 = vx*vx+vy*vy+vz*vz;
      
      float w = (float)cell->weight(sqrt(d2), R_laf);
      totalW += w;
      
      vx = p[0] - o[0];
      vy = p[1] - o[1];
      vz = p[2] - o[2];
      
      float u = t1[0]*vx + t1[1]*vy + t1[2]*vz;
      float v = t2[0]*vx + t2[1]*vy + t2[2]*vz;
      float g = w*(n[0]*vx + n[1]*vy + n[2]*vz);
      
      float tmp[6];
      tmp[0] = w;
      tmp[1] = w*u;
      tmp[2] = w*v;
      tmp[3] = w*u*u;
      tmp[4] = w*u*v;
      tmp[5] = w*v*v;
      
      for(int j=0; j<6; j++){
        for(int k=j; k<6; k++)
          A[j][k] += tmp[j]*tmp[k];
        b[j] += tmp[j]*g;
      }
    }
    /*
    for(i=1; i<7; i++){
      for(int j=i; j<7; j++)
        A[i][j] /= (float)totalW;
      b[i] /= (float)totalW;
    }*/
    
    for(int i=1; i<6; i++)
      for(int j=0; j<i; j++)
        A[i][j] = A[j][i];
    
    vnl_svd< float > svd( A );
    svd.zero_out_relative( SVD_T_O );

    vnl_vector< float > x = svd.solve( b );
    

  float c0 = x[0];
  float cu = x[1];
  float cv = x[2];
  float cuu = x[3];
  float cuv = x[4];
  float cvv = x[5];
  
  //convert into world coordinates (u, v, w) -> (x, y, z)
  _cxx = -(cuu*t1[0]*t1[0] + cvv*t2[0]*t2[0] + cuv*t1[0]*t2[0]);
  _cyy = -(cuu*t1[1]*t1[1] + cvv*t2[1]*t2[1] + cuv*t1[1]*t2[1]);
  _czz = -(cuu*t1[2]*t1[2] + cvv*t2[2]*t2[2] + cuv*t1[2]*t2[2]);
  
  _cxy = -(2.0f*(cuu*t1[0]*t1[1] + cvv*t2[0]*t2[1]) + cuv*(t1[0]*t2[1] + t1[1]*t2[0]));
  _cyz = -(2.0f*(cuu*t1[1]*t1[2] + cvv*t2[1]*t2[2]) + cuv*(t1[1]*t2[2] + t1[2]*t2[1]));
  _czx = -(2.0f*(cuu*t1[2]*t1[0] + cvv*t2[2]*t2[0]) + cuv*(t1[2]*t2[0] + t1[0]*t2[2]));
  
  float cx = n[0] - cu*t1[0] - cv*t2[0];
  float cy = n[1] - cu*t1[1] - cv*t2[1];
  float cz = n[2] - cu*t1[2] - cv*t2[2];
  
  _cx = cx - _cxy*o[1] - _czx*o[2] - 2.0f*_cxx*o[0];
  _cy = cy - _cyz*o[2] - _cxy*o[0] - 2.0f*_cyy*o[1];
  _cz = cz - _czx*o[0] - _cyz*o[1] - 2.0f*_czz*o[2];
  
  _c0 = -c0 - cx*o[0] - cy*o[1] - cz*o[2] 
    + _cxy*o[0]*o[1] + _cyz*o[1]*o[2] + _czx*o[2]*o[0]
      + _cxx*o[0]*o[0] + _cyy*o[1]*o[1] + _czz*o[2]*o[2];
  
  if(R_error == 0){
    cell->_error = 0;
    return;
  }
  
  float error = 0;
  //max error
  bool* bound = ps->_bound;
  float R2 = R_error*R_error;
  for(int i=0; i<listN; i++){
    int j = index_list[i];
    
    //if(bound[j])
      //continue;
    
    float* p = point[j];
    
    float d = (p[0]-c[0])*(p[0]-c[0]) + (p[1]-c[1])*(p[1]-c[1]) + (p[2]-c[2])*(p[2]-c[2]);
    if(d > R2)
      continue;
    
    float f = (float)fabs(value(p[0], p[1], p[2]));
    float g[3];
    gradient(g, p[0], p[1], p[2]);
    float l = (float)sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);
    if(l != 0)
      f /= l;
    
    if(f > error)
      error = f;
  }
  cell->_error = error;	
}
