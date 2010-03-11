#include "QuadricOCorner.h"
#include "ImplicitOctTree.h"

#include <vnl/algo/vnl_svd.h>
#include <vnl/vnl_math.h>

#define SVD_T_C 0.001f

QuadricOCorner::QuadricOCorner(PointSet* ps, float R_error, float R_laf, ImplicitOctCell *cell,
                               int* index_list, int listN, int index1, int index2, int index3, bool &OK, float edgeT){
  float c[3];
  cell->cellCenter(c);
  
  float (*point)[3] = ps->_point;
  float (*normal)[3] = ps->_normal;
  
  
  //separate point set
  int *id = new int[listN];
  for(int i=0; i<listN; i++)
    id[i] = -1;
  
  bool flag = true;
  
  float *n1 = normal[index1];
  float *n2 = normal[index2];
  float n3[3]; // = normal[index3];
  n3[0] = n1[1]*n2[2] - n1[2]*n2[1];
  n3[1] = n1[2]*n2[0] - n1[0]*n2[2];
  n3[2] = n1[0]*n2[1] - n1[1]*n2[0];
  float len = (float)sqrt(n3[0]*n3[0] + n3[1]*n3[1] + n3[2]*n3[2]);
  if(len != 0){
    n3[0] /= len;
    n3[1] /= len;
    n3[2] /= len;
  }
  
  
  for(int i=0; i<listN; i++){
    float*m = normal[index_list[i]];
    float dot1 = n1[0]*m[0] + n1[1]*m[1] + n1[2]*m[2];
    float dot2 = n2[0]*m[0] + n2[1]*m[1] + n2[2]*m[2];
    float dot3 = (float)fabs(n3[0]*m[0] + n3[1]*m[1] + n3[2]*m[2]);
    
    if(dot1 > dot2){
      if(fabs(dot1) > dot3)
        id[i] = 0;
      else
        id[i] = 2;
    }
    else{
      if(fabs(dot2) > dot3)
        id[i] = 1;
      else
        id[i] = 2;
    }
  }
  
  //edge check again in the set id = 2
  int mini, minj;
  float min = 1;
  for(int i=0; i<listN; i++){
    if(id[i] != 2)
      continue;
    float* ni = normal[index_list[i]];
    for(int j=i+1; j<listN; j++){
      if(id[j] != 2)
        continue;
      float* nj = normal[index_list[j]];
      float v = ni[0]*nj[0] + ni[1]*nj[1] + ni[2]*nj[2];
      if(v < min){
        min = v;
        mini = index_list[i];
        minj = index_list[j];
      }
    }
  }
  if(min < edgeT){
    _deg = 4;
    n1 = normal[mini];
    n2 = normal[minj];
    for(int i=0; i<listN; i++){
      if(id[i] != 2)
        continue;
      float*m = normal[index_list[i]];
      if(n1[0]*m[0] + n1[1]*m[1] + n1[2]*m[2] < n2[0]*m[0] + n2[1]*m[1] + n2[2]*m[2])
        id[i] = 3;
    }
  }
  else
    _deg = 3;
  
  int* split = new int[_deg+1];
  split[0] = 0;
  for(int i=0; i<_deg-1; i++){
    int count = split[i];
    for(int j=split[i]; j<listN; j++){
      if(id[j] == i){
        int tmp = index_list[j];
        index_list[j] = index_list[count];
        index_list[count] = tmp;
        
        tmp = id[j];
        id[j] = id[count];
        id[count] = tmp;
        
        count++;
      }
    }
    split[i+1] = count;
  }
  split[_deg] = listN;
  delete[] id;
  
  //compute functions
  _c0 = new float[_deg];
  _cx = new float[_deg]; _cy = new float[_deg]; _cz = new float[_deg];
  _cxx = new float[_deg]; _cyy = new float[_deg]; _czz = new float[_deg];
  _cxy = new float[_deg]; _cyz = new float[_deg]; _czx = new float[_deg];
  
  float (*o)[3] = new float[_deg][3];
  float (*m)[3] = new float[_deg][3];
  for(int i=0; i<_deg; i++)
    computeCoefficient(ps, R_laf, cell, index_list, split[i], split[i+1], i, o[i], m[i]);
  
  delete[] split;
  
  //check convex or concave
  _op = new bool[_deg-1];
  if(_deg == 3){
    OK = true;
    bool op_tmp[3];
    for(int i=0; i<3; i++){
      int j = (i+1)%3;
      float vx = o[i][0] - o[j][0];
      float vy = o[i][1] - o[j][1];
      float vz = o[i][2] - o[j][2];
      double len = sqrt(vx*vx + vy*vy + vz*vz);
      
      double dot = (vx*m[i][0] + vy*m[i][1] + vz*m[i][2])/len;
      if(dot > 1.0)
        dot = 1.0;
      else if(dot < -1.0)
        dot = -1.0;
      double angle1 = acos(dot);
      
      dot = -(vx*m[j][0] + vy*m[j][1] + vz*m[j][2])/len;
      if(dot > 1.0)
        dot = 1.0;
      else if(dot < -1.0)
        dot = -1.0;
      double angle2 = acos(dot);
      
      op_tmp[i] = (angle1 + angle2 < vnl_math::pi );
    }
    
    if(op_tmp[0] && op_tmp[1] && op_tmp[2])
      _op[0] = _op[1] = true;
    else if(!op_tmp[0] && !op_tmp[1] && !op_tmp[2])
      _op[0] = _op[1] = false;
    
    //saddle cases
    else if(op_tmp[0] && !op_tmp[1] && !op_tmp[2]){
      _op[0] = true;
      _op[1] = false;
    }
    else if(!op_tmp[0] && op_tmp[1] && op_tmp[2]){
      _op[0] = false;
      _op[1] = true;
    }
    else if(!op_tmp[0] && !op_tmp[1] && op_tmp[2]){
      swapCoefficient(1, 2);
      _op[0] = true;
      _op[1] = false;
    }
    else if(!op_tmp[0] && op_tmp[1] && !op_tmp[2]){
      swapCoefficient(0, 2);
      _op[0] = true;
      _op[1] = false;
    }
    else if(op_tmp[0] && op_tmp[1] && !op_tmp[2]){
      swapCoefficient(1, 2);
      _op[0] = false;
      _op[1] = true;
    }
    else if(op_tmp[0] && !op_tmp[1] && op_tmp[2]){
      swapCoefficient(0, 2);
      _op[0] = false;
      _op[1] = true;
    }
  }
  //degree = 4
  //If saddle corner, failed.
  else{
    bool op_tmp[6];
    int k = 0;
    for(int i=0; i<3; i++){
      for(int j=i+1; j<4; j++){
        float vx = o[i][0] - o[j][0];
        float vy = o[i][1] - o[j][1];
        float vz = o[i][2] - o[j][2];
        double len = sqrt(vx*vx + vy*vy + vz*vz);
        
        double dot = (vx*m[i][0] + vy*m[i][1] + vz*m[i][2])/len;
        if(dot > 1.0)
          dot = 1.0;
        else if(dot < -1.0)
          dot = -1.0;
        double angle1 = acos(dot);
        
        dot = -(vx*m[j][0] + vy*m[j][1] + vz*m[j][2])/len;
        if(dot > 1.0)
          dot = 1.0;
        else if(dot < -1.0)
          dot = -1.0;
        double angle2 = acos(dot);
        
        op_tmp[k++] = (angle1 + angle2 < vnl_math::pi);
      }
    }
    OK = true;
    for(int i=1; i<k; i++){
      if(op_tmp[0] != op_tmp[i]){
        OK = false;
        break;
      }
    }
    
    if(OK){
      _op[0] = _op[1] = _op[2] = op_tmp[0];
    }
  }
  
  delete[] o;
  delete[] m;
  
  if(!OK)
    return;
  
  if(R_error == 0){
    cell->_error = 0;
    return;
  }
  
  //max error
  float error = 0;
  
  float R2 = R_error*R_error;
  bool* bound = ps->_bound;
  for(int i=0; i<listN; i++){
    int j = index_list[i];
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

void QuadricOCorner::computeCoefficient(PointSet* ps, float R, ImplicitOctCell *cell,
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
    _cxx[f_index] = _cyy[f_index] = _czz[f_index] = 0;
    _cxy[f_index] = _cyz[f_index] = _czx[f_index] = 0;
    _cx[f_index] = n[0];
    _cy[f_index] = n[1];
    _cz[f_index] = n[2];
    _c0[f_index] = -(n[0]*o[0] + n[1]*o[1] + n[2]*o[2]);
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
      for(int k=j; k<6; k++)
        A[j][k] += tmp[j]*tmp[k];
      b[j] += w*tmp[j]*g;
    }
  }
  for(int i=1; i<6; i++)
    for(int j=0; j<i; j++)
      A[i][j] = A[j][i];
  
  /*
  for(int j=0; j<6; j++){
      for(int k=0; k<6; k++)
        A[j+1][k+1] /= (float)totalW;
      b[j+1] /= (float)totalW;
  }*/
  
  vnl_svd< float > svd( A );
  svd.zero_out_relative( SVD_T_C );

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
  
  _cxx[f_index] = -(cuu*t1[0]*t1[0] + cvv*t2[0]*t2[0] + cuv*t1[0]*t2[0]);
  _cyy[f_index] = -(cuu*t1[1]*t1[1] + cvv*t2[1]*t2[1] + cuv*t1[1]*t2[1]);
  _czz[f_index] = -(cuu*t1[2]*t1[2] + cvv*t2[2]*t2[2] + cuv*t1[2]*t2[2]);
  
  _cxy[f_index] = -(2.0f*(cuu*t1[0]*t1[1] + cvv*t2[0]*t2[1]) + cuv*(t1[0]*t2[1] + t1[1]*t2[0]));
  _cyz[f_index] = -(2.0f*(cuu*t1[1]*t1[2] + cvv*t2[1]*t2[2]) + cuv*(t1[1]*t2[2] + t1[2]*t2[1]));
  _czx[f_index] = -(2.0f*(cuu*t1[2]*t1[0] + cvv*t2[2]*t2[0]) + cuv*(t1[2]*t2[0] + t1[0]*t2[2]));
  
  _cx[f_index] = cx - _cxy[f_index]*o[1] - _czx[f_index]*o[2] - 2.0f*_cxx[f_index]*o[0];
  _cy[f_index] = cy - _cyz[f_index]*o[2] - _cxy[f_index]*o[0] - 2.0f*_cyy[f_index]*o[1];
  _cz[f_index] = cz - _czx[f_index]*o[0] - _cyz[f_index]*o[1] - 2.0f*_czz[f_index]*o[2];
  
  _c0[f_index] = -c0 - cx*o[0] - cy*o[1] - cz*o[2] 
    + _cxy[f_index]*o[0]*o[1] + _cyz[f_index]*o[1]*o[2] + _czx[f_index]*o[2]*o[0]
      + _cxx[f_index]*o[0]*o[0] + _cyy[f_index]*o[1]*o[1] + _czz[f_index]*o[2]*o[2];
}
