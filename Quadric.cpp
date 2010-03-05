#include "Quadric.h"
#include "ImplicitOctTree.h"
#include "SVD.h"

#define SVD_T 0.0000001f
#define NEAR 6
//#define N_SMALL_Q 20

Quadric::Quadric(PointSet* ps, float R_error, float R_laf, ImplicitOctCell *cell, 
                 int* index_list, int listN, float p_ave[3]){
  //if(listN > N_SMALL_Q)
    computePolySVD(ps, R_laf, cell, index_list, listN, p_ave);
  //else
    //computePolySVD2(ps, R_laf, cell, index_list, listN, p_ave);
  
  if(R_error != 0)
    cell->_error = computeMaxError(ps, R_error, cell, index_list, listN);
  else
    cell->_error = 0;
}


//This function is called when the number of points to be fitted is not large.
void Quadric::computePolySVD2(PointSet* ps, float R, ImplicitOctCell *cell, 
                              int* index_list, int listN, float p_ave[3]){
  float (*point)[3] = ps->_point;
  float (*normal)[3] = ps->_normal;
  
  float o[3];
  o[0] = p_ave[0];
  o[1] = p_ave[1];
  o[2] = p_ave[2];
  
  float c[3];
  cell->cellCenter(c);
  
  //additional points
  int adN = 9;
  float (*ad_point)[3] = new float[adN][3];
  int (*ad_i)[NEAR] = new int[adN][NEAR];
  float (*ad_d)[NEAR] = new float[adN][NEAR];
  for(int i=0; i<adN; i++){
    for(int j=0; j<NEAR; j++){
      ad_i[i][j] = -1;
      ad_d[i][j] = std::numeric_limits<float>::infinity();
    }
  }
  
  //For SVD
  float** A= new float*[listN+adN+1];
  float* b = new float[listN+adN+1];
  
  float size = cell->_size;
  
  ad_point[0][0] = c[0] - size;
  ad_point[0][1] = c[1] - size;
  ad_point[0][2] = c[2] - size;
  
  ad_point[1][0] = c[0] + size;
  ad_point[1][1] = c[1] - size;
  ad_point[1][2] = c[2] - size;
  
  ad_point[2][0] = c[0] - size;
  ad_point[2][1] = c[1] + size;
  ad_point[2][2] = c[2] - size;
  
  ad_point[3][0] = c[0] + size;
  ad_point[3][1] = c[1] + size;
  ad_point[3][2] = c[2] - size;
  
  
  ad_point[4][0] = c[0] - size;
  ad_point[4][1] = c[1] - size;
  ad_point[4][2] = c[2] + size;
  
  ad_point[5][0] = c[0] + size;
  ad_point[5][1] = c[1] - size;
  ad_point[5][2] = c[2] + size;
  
  ad_point[6][0] = c[0] - size;
  ad_point[6][1] = c[1] + size;
  ad_point[6][2] = c[2] + size;
  
  ad_point[7][0] = c[0] + size;
  ad_point[7][1] = c[1] + size;
  ad_point[7][2] = c[2] + size;
  
  ad_point[8][0] = c[0];
  ad_point[8][1] = c[1];
  ad_point[8][2] = c[2];
  
  double totalW = 0;
  for(int i=0; i<listN; i++){
    int in = index_list[i];
    float* p = point[in];
    float* n = normal[in];
    
    float vx = p[0] - c[0];
    float vy = p[1] - c[1];
    float vz = p[2] - c[2];
    float w = (float)cell->weight(sqrt(vx*vx+vy*vy+vz*vz), R);
    totalW += w;
    
    vx = p[0] - o[0];
    vy = p[1] - o[1];
    vz = p[2] - o[2];
    
    float* Ai = A[i+1] = new float[11];
    Ai[1] = w;
    Ai[2] = w*vx;
    Ai[3] = w*vy;
    Ai[4] = w*vz;
    Ai[5] = w*vx*vx;
    Ai[6] = w*vy*vy;
    Ai[7] = w*vz*vz;
    Ai[8] = w*vx*vy;
    Ai[9] = w*vy*vz;
    Ai[10] = w*vz*vx;
    
    b[i+1] = 0;
    
    //near points
    for(int j=0; j<adN; j++){
      float* q = ad_point[j];
      float d = (p[0]-q[0])*(p[0]-q[0]) + 
        (p[1]-q[1])*(p[1]-q[1]) +
          (p[2]-q[2])*(p[2]-q[2]);
      int insert = -1;
      for(int k=0; k<NEAR; k++){
        if(ad_d[j][k] > d)
          insert = k;
        else
          break;
      }
      if(insert < 0)
        continue;
      for(int k=0; k<insert; k++){
        ad_d[j][k] = ad_d[j][k+1];
        ad_i[j][k] = ad_i[j][k+1];
      }
      ad_d[j][insert] = d;
      ad_i[j][insert] = in;
    }
  }
  
  for(int i=1; i<listN+1; i++){
    for(int j=1; j<11; j++)
      A[i][j] /= (float)totalW;
  }
  
  //Extra points
  int count = 0;
  for(int i=0; i<adN; i++){
    float* p = ad_point[i];
    
    //distance (inner product with normal)
    double v = 0;
    for(int j=0; j<NEAR; j++){
      int in = ad_i[i][j];
      float *q = point[in];
      float *n = normal[in];
      ad_d[i][j]= n[0]*(p[0]-q[0]) + n[1]*(p[1]-q[1]) + n[2]*(p[2]-q[2]);
      v += ad_d[i][j];
    }
    v /= NEAR;
    
    //sign check
    bool flag = true;
    for(int j=1; j<NEAR; j++){
      if(ad_d[i][0]*ad_d[i][j] <= 0){
        flag = false;
        break;
      }
    }
    if(!flag)
      continue;
    
    count++;
    float* Ai = A[count+listN] = new float[11];
    
    float w = 1.0f/adN;
    
    float vx = p[0] - o[0];
    float vy = p[1] - o[1];
    float vz = p[2] - o[2];
    
    Ai[1] = w;
    Ai[2] = w*vx;
    Ai[3] = w*vy;
    Ai[4] = w*vz;
    Ai[5] = w*vx*vx;
    Ai[6] = w*vy*vy;
    Ai[7] = w*vz*vz;
    Ai[8] = w*vx*vy;
    Ai[9] = w*vy*vz;
    Ai[10] = w*vz*vx;
    
    b[count+listN] = (float)v*w;
  }
  delete[] ad_d;
  delete[] ad_point;
  delete[] ad_i;
  
  float w[11];
  float **v = new float*[11];
  for(int i=1; i<11; i++)
    v[i] = new float[11];
  SVD::svdcmp(A, listN+count-1, 10, w, v);
  
  float wmax=0.0f;
  for (int k=1;k<11;k++)
    if (fabs(w[k]) > wmax) wmax=(float)fabs(w[k]);
  
  if(wmax < 0.000000000001f || count == 0){
    for(int i=1; i<listN+count; i++)
      delete[] A[i];
    delete[] A;
    
    for(int i=1; i<11; i++)
      delete[] v[i];
    delete[] v;
    
    delete[] b;
    
    return;
  }
  
  float wmin=wmax*SVD_T;
  for (int k=1;k<11;k++){
    if (fabs(w[k]) < wmin) 
      w[k]=0.0;
  }
  
  float x[11];
  SVD::svbksb(A, w, v, listN+count-1, 10, b, x);
  
  _cxx = x[5];
  _cyy = x[6];
  _czz = x[7];
  
  _cxy = x[8];
  _cyz = x[9];
  _czx = x[10];
  
  _cx = x[2] - _cxy*o[1] - _czx*o[2] - 2.0f*_cxx*o[0];
  _cy = x[3] - _cyz*o[2] - _cxy*o[0] - 2.0f*_cyy*o[1];
  _cz = x[4] - _czx*o[0] - _cyz*o[1] - 2.0f*_czz*o[2];
  
  _c0 = x[1] - x[2]*o[0] - x[3]*o[1] - x[4]*o[2] 
        + _cxy*o[0]*o[1] + _cyz*o[1]*o[2] + _czx*o[2]*o[0]
        + _cxx*o[0]*o[0] + _cyy*o[1]*o[1] + _czz*o[2]*o[2];
  
  for(int i=1; i<listN+count; i++)
    delete[] A[i];
  delete[] A;
  
  
  for(int i=1; i<11; i++)
    delete[] v[i];
  delete[] v;
  
  delete[] b;
}

void Quadric::computePolySVD(PointSet* ps, float R, ImplicitOctCell *cell, 
                             int* index_list, int listN, float p_ave[3]){
  float (*point)[3] = ps->_point;
  float (*normal)[3] = ps->_normal;
  
  float o[3];
  o[0] = p_ave[0];
  o[1] = p_ave[1];
  o[2] = p_ave[2];
  
  float c[3];
  cell->cellCenter(c);
  
  //For SVD
  float** A= new float*[11];
  for(int i=1; i<11; i++){
    A[i] = new float[11];
    A[i][1] = A[i][2] = A[i][3] = 
      A[i][4] = A[i][5] = A[i][6] =
        A[i][7] = A[i][8] = A[i][9] = A[i][10] = 0;
  }
  float b[11];
  b[1] = b[2] = b[3] = b[4] = b[5] = b[6] = b[7] = b[8] = b[9] = b[10] = 0;
  
  //additional points
  int adN = 9;
  float (*ad_point)[3] = new float[adN][3];
  int (*ad_i)[NEAR] = new int[adN][NEAR];
  float (*ad_d)[NEAR] = new float[adN][NEAR];
  for(int i=0; i<adN; i++){
    for(int j=0; j<NEAR; j++){
      ad_i[i][j] = -1;
      ad_d[i][j] = std::numeric_limits<float>::infinity();
    }
  }
  
  float size = cell->_size; //R/(float)sqrt(3.0); //cell->_size;
  
  ad_point[0][0] = c[0] - size;
  ad_point[0][1] = c[1] - size;
  ad_point[0][2] = c[2] - size;
  
  ad_point[1][0] = c[0] + size;
  ad_point[1][1] = c[1] - size;
  ad_point[1][2] = c[2] - size;
  
  ad_point[2][0] = c[0] - size;
  ad_point[2][1] = c[1] + size;
  ad_point[2][2] = c[2] - size;
  
  ad_point[3][0] = c[0] + size;
  ad_point[3][1] = c[1] + size;
  ad_point[3][2] = c[2] - size;
  
  
  ad_point[4][0] = c[0] - size;
  ad_point[4][1] = c[1] - size;
  ad_point[4][2] = c[2] + size;
  
  ad_point[5][0] = c[0] + size;
  ad_point[5][1] = c[1] - size;
  ad_point[5][2] = c[2] + size;
  
  ad_point[6][0] = c[0] - size;
  ad_point[6][1] = c[1] + size;
  ad_point[6][2] = c[2] + size;
  
  ad_point[7][0] = c[0] + size;
  ad_point[7][1] = c[1] + size;
  ad_point[7][2] = c[2] + size;
  
  ad_point[8][0] = c[0];
  ad_point[8][1] = c[1];
  ad_point[8][2] = c[2];
  
  double totalW = 0;
  for(int i=0; i<listN; i++){
    int in = index_list[i];
    float* p = point[in];
    float* n = normal[in];
    
    float vx = p[0] - c[0];
    float vy = p[1] - c[1];
    float vz = p[2] - c[2];
    float w = (float)cell->weight(sqrt(vx*vx+vy*vy+vz*vz), R);
    totalW += w;
    
    vx = p[0] - o[0];
    vy = p[1] - o[1];
    vz = p[2] - o[2];
    
    float x[10];
    x[0] = w;
    x[1] = w*vx;
    x[2] = w*vy;
    x[3] = w*vz;
    x[4] = w*vx*vx;
    x[5] = w*vy*vy;
    x[6] = w*vz*vz;
    x[7] = w*vx*vy;
    x[8] = w*vy*vz;
    x[9] = w*vz*vx;
    
    for(int j=0; j<10; j++){
      for(int k=j; k<10; k++)
        A[j+1][k+1] += x[j]*x[k];
    }
    
    //near points
    for(int j=0; j<adN; j++){
      float* q = ad_point[j];
      float d = (p[0]-q[0])*(p[0]-q[0]) + 
        (p[1]-q[1])*(p[1]-q[1]) +
          (p[2]-q[2])*(p[2]-q[2]);
      int insert = -1;
      for(int k=0; k<NEAR; k++){
        if(ad_d[j][k] > d)
          insert = k;
        else
          break;
      }
      if(insert < 0)
        continue;
      for(int k=0; k<insert; k++){
        ad_d[j][k] = ad_d[j][k+1];
        ad_i[j][k] = ad_i[j][k+1];
      }
      ad_d[j][insert] = d;
      ad_i[j][insert] = in;
    }
  }
  
  for(int i=1; i<11; i++){
    for(int j=i; j<11; j++)
      A[i][j] /= (float)totalW;
  }
  
  //Extra points
  int count = 0;
  for(int i=0; i<adN; i++){
    float* p = ad_point[i];
    
    //distance (inner product with normal)
    double v = 0;
    for(int j=0; j<NEAR; j++){
      int in = ad_i[i][j];
      float *q = point[in];
      float *n = normal[in];
      ad_d[i][j]= n[0]*(p[0]-q[0]) + n[1]*(p[1]-q[1]) + n[2]*(p[2]-q[2]);
      v += ad_d[i][j];
    }
    v /= NEAR;
    
    //sign check
    bool flag = true;
    for(int j=1; j<NEAR; j++){
      if(ad_d[i][0]*ad_d[i][j] <= 0){
        flag = false;
        break;
      }
    }
    if(!flag)
      continue;
	count++;
    
    float x[10];
    
    float w = 1.0f/adN;
    
    float vx = p[0] - o[0];
    float vy = p[1] - o[1];
    float vz = p[2] - o[2];
    
    x[0] = w;
    x[1] = w*vx;
    x[2] = w*vy;
    x[3] = w*vz;
    x[4] = w*vx*vx;
    x[5] = w*vy*vy;
    x[6] = w*vz*vz;
    x[7] = w*vx*vy;
    x[8] = w*vy*vz;
    x[9] = w*vz*vx;
    
    for(int j=0; j<10; j++){
      for(int k=j; k<10; k++)
        A[j+1][k+1] += x[j]*x[k];
      b[j+1] += (float)(x[j]*v*w);
    }
  }
  delete[] ad_d;
  delete[] ad_point;
  delete[] ad_i;
  
  for(int i=2; i<11; i++)
    for(int j=1; j<i; j++)
      A[i][j] = A[j][i];
  
  float w[11];
  float **v = new float*[11];
  for(int i=1; i<11; i++)
    v[i] = new float[11];
  SVD::svdcmp(A, 10, 10, w, v);
  
  float wmax=0.0f;
  for (int k=1;k<11;k++)
    if (fabs(w[k]) > wmax) wmax=(float)fabs(w[k]);
  
  if(wmax < 0.000000000001f || count == 0){
    _cxx = _cyy = _czz = _cxy = _cyz = _czx = _cx = _cy = _cz = 0;
	_c0 = 10000000;
    for(int i=1; i<11; i++){
      delete[] A[i];
      delete[] v[i];
    }
    delete[] A;
    delete[] v;
    return;
  }
  
  float wmin=wmax*0.0000001f;
  for (int k=1;k<11;k++){
    if (fabs(w[k]) < wmin) 
      w[k]=0.0;
  }
  
  float x[11];
  SVD::svbksb(A, w, v, 10, 10, b, x);
  
  _cxx = x[5];
  _cyy = x[6];
  _czz = x[7];
  
  _cxy = x[8];
  _cyz = x[9];
  _czx = x[10];
  
  _cx = x[2] - _cxy*o[1] - _czx*o[2] - 2.0f*_cxx*o[0];
  _cy = x[3] - _cyz*o[2] - _cxy*o[0] - 2.0f*_cyy*o[1];
  _cz = x[4] - _czx*o[0] - _cyz*o[1] - 2.0f*_czz*o[2];
  
  _c0 = x[1] - x[2]*o[0] - x[3]*o[1] - x[4]*o[2] 
        + _cxy*o[0]*o[1] + _cyz*o[1]*o[2] + _czx*o[2]*o[0]
        + _cxx*o[0]*o[0] + _cyy*o[1]*o[1] + _czz*o[2]*o[2];
  
  
  for(int i=1; i<11; i++){
    delete[] A[i];
    delete[] v[i];
  }
  delete[] A;
  delete[] v;
}

float Quadric::computeMaxError(PointSet* ps, float R, ImplicitOctCell *cell, 
                               int* index_list, int listN){
  float error = 0;
  
  float c[3];
  cell->cellCenter(c);
  float (*point)[3] = ps->_point;
  bool *bound = ps->_bound;
  
  float R2 = R*R;
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
    
    float f = (float)fabs(value(p[0], p[1], p[2]));
    float g[3];
    gradient(g, p[0], p[1], p[2]);
    float e = f/(float)sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);
    if(e > error)
      error = e;
  }
  return error;
}