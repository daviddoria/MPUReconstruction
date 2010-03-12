#include "ImplicitOctTree.h"
#include <limits>

ImplicitOctCell::ImplicitOctCell(ImplicitOctTree* tree, unsigned char level, 
                                 float cx, float cy, float cz, float size){
  _tree = tree;
  
  _sub_cell = NULL;
  _level = level;
  _leaf = true;
  
  _cx = cx;
  _cy = cy;
  _cz = cz;
  _size = size;
  
  _laf = NULL;
  _checked = false;
  _error = std::numeric_limits<float>::infinity();
  
  _tree->_cell_count++;
}

ImplicitOctCell::~ImplicitOctCell(){
  if(!_leaf){
    for(int i=0; i<8; i++)
      delete _sub_cell[i];
    delete[] _sub_cell;
  }
  if(_laf != NULL)
    delete _laf;
  
  _tree->_cell_count--;
}

float ImplicitOctCell::supportSize(){
  return _tree->_support*TWOSQRT3*_size;
}

void ImplicitOctCell::evaluateFunction(float x, float y, float z, 
                                       int until, float e, double &wf, double &w){
  double d = sqrt((x-_cx)*(x-_cx) + (y-_cy)*(y-_cy) + (z-_cz)*(z-_cz));
  float R = supportSize();
  if(d > R)
    return;
  
  if(!_checked){
    computeLA();
    /*
    if(_error > e && _level < until){
      delete _laf;
      _laf = NULL;
    }*/
  }
  
  if(_error > e && _level < until){
    if(_laf != NULL){
       delete _laf;
      _laf = NULL;
    }
    
    if(_leaf)
      split();
    
    for(int i=0; i<8; i++)
      _sub_cell[i]->evaluateFunction(x, y, z, until, e, wf, w);
  }	
  else{
    if(_laf == NULL)
      computeLA();
    
    //evaluate LA
    double w1 = weight(d, R);
    wf += w1* _laf->value(x,y,z);
    w += w1;
  }
}

void ImplicitOctCell::buildFunction(int until, float e){
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
      _sub_cell[i]->buildFunction(until, e);
  }
}

void ImplicitOctCell::evaluateGradient(float x, float y, float z, int until, float e, 
                                       double &wf, double &w, double gwf[], double gw[]){
  float vx = x - _cx;
  float vy = y - _cy;
  float vz = z - _cz;
  double d = sqrt(vx*vx + vy*vy + vz*vz);
  float R = supportSize();
  if(d > R)
    return;
  
  if(!_checked){
    computeLA();
    /*
    if(_error > e && _level < until){
      delete _laf;
      _laf = NULL;
    }*/
  }
  
  if(_error > e && _level < until){
    if(_laf != NULL){
       delete _laf;
      _laf = NULL;
    }
    
    if(_leaf)
      split();
    
    for(int i=0; i<8; i++)
      _sub_cell[i]->evaluateGradient(x, y, z, until, e, wf, w, gwf, gw);
  }
  else{
    if(_laf == NULL)
      computeLA();
    
    double w1 = weight(d, R);
    float gf[3];
    float f = _laf->valueAndGradient(gf, x, y, z);
    double gw1[3];
    weightG(gw1, vx, vy, vz, d, R);
    wf += w1*f;
    
    w += w1;
    
    gwf[0] += (gw1[0]*f + w1*gf[0]);
    gwf[1] += (gw1[1]*f + w1*gf[1]);
    gwf[2] += (gw1[2]*f + w1*gf[2]);
    
    gw[0] += gw1[0];
    gw[1] += gw1[1];
    gw[2] += gw1[2];
  }
}

//Is not used
void ImplicitOctCell::smoothNormal(int level, float (*temp_normal)[3]){
  /*
  if(_end == _start)
    return;
  if(_level < level){
    if(_leaf)
      split();
    if(_leaf)
      return;
    for(int i=0; i<8; i++)
      _sub_cell[i]->smoothNormal(level, temp_normal);
  }
  else{
    float R = supportSize();
    float c[3];
    cellCenter(c);
    
    int* index_list;
    int listN = _tree->collectPointsInSphere(index_list, c, R);
    
    if(listN == 0)
      return;
    
    float (*point)[3] = _tree->_ps->_point;
    float (*normal)[3] = _tree->_ps->_normal;
    
    float o[] = {0,0,0};
    float n[] = {0,0,0};
    float totalW = 0;
    for(int i=0; i<listN; i++){
      int j = index_list[i];
      float* q = point[j];
      float* m = normal[j];
      
      float w = (float)theWeight(q[0], q[1], q[2]);
      if(w == 0)
        continue;
      totalW += w;
      
      o[0] += w*q[0];
      o[1] += w*q[1];
      o[2] += w*q[2];
      
      n[0] += w*m[0];
      n[1] += w*m[1];
      n[2] += w*m[2];
    }
    if(totalW == 0)
      return;
    
    double l = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    if((float)l == 0)
      return;
    n[0] = (float)(n[0]/l);
    n[1] = (float)(n[1]/l);
    n[2] = (float)(n[2]/l);
    
    
    for(i=0; i<listN; i++){
      int j = index_list[i];
      float* q = point[j];
      
      float w = (float)theWeight(q[0], q[1], q[2]);
      temp_normal[j][0] += w*n[0];
      temp_normal[j][1] += w*n[1];
      temp_normal[j][2] += w*n[2];
    }
  }*/
}

//Is not used
void ImplicitOctCell::computeConfidence(int target_level, float T){
  /*
  if(_end == _start)
    return;
  if(_level < target_level){
    if(_leaf)
      split();
    if(_leaf)
      return;
    for(int i=0; i<8; i++)
      _sub_cell[i]->computeConfidence(target_level, T);
  }
  else{
    float R = supportSize();
    float c[3];
    cellCenter(c);
    
    int* index_list;
    int listN = _tree->collectPointsInSphere(index_list, c, R);
    
    if(listN == 0)
      return;
    
    float (*point)[3] = _tree->_ps->_point;
    float (*normal)[3] = _tree->_ps->_normal;
    
    float o[] = {0,0,0};
    float n[] = {0,0,0};
    float totalW = 0;
    for(int i=0; i<listN; i++){
      int j = index_list[i];
      float* q = point[j];
      float* m = normal[j];
      
      float w = (float)theWeight(q[0], q[1], q[2]);
      if(w == 0)
        continue;
      totalW += w;
      
      o[0] += w*q[0];
      o[1] += w*q[1];
      o[2] += w*q[2];
      
      n[0] += w*m[0];
      n[1] += w*m[1];
      n[2] += w*m[2];
    }
    if(totalW == 0)
      return;
    
    o[0] /= totalW;
    o[1] /= totalW;
    o[2] /= totalW;
    
    double l = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    if((float)l == 0)
      return;
    n[0] = (float)(n[0]/l);
    n[1] = (float)(n[1]/l);
    n[2] = (float)(n[2]/l);
    
    float dot = n[0]*(o[0]-c[0]) + n[1]*(o[1]-c[1]) + n[2]*(o[2]-c[2]);
    c[0] += dot*n[0];
    c[1] += dot*n[1];
    c[2] += dot*n[2];
    
    float v = (float)(sqrt((c[0]-o[0])*(c[0]-o[0]) + 
                           (c[1]-o[1])*(c[1]-o[1]) + 
                           (c[2]-o[2])*(c[2]-o[2]))/sqrt(R*R - dot*dot));		
    
    bool flag = false;
    if(1.0-v < T)
      flag = true;
    
    
    for(i=0; i<listN; i++){
      int j = index_list[i];
      _tree->_ps->_bound[j] |= flag;
    }
  }*/
}

void ImplicitOctCell::computeLA(){
  _checked = true;
  
  float R = supportSize();
  
  float c[3];
  cellCenter(c);
  float R1 = R;
  
  PointSet* _ps = _tree->_ps;
  float (*point)[3] = _ps->_point;
  float (*normal)[3] = _ps->_normal;
//  bool *bound = _ps->_bound;
  
  //if(_level < 3){
    //return;
  //}
  
  int listN, *index_list;
  while(true){
    listN = _ps->_tree->collectPointsInSphere(index_list, c, R1);
    if(listN > _tree->_Nmin)
      break;
    else
      R1 += R*_tree->_lambda;
  }
  
  //if(R1 != R && listN > _ps->_pointN*0.1)
    //return;
  
  int i;
  //averaged tangent plane
  float n[] = {0,0,0};
  float o[] = {0,0,0};
  float totalW = 0;
  for(i=0; i<listN; i++){
    int j = index_list[i];
    float* q = point[j];
    float* m = normal[j];
    
    double d = sqrt((q[0]-c[0])*(q[0]-c[0]) + 
                    (q[1]-c[1])*(q[1]-c[1]) + 
                    (q[2]-c[2])*(q[2]-c[2]));
    
    float w = (float)weight(d, R1);
    totalW += w;
    
    o[0] += w*q[0];
    o[1] += w*q[1];
    o[2] += w*q[2];
    
    n[0] += w*m[0];
    n[1] += w*m[1];
    n[2] += w*m[2];
  }
  o[0] /= totalW;
  o[1] /= totalW;
  o[2] /= totalW;
  
  double l = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  if((float)l == 0)
    return;
  n[0] = (float)(n[0]/l);
  n[1] = (float)(n[1]/l);
  n[2] = (float)(n[2]/l);
  
  //deciding type of polynomial approximation
  if(!_tree->_sharp || listN > 2*_tree->_Nmin){
    //general quadric or normal oriented quadric ?
    float T = 0.0f;
    bool flag = false;
    for(int i=0; i<listN; i++){
      float* m = normal[index_list[i]];
      if(n[0]*m[0] + n[1]*m[1] + n[2]*m[2] < T){
        flag = true;
        break;
      }
    }
    if(flag)
      _laf = new Quadric(_ps, R, R1, this, index_list, listN, o);
    else
      _laf = new QuadricO(_ps, R, R1, this, index_list, listN, o, n);
  }
  else{
    //normal oriented quadric or sharp feature function ?
    
    //From Kobbelt SIG'01
    int mini, minj;
    float min = 1;
    for(i=0; i<listN; i++){
      float* ni = normal[index_list[i]];
      for(int j=i+1; j<listN; j++){
        float* nj = normal[index_list[j]];
        float v = ni[0]*nj[0] + ni[1]*nj[1] + ni[2]*nj[2];
        if(v < min){
          min = v;
          mini = index_list[i];
          minj = index_list[j];
        }
      }
    }
    
    if(min < _tree->_edgeT){
      float* ni = normal[mini];
      float* nj = normal[minj];
      float nx = ni[1]*nj[2] - ni[2]*nj[1];
      float ny = ni[2]*nj[0] - ni[0]*nj[2];
      float nz = ni[0]*nj[1] - ni[1]*nj[0];
      float max = 0;
      int maxi;
      for(i=0; i<listN; i++){
        float* m = normal[index_list[i]];
        float v = (float)fabs(nx*m[0] + ny*m[1] + nz*m[2]);
        if(max < v){
          max = v;
          maxi = index_list[i];
        }
      }
      if(max < _tree->_cornerT)
        _laf = new QuadricOEdge(_ps, R, R1, this, index_list, listN, mini, minj);
      else{
        _laf = new Quadric(_ps, R, R1, this, index_list, listN, o);
        bool OK;
        _laf = new QuadricOCorner(_ps, R, R1, this, index_list, listN, mini, minj, maxi, OK, _tree->_edgeT);
        
        //If failed, quadric is used.
        if(!OK){
          delete _laf;
          //_laf = new QuadricO(_ps, R, R1, this, index_list, listN, o, n);
          _laf = new Quadric(_ps, R, R1, this, index_list, listN, o);
        }
      }
    }
    else
      _laf = new QuadricO(_ps, R, R1, this, index_list, listN, o, n);
  }
}

//Split cell 
void ImplicitOctCell::split(){
  _sub_cell = new ImplicitOctCell*[8];
  
  _leaf = false;
  
  float s = 0.5f*_size;
  
  _sub_cell[0] = new ImplicitOctCell(_tree, _level+1, _cx-s, _cy-s, _cz-s, s);
  _sub_cell[1] = new ImplicitOctCell(_tree, _level+1, _cx+s, _cy-s, _cz-s, s);
  _sub_cell[2] = new ImplicitOctCell(_tree, _level+1, _cx-s, _cy+s, _cz-s, s);
  _sub_cell[3] = new ImplicitOctCell(_tree, _level+1, _cx+s, _cy+s, _cz-s, s);
  _sub_cell[4] = new ImplicitOctCell(_tree, _level+1, _cx-s, _cy-s, _cz+s, s);
  _sub_cell[5] = new ImplicitOctCell(_tree, _level+1, _cx+s, _cy-s, _cz+s, s);                    
  _sub_cell[6] = new ImplicitOctCell(_tree, _level+1, _cx-s, _cy+s, _cz+s, s);
  _sub_cell[7] = new ImplicitOctCell(_tree, _level+1, _cx+s, _cy+s, _cz+s, s);
}
