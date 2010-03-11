// Polygonizer.cpp: Polygonizer �N���X�̃C���v�������e�[�V����
//
//////////////////////////////////////////////////////////////////////

//#include "stdafx.h"
#include "Polygonizer.h"
#include <vnl/algo/vnl_svd.h>
#include "PolygonalMesh.h"
#include "Bloomenthal.h"

//////////////////////////////////////////////////////////////////////
// �\�z/���
//////////////////////////////////////////////////////////////////////

Polygonizer::Polygonizer()
{

}

Polygonizer::~Polygonizer()
{

}

void Polygonizer::setDim(int dimX, int dimY, int dimZ)
{
	this->dimX = dimX;
	this->dimY = dimY;
	this->dimZ = dimZ;
}

void Polygonizer::setOrigin(float x, float y, float z)
{
	originX = x;
	originY = y;
	originZ = z;
}

void Polygonizer::setSpace(float x, float y, float z)
{
	spaceX = x;
	spaceY = y;
	spaceZ = z;
}

PolygonalMesh* Polygonizer::computeSurfaceNetSIG02(float epsilon, float tau)
{
	float p[3];
	bool ***isIn = new bool**[dimZ+1];
	bool ***isValid = new bool**[dimZ+1];
	for(int i=0; i<dimZ+1; i++){
		isIn[i] = new bool*[dimY+1];
		isValid[i] = new bool*[dimY+1];
		p[2] = originZ + i*spaceZ;
		for(int j=0; j<dimY+1; j++){
			isIn[i][j] = new bool[dimX+1];
			isValid[i][j] = new bool[dimX+1];
			p[1] = originY + j*spaceY;
			for(int k=0; k<dimX+1; k++){
				p[0] = originX + k*spaceX;
				isIn[i][j][k] 
					= (func->value(p[0], p[1], p[2], isValid[i][j][k]) > 0);
			}
		}
	}

	int ***index = new int**[dimZ+1];
	for(int i=0; i<dimZ+1; i++){
		index[i] = new int*[dimY+1];
		for(int j=0; j<dimY+1; j++){
			index[i][j] = new int[dimX+1];
			for(int k=0; k<dimX+1; k++)
				index[i][j][k] = -1;
		}
	}

	int current = 0;
	int face_N = 0;
	for(int i=0; i<dimZ; i++)
		for(int j=0; j<dimY; j++)
			for(int k=0; k<dimX-1; k++){
				if(!isValid[i][j][k] || !isValid[i][j][k+1])
					continue;
				if((!isIn[i][j][k] && isIn[i][j][k+1]) || 
				    (isIn[i][j][k] && !isIn[i][j][k+1])){
					face_N++;
					if(index[i][j][k+1] < 0){
						index[i][j][k+1] = current;
						current++;
					}
					if(index[i][j+1][k+1] < 0){
						index[i][j+1][k+1] = current;
						current++;
					}
					if(index[i+1][j][k+1] < 0){
						index[i+1][j][k+1] = current;
						current++;
					}
					if(index[i+1][j+1][k+1] < 0){
						index[i+1][j+1][k+1] = current;
						current++;
					}
				}
			}

	for(int i=0; i<dimX; i++)
		for(int j=0; j<dimZ; j++)
			for(int k=0; k<dimY-1; k++){
				if(!isValid[j][k][i] || !isValid[j][k+1][i])
					continue;
				if((!isIn[j][k][i] && isIn[j][k+1][i]) || 
				    (isIn[j][k][i] && !isIn[j][k+1][i])){
					face_N++;
					if(index[j][k+1][i] < 0){
						index[j][k+1][i] = current;
						current++;
					}
					if(index[j+1][k+1][i] < 0){
						index[j+1][k+1][i] = current;
						current++;
					}
					if(index[j][k+1][i+1] < 0){
						index[j][k+1][i+1] = current;
						current++;
					}
					if(index[j+1][k+1][i+1] < 0){
						index[j+1][k+1][i+1] = current;
						current++;
					}
				}
			}

	for(int i=0; i<dimY; i++)
		for(int j=0; j<dimX; j++)
			for(int k=0; k<dimZ-1; k++){
				if(!isValid[k][i][j] || !isValid[k+1][i][j])
					continue;
				if((!isIn[k][i][j] && isIn[k+1][i][j]) || 
				    (isIn[k][i][j] && !isIn[k+1][i][j])){
					face_N++;
					if(index[k+1][i][j] < 0){
						index[k+1][i][j] = current;
						current++;
					}
					if(index[k+1][i+1][j] < 0){
						index[k+1][i+1][j] = current;
						current++;
					}
					if(index[k+1][i][j+1] < 0){
						index[k+1][i][j+1] = current;
						current++;
					}
					if(index[k+1][i+1][j+1] < 0){
						index[k+1][i+1][j+1] = current;
						current++;
					}
				}
			}
	
	PolygonalMesh* mesh = new PolygonalMesh;
	int vertex_N = current;
	mesh->setVertexCount(vertex_N);
	float (*vertex)[3] = mesh->vertex;
	int* degree = mesh->degree_f = new int[vertex_N];
	current = 0;
	for(int i=0; i<vertex_N; i++){
		vertex[i][0] = vertex[i][1] = vertex[i][2] = 0;
		degree[i] = 0;
	}
	/*
	for(i=0; i<dimZ+1; i++)
		for(int j=0; j<dimY+1; j++)
			for(int k=0; k<dimX+1; k++)
				if(index[i][j][k] >= 0){
					index[i][j][k] = current;
					vertex[current][0] = spaceX*((float)k-0.5f) + originX;
					vertex[current][1] = spaceY*((float)j-0.5f) + originY;
					vertex[current][2] = spaceZ*((float)i-0.5f) + originZ;
					current++;
				}
				*/


	mesh->setFaceCount(face_N);
	double (*Q)[10] = new double[vertex_N][10];
	for(int i=0; i<vertex_N; i++)
		MAT_INIT(Q[i]);
	for(int i=0; i<face_N; i++)
		mesh->setPolygonCount(i, 4);
	face_N = 0;
	int **face = mesh->face;
	bool flag = false;
	for(int i=0; i<dimZ; i++)
		for(int j=0; j<dimY; j++)
			for(int k=0; k<dimX-1; k++){
				if(!isValid[i][j][k] || !isValid[i][j][k+1])
					continue;
				if(isIn[i][j][k] && !isIn[i][j][k+1]){
					face[face_N][0] = index[i][j][k+1];
					face[face_N][1] = index[i][j+1][k+1];
					face[face_N][2] = index[i+1][j+1][k+1];
					face[face_N][3] = index[i+1][j][k+1];
					face_N++;
					flag = true;
				}
				else if(!isIn[i][j][k] && isIn[i][j][k+1]){
					face[face_N][0] = index[i][j][k+1];
					face[face_N][1] = index[i+1][j][k+1];
					face[face_N][2] = index[i+1][j+1][k+1];
					face[face_N][3] = index[i][j+1][k+1];
					face_N++;
					flag = true;
				}
				if(!flag)
					continue;
				flag = false;
				float p[3], s[3], e[3];
				s[0] = originX + k*spaceX;
				e[0] = originX + (k+1)*spaceX;
				s[1] = e[1] = originY + j*spaceY;
				s[2] = e[2] = originZ + i*spaceZ;
				bisection(p, s, e, epsilon);
				
				float g[3];
				func->gradient(g, p[0], p[1], p[2]);
				double len = PolygonalMesh::LENGTH(g);
				//if((float)len == 0)
					//continue;
				double nor[3];
				nor[0] = g[0]/len;
				nor[1] = g[1]/len;
				nor[2] = g[2]/len;

				double d = -PolygonalMesh::DOT(nor, p);
				double Q_tmp[10];
				MATRIX(Q_tmp, nor, d);

				int i0 = index[i][j][k+1];
				MAT_SUM(Q[i0], Q_tmp);
				vertex[i0][0] += p[0];
				vertex[i0][1] += p[1];
				vertex[i0][2] += p[2];
				degree[i0]++;

				i0 = index[i][j+1][k+1];
				MAT_SUM(Q[i0], Q_tmp);
				vertex[i0][0] += p[0];
				vertex[i0][1] += p[1];
				vertex[i0][2] += p[2];
				degree[i0]++;

				i0 = index[i+1][j+1][k+1];
				MAT_SUM(Q[i0], Q_tmp);
				vertex[i0][0] += p[0];
				vertex[i0][1] += p[1];
				vertex[i0][2] += p[2];
				degree[i0]++;

				i0 = index[i+1][j][k+1];
				MAT_SUM(Q[i0], Q_tmp);
				vertex[i0][0] += p[0];
				vertex[i0][1] += p[1];
				vertex[i0][2] += p[2];
				degree[i0]++;
			}

	for(int i=0; i<dimX; i++)
		for(int j=0; j<dimZ; j++)
			for(int k=0; k<dimY-1; k++){
				if(!isValid[j][k][i] || !isValid[j][k+1][i])
					continue;
				if(isIn[j][k][i] && !isIn[j][k+1][i]){
					face[face_N][0] = index[j][k+1][i];
					face[face_N][1] = index[j+1][k+1][i];
					face[face_N][2] = index[j+1][k+1][i+1];
					face[face_N][3] = index[j][k+1][i+1];
					face_N++;
					flag = true;
				}
				else if(!isIn[j][k][i] && isIn[j][k+1][i]){
					face[face_N][0] = index[j][k+1][i];
					face[face_N][1] = index[j][k+1][i+1];
					face[face_N][2] = index[j+1][k+1][i+1];
					face[face_N][3] = index[j+1][k+1][i];
					face_N++;
					flag = true;
				}
				if(!flag)
					continue;
				flag = false;
				float p[3], s[3], e[3];
				s[1] = originY + k*spaceY;
				e[1] = originY + (k+1)*spaceY;
				s[2] = e[2] = originZ + j*spaceZ;
				s[0] = e[0] = originX + i*spaceX;
				bisection(p, s, e, epsilon);
				
				float g[3];
				func->gradient(g, p[0], p[1], p[2]);
				double len = PolygonalMesh::LENGTH(g);
				//if((float)len == 0)
					//continue;
				double nor[3];
				nor[0] = g[0]/len;
				nor[1] = g[1]/len;
				nor[2] = g[2]/len;

				double d = -PolygonalMesh::DOT(nor, p);
				double Q_tmp[10];
				MATRIX(Q_tmp, nor, d);
				
				int i0 = index[j][k+1][i];
				MAT_SUM(Q[i0], Q_tmp);
				vertex[i0][0] += p[0];
				vertex[i0][1] += p[1];
				vertex[i0][2] += p[2];
				degree[i0]++;

				i0 = index[j+1][k+1][i];
				MAT_SUM(Q[i0], Q_tmp);
				vertex[i0][0] += p[0];
				vertex[i0][1] += p[1];
				vertex[i0][2] += p[2];
				degree[i0]++;

				i0 = index[j+1][k+1][i+1];
				MAT_SUM(Q[i0], Q_tmp);
				vertex[i0][0] += p[0];
				vertex[i0][1] += p[1];
				vertex[i0][2] += p[2];
				degree[i0]++;

				i0 = index[j][k+1][i+1];
				MAT_SUM(Q[i0], Q_tmp);
				vertex[i0][0] += p[0];
				vertex[i0][1] += p[1];
				vertex[i0][2] += p[2];
				degree[i0]++;
			}

	for(int i=0; i<dimY; i++)
		for(int j=0; j<dimX; j++)
			for(int k=0; k<dimZ-1; k++){
				if(!isValid[k][i][j] || !isValid[k+1][i][j])
					continue;
				if(isIn[k][i][j] && !isIn[k+1][i][j]){
					face[face_N][0] = index[k+1][i][j];
					face[face_N][1] = index[k+1][i][j+1];
					face[face_N][2] = index[k+1][i+1][j+1];
					face[face_N][3] = index[k+1][i+1][j];
					face_N++;
					flag = true;
				}
				else if(!isIn[k][i][j] && isIn[k+1][i][j]){
					face[face_N][0] = index[k+1][i][j];
					face[face_N][1] = index[k+1][i+1][j];
					face[face_N][2] = index[k+1][i+1][j+1];
					face[face_N][3] = index[k+1][i][j+1];
					face_N++;
					flag = true;
				}
				if(!flag)
					continue;
				flag = false;
				float p[3], s[3], e[3];
				s[2] = originZ + k*spaceZ;
				e[2] = originZ + (k+1)*spaceZ;
				s[0] = e[0] = originX + j*spaceX;
				s[1] = e[1] = originY + i*spaceY;
				bisection(p, s, e, epsilon);
				
				float g[3];
				func->gradient(g, p[0], p[1], p[2]);
				double len = PolygonalMesh::LENGTH(g);
				//if((float)len == 0)
					//continue;
				double nor[3];
				nor[0] = g[0]/len;
				nor[1] = g[1]/len;
				nor[2] = g[2]/len;

				double d = -PolygonalMesh::DOT(nor, p);
				double Q_tmp[10];
				MATRIX(Q_tmp, nor, d);

				int i0 = index[k+1][i][j];
				MAT_SUM(Q[i0], Q_tmp);
				vertex[i0][0] += p[0];
				vertex[i0][1] += p[1];
				vertex[i0][2] += p[2];
				degree[i0]++;

				i0 = index[k+1][i][j+1];
				MAT_SUM(Q[i0], Q_tmp);
				vertex[i0][0] += p[0];
				vertex[i0][1] += p[1];
				vertex[i0][2] += p[2];
				degree[i0]++;

				i0 = index[k+1][i+1][j+1];
				MAT_SUM(Q[i0], Q_tmp);
				vertex[i0][0] += p[0];
				vertex[i0][1] += p[1];
				vertex[i0][2] += p[2];
				degree[i0]++;

				i0 = index[k+1][i+1][j];
				MAT_SUM(Q[i0], Q_tmp);
				vertex[i0][0] += p[0];
				vertex[i0][1] += p[1];
				vertex[i0][2] += p[2];
				degree[i0]++;
			}
	
	//FOR SVD
			
	vnl_matrix< float > A( 3, 3, 0. ); 
	vnl_vector< float > b(3, 0.), x(3, 0.);

	for(int i=0; i<vertex_N; i++){
		if(degree[i] == 0)
			continue;
		vertex[i][0] /= degree[i];
		vertex[i][1] /= degree[i];
		vertex[i][2] /= degree[i];
continue;
		A[0][0] = (float)Q[i][0];
		A[1][0] = A[0][1] = (float)Q[i][1];
		A[2][0] = A[0][2] = (float)Q[i][2];
		A[1][1] = (float)Q[i][3];
		A[1][2] = A[2][1] = (float)Q[i][4];
		A[2][2] = (float)Q[i][5];

		float Av[3];
		MAT_BY_VEC(Av, Q[i], vertex[i]);
		b[0] = -(float)Q[i][6] - Av[0];
		b[1] = -(float)Q[i][7] - Av[1];
		b[2] = -(float)Q[i][8] - Av[2];

    vnl_svd<float> svd( A );
    svd.zero_out_absolute( 0.0000001 );

    x = svd.solve( b );

    if(fabs(x[1]) > spaceX || fabs(x[2]) > spaceY || fabs(x[3]) > spaceZ)
			continue;

		mesh->vertex[i][0] += x[1];
		mesh->vertex[i][1] += x[2];
		mesh->vertex[i][2] += x[3];
	}

	for(int i=0; i<dimZ+1; i++){
		for(int j=0; j<dimY+1; j++){
			delete[] isIn[i][j];
			delete[] index[i][j];
		}
		delete[] isIn[i];
		delete[] index[i];
	}
	delete[] isIn;
	delete[] index;
	delete[] Q;
	delete[] isValid;

	return mesh;
}

void Polygonizer::bisection(float p[], float start[], float end[], float e)
{
	float p1[3], p2[3], p3[3];
	float f1, f2, f3;
	p1[0] = start[0];
	p1[1] = start[1];
	p1[2] = start[2];
	p2[0] = end[0];
	p2[1] = end[1];
	p2[2] = end[2];
	f1 = func->value(p1[0], p1[1], p1[2]);
	f2 = func->value(p2[0], p2[1], p2[2]);
	for(int j=0; j<20; j++){
		if(PolygonalMesh::DIST(p1, p2) < 0.001){
			p[0] = p3[0];
			p[1] = p3[1];
			p[2] = p3[2];
		}

		p3[0] = 0.5f*(p1[0] + p2[0]);
		p3[1] = 0.5f*(p1[1] + p2[1]);
		p3[2] = 0.5f*(p1[2] + p2[2]);
		f3 = func->value(p3[0], p3[1], p3[2]);
		if(fabs(f3) < 0.00001){
			p[0] = p3[0];
			p[1] = p3[1];
			p[2] = p3[2];
		}
		else if(f1*f3 >= 0){
			p1[0] = p3[0];
			p1[1] = p3[1];
			p1[2] = p3[2];
			f1 = f3;
		}
		else{
			p2[0] = p3[0];
			p2[1] = p3[1];
			p2[2] = p3[2];
			f2 = f3;
		}
	}
}

PolygonalMesh* Polygonizer::computeSurfaceNet()
{
	float p[3];
	bool ***isIn = new bool**[dimZ+1];
	bool ***isValid = new bool**[dimZ+1];
	for(int i=0; i<dimZ+1; i++){
		isIn[i] = new bool*[dimY+1];
		isValid[i] = new bool*[dimY+1];
		p[2] = originZ + i*spaceZ;
		for(int j=0; j<dimY+1; j++){
			isIn[i][j] = new bool[dimX+1];
			isValid[i][j] = new bool[dimX+1];
			p[1] = originY + j*spaceY;
			for(int k=0; k<dimX+1; k++){
				p[0] = originX + k*spaceX;
				isIn[i][j][k] 
					= (func->value(p[0], p[1], p[2], isValid[i][j][k]) > 0);
			}
		}
	}

	int ***index = new int**[dimZ+1];
	for(int i=0; i<dimZ+1; i++){
		index[i] = new int*[dimY+1];
		for(int j=0; j<dimY+1; j++){
			index[i][j] = new int[dimX+1];
			for(int k=0; k<dimX+1; k++)
				index[i][j][k] = -1;
		}
	}

	int current = 0;
	int face_N = 0;
	for(int i=0; i<dimZ; i++)
		for(int j=0; j<dimY; j++)
			for(int k=0; k<dimX-1; k++){
				if(!isValid[i][j][k] || !isValid[i][j][k+1])
					continue;
				if((!isIn[i][j][k] && isIn[i][j][k+1]) || 
				    (isIn[i][j][k] && !isIn[i][j][k+1])){
					face_N++;
					if(index[i][j][k+1] < 0){
						index[i][j][k+1] = current;
						current++;
					}
					if(index[i][j+1][k+1] < 0){
						index[i][j+1][k+1] = current;
						current++;
					}
					if(index[i+1][j][k+1] < 0){
						index[i+1][j][k+1] = current;
						current++;
					}
					if(index[i+1][j+1][k+1] < 0){
						index[i+1][j+1][k+1] = current;
						current++;
					}
				}
			}

	for(int i=0; i<dimX; i++)
		for(int j=0; j<dimZ; j++)
			for(int k=0; k<dimY-1; k++){
				if(!isValid[j][k][i] || !isValid[j][k+1][i])
					continue;
				if((!isIn[j][k][i] && isIn[j][k+1][i]) || 
				    (isIn[j][k][i] && !isIn[j][k+1][i])){
					face_N++;
					if(index[j][k+1][i] < 0){
						index[j][k+1][i] = current;
						current++;
					}
					if(index[j+1][k+1][i] < 0){
						index[j+1][k+1][i] = current;
						current++;
					}
					if(index[j][k+1][i+1] < 0){
						index[j][k+1][i+1] = current;
						current++;
					}
					if(index[j+1][k+1][i+1] < 0){
						index[j+1][k+1][i+1] = current;
						current++;
					}
				}
			}

	for(int i=0; i<dimY; i++)
		for(int j=0; j<dimX; j++)
			for(int k=0; k<dimZ-1; k++){
				if(!isValid[k][i][j] || !isValid[k+1][i][j])
					continue;
				if((!isIn[k][i][j] && isIn[k+1][i][j]) || 
				    (isIn[k][i][j] && !isIn[k+1][i][j])){
					face_N++;
					if(index[k+1][i][j] < 0){
						index[k+1][i][j] = current;
						current++;
					}
					if(index[k+1][i+1][j] < 0){
						index[k+1][i+1][j] = current;
						current++;
					}
					if(index[k+1][i][j+1] < 0){
						index[k+1][i][j+1] = current;
						current++;
					}
					if(index[k+1][i+1][j+1] < 0){
						index[k+1][i+1][j+1] = current;
						current++;
					}
				}
			}
	
	PolygonalMesh* mesh = new PolygonalMesh;
	int vertex_N = current;
	mesh->setVertexCount(vertex_N);
	float (*vertex)[3] = mesh->vertex;
	current = 0;
	for(int i=0; i<dimZ+1; i++)
		for(int j=0; j<dimY+1; j++)
			for(int k=0; k<dimX+1; k++)
				if(index[i][j][k] >= 0){
					index[i][j][k] = current;
					vertex[current][0] = spaceX*((float)k-0.5f) + originX;
					vertex[current][1] = spaceY*((float)j-0.5f) + originY;
					vertex[current][2] = spaceZ*((float)i-0.5f) + originZ;
					current++;
				}


	mesh->setFaceCount(face_N);
	for(int i=0; i<face_N; i++)
		mesh->setPolygonCount(i, 4);
	face_N = 0;
	int **face = mesh->face;
	for(int i=0; i<dimZ; i++)
		for(int j=0; j<dimY; j++)
			for(int k=0; k<dimX-1; k++){
				if(!isValid[i][j][k] || !isValid[i][j][k+1])
					continue;
				if(isIn[i][j][k] && !isIn[i][j][k+1]){
					face[face_N][0] = index[i][j][k+1];
					face[face_N][1] = index[i][j+1][k+1];
					face[face_N][2] = index[i+1][j+1][k+1];
					face[face_N][3] = index[i+1][j][k+1];
					face_N++;
				}
				else if(!isIn[i][j][k] && isIn[i][j][k+1]){
					face[face_N][0] = index[i][j][k+1];
					face[face_N][1] = index[i+1][j][k+1];
					face[face_N][2] = index[i+1][j+1][k+1];
					face[face_N][3] = index[i][j+1][k+1];
					face_N++;
				}
			}

	for(int i=0; i<dimX; i++)
		for(int j=0; j<dimZ; j++)
			for(int k=0; k<dimY-1; k++){
				if(!isValid[j][k][i] || !isValid[j][k+1][i])
					continue;
				if(isIn[j][k][i] && !isIn[j][k+1][i]){
					face[face_N][0] = index[j][k+1][i];
					face[face_N][1] = index[j+1][k+1][i];
					face[face_N][2] = index[j+1][k+1][i+1];
					face[face_N][3] = index[j][k+1][i+1];
					face_N++;
				}
				else if(!isIn[j][k][i] && isIn[j][k+1][i]){
					face[face_N][0] = index[j][k+1][i];
					face[face_N][1] = index[j][k+1][i+1];
					face[face_N][2] = index[j+1][k+1][i+1];
					face[face_N][3] = index[j+1][k+1][i];
					face_N++;
				}
			}

	for(int i=0; i<dimY; i++)
		for(int j=0; j<dimX; j++)
			for(int k=0; k<dimZ-1; k++){
				if(!isValid[k][i][j] || !isValid[k+1][i][j])
					continue;
				if(isIn[k][i][j] && !isIn[k+1][i][j]){
					face[face_N][0] = index[k+1][i][j];
					face[face_N][1] = index[k+1][i][j+1];
					face[face_N][2] = index[k+1][i+1][j+1];
					face[face_N][3] = index[k+1][i+1][j];
					face_N++;
				}
				else if(!isIn[k][i][j] && isIn[k+1][i][j]){
					face[face_N][0] = index[k+1][i][j];
					face[face_N][1] = index[k+1][i+1][j];
					face[face_N][2] = index[k+1][i+1][j+1];
					face[face_N][3] = index[k+1][i][j+1];
					face_N++;
				}
			}
	
	for(int i=0; i<dimZ+1; i++){
		for(int j=0; j<dimY+1; j++){
			delete[] isIn[i][j];
			delete[] index[i][j];
		}
		delete[] isIn[i];
		delete[] index[i];
	}
	delete[] isIn;
	delete[] index;

	return mesh;
}

void Polygonizer::cutQuad(PolygonalMesh* mesh)
{
	int face_N = mesh->face_N;
	int **face = mesh->face;
	//float (*vertex)[3] = mesh->vertex;
	int **new_face = new int*[2*face_N];
	delete[] mesh->poly_N;
	int *poly_N = mesh->poly_N = new int[2*face_N];
	int *degree = mesh->degree_f;
	for(int i=0; i<face_N; i++){
		poly_N[2*i] = poly_N[2*i+1] = 3;
		int *nf1 = new_face[2*i] = new int[3];
		int *nf2 = new_face[2*i+1] = new int[3];

		int *f = face[i];
		int i0 = f[0];
		int i1 = f[1];
		int i2 = f[2];
		int i3 = f[3];

		if(degree[i0] + degree[i2] < degree[i1] + degree[i3]){
			nf1[0] = i0;
			nf1[1] = i1;
			nf1[2] = i2;

			nf2[0] = i0;
			nf2[1] = i2;
			nf2[2] = i3;
		}
		else{
			nf1[0] = i0;
			nf1[1] = i1;
			nf1[2] = i3;

			nf2[0] = i1;
			nf2[1] = i2;
			nf2[2] = i3;
		}

		delete[] f;
	}
	delete[] face;
	mesh->face = new_face;
	mesh->face_N = 2*face_N;
	delete[] mesh->normal_f;
	mesh->normal_f = NULL;
	delete[] mesh->normal;
	mesh->normal = NULL;
}

PolygonalMesh* Polygonizer::computeSurfaceNerLinear(float epsilon, float tau)
{
	bool ***isIn = new bool**[dimZ+1];
	float ***value = new float**[dimZ+1];
	for(int i=0; i<dimZ+1; i++){
		isIn[i] = new bool*[dimY+1];
		value[i] = new float*[dimY+1];

		for(int j=0; j<dimY+1; j++){
			isIn[i][j] = new bool[dimX+1];
			value[i][j] = new float[dimX+1];

			for(int k=0; k<dimX+1; k++){
				isIn[i][j][k] = false;
				value[i][j][k] = 0;
			}
		}
	}

	float o[3];
	o[0] = originX;
	o[1] = originY;
	o[2] = originZ;
	float space[3];
	space[0] = spaceX;
	space[1] = spaceY;
	space[2] = spaceZ;
	int dim[3];
	dim[0] = dimX;
	dim[1] = dimY;
	dim[2] = dimZ;

	func->asignValueToVoxels(value, isIn, o, space, dim);

	int ***index = new int**[dimZ+1];
	for(int i=0; i<dimZ+1; i++){
		index[i] = new int*[dimY+1];
		for(int j=0; j<dimY+1; j++){
			index[i][j] = new int[dimX+1];
			for(int k=0; k<dimX+1; k++)
				index[i][j][k] = -1;
		}
	}

	int current = 0;
	int face_N = 0;
	for(int i=0; i<dimZ; i++)
		for(int j=0; j<dimY; j++)
			for(int k=0; k<dimX-1; k++){
				if(!isIn[i][j][k] || !isIn[i][j][k+1])
					continue;
				if((value[i][j][k] > 0 && value[i][j][k+1] <= 0) || 
				    (value[i][j][k] <= 0 && value[i][j][k+1] > 0)){
					face_N++;
					if(index[i][j][k+1] < 0){
						index[i][j][k+1] = current;
						current++;
					}
					if(index[i][j+1][k+1] < 0){
						index[i][j+1][k+1] = current;
						current++;
					}
					if(index[i+1][j][k+1] < 0){
						index[i+1][j][k+1] = current;
						current++;
					}
					if(index[i+1][j+1][k+1] < 0){
						index[i+1][j+1][k+1] = current;
						current++;
					}
				}
			}

	for(int i=0; i<dimX; i++)
		for(int j=0; j<dimZ; j++)
			for(int k=0; k<dimY-1; k++){
				if(!isIn[j][k][i] || !isIn[j][k+1][i])
					continue;
				if((value[j][k][i] > 0 && value[j][k+1][i] <= 0) || 
				    (value[j][k][i] <= 0&& value[j][k+1][i] > 0)){
					face_N++;
					if(index[j][k+1][i] < 0){
						index[j][k+1][i] = current;
						current++;
					}
					if(index[j+1][k+1][i] < 0){
						index[j+1][k+1][i] = current;
						current++;
					}
					if(index[j][k+1][i+1] < 0){
						index[j][k+1][i+1] = current;
						current++;
					}
					if(index[j+1][k+1][i+1] < 0){
						index[j+1][k+1][i+1] = current;
						current++;
					}
				}
			}

	for(int i=0; i<dimY; i++)
		for(int j=0; j<dimX; j++)
			for(int k=0; k<dimZ-1; k++){
				if(!isIn[k][i][j] || !isIn[k+1][i][j])
					continue;
				if((value[k][i][j] > 0  && value[k+1][i][j] <= 0) || 
				    (value[k][i][j] <= 0 && value[k+1][i][j] > 0)){
					face_N++;
					if(index[k+1][i][j] < 0){
						index[k+1][i][j] = current;
						current++;
					}
					if(index[k+1][i+1][j] < 0){
						index[k+1][i+1][j] = current;
						current++;
					}
					if(index[k+1][i][j+1] < 0){
						index[k+1][i][j+1] = current;
						current++;
					}
					if(index[k+1][i+1][j+1] < 0){
						index[k+1][i+1][j+1] = current;
						current++;
					}
				}
			}

	PolygonalMesh* mesh = new PolygonalMesh;
	int vertex_N = current;
	mesh->setVertexCount(vertex_N);
	float (*vertex)[3] = mesh->vertex;
	int* degree = mesh->degree_f = new int[vertex_N];
	current = 0;
	for(int i=0; i<vertex_N; i++){
		vertex[i][0] = vertex[i][1] = vertex[i][2] = 0;
		degree[i] = 0;
	}

	mesh->setFaceCount(face_N);
	double (*Q)[10] = new double[vertex_N][10];
	for(int i=0; i<vertex_N; i++)
		MAT_INIT(Q[i]);
	for(int i=0; i<face_N; i++)
		mesh->setPolygonCount(i, 4);
	face_N = 0;
	int **face = mesh->face;
	bool flag = false;
	for(int i=0; i<dimZ; i++)
		for(int j=0; j<dimY; j++)
			for(int k=0; k<dimX-1; k++){
				if(!isIn[i][j][k] || !isIn[i][j][k+1])
					continue;
				if(value[i][j][k] > 0 && value[i][j][k+1] <= 0){
					face[face_N][0] = index[i][j][k+1];
					face[face_N][1] = index[i][j+1][k+1];
					face[face_N][2] = index[i+1][j+1][k+1];
					face[face_N][3] = index[i+1][j][k+1];
					face_N++;
					flag = true;
				}
				else if(value[i][j][k] <= 0 && value[i][j][k+1] > 0){
					face[face_N][0] = index[i][j][k+1];
					face[face_N][1] = index[i+1][j][k+1];
					face[face_N][2] = index[i+1][j+1][k+1];
					face[face_N][3] = index[i][j+1][k+1];
					face_N++;
					flag = true;
				}
				if(!flag)
					continue;
				flag = false;
				float p[3], s[3], e[3];
				s[0] = originX + k*spaceX;
				e[0] = originX + (k+1)*spaceX;
				s[1] = e[1] = originY + j*spaceY;
				s[2] = e[2] = originZ + i*spaceZ;
				float f1 = value[i][j][k];
				float f2 = value[i][j][k+1];
				searchZero(p, s, e, f1, f2, epsilon);
				
				float g[3];
				//func->gradient(g, p[0], p[1], p[2]);
				double len = PolygonalMesh::LENGTH(g);
			
				double nor[3];
				nor[0] = g[0]/len;
				nor[1] = g[1]/len;
				nor[2] = g[2]/len;

				double d = -PolygonalMesh::DOT(nor, p);
				double Q_tmp[10];
				MATRIX(Q_tmp, nor, d);

				int i0 = index[i][j][k+1];
				MAT_SUM(Q[i0], Q_tmp);
				vertex[i0][0] += p[0];
				vertex[i0][1] += p[1];
				vertex[i0][2] += p[2];
				degree[i0]++;

				i0 = index[i][j+1][k+1];
				MAT_SUM(Q[i0], Q_tmp);
				vertex[i0][0] += p[0];
				vertex[i0][1] += p[1];
				vertex[i0][2] += p[2];
				degree[i0]++;

				i0 = index[i+1][j+1][k+1];
				MAT_SUM(Q[i0], Q_tmp);
				vertex[i0][0] += p[0];
				vertex[i0][1] += p[1];
				vertex[i0][2] += p[2];
				degree[i0]++;

				i0 = index[i+1][j][k+1];
				MAT_SUM(Q[i0], Q_tmp);
				vertex[i0][0] += p[0];
				vertex[i0][1] += p[1];
				vertex[i0][2] += p[2];
				degree[i0]++;
			}

	for(int i=0; i<dimX; i++)
		for(int j=0; j<dimZ; j++)
			for(int k=0; k<dimY-1; k++){
				if(!isIn[j][k][i] || !isIn[j][k+1][i])
					continue;
				if(value[j][k][i] > 0 && value[j][k+1][i] <= 0){
					face[face_N][0] = index[j][k+1][i];
					face[face_N][1] = index[j+1][k+1][i];
					face[face_N][2] = index[j+1][k+1][i+1];
					face[face_N][3] = index[j][k+1][i+1];
					face_N++;
					flag = true;
				}
				else if(value[j][k][i] <= 0 && value[j][k+1][i] > 0){
					face[face_N][0] = index[j][k+1][i];
					face[face_N][1] = index[j][k+1][i+1];
					face[face_N][2] = index[j+1][k+1][i+1];
					face[face_N][3] = index[j+1][k+1][i];
					face_N++;
					flag = true;
				}
				if(!flag)
					continue;
				flag = false;
				float p[3], s[3], e[3];
				s[1] = originY + k*spaceY;
				e[1] = originY + (k+1)*spaceY;
				s[2] = e[2] = originZ + j*spaceZ;
				s[0] = e[0] = originX + i*spaceX;
				float f1 = value[j][k][i];
				float f2 = value[j][k+1][i];
				searchZero(p, s, e, f1, f2, epsilon);
				
				float g[3];
				//func->gradient(g, p[0], p[1], p[2]);
				double len = PolygonalMesh::LENGTH(g);
			
				double nor[3];
				nor[0] = g[0]/len;
				nor[1] = g[1]/len;
				nor[2] = g[2]/len;

				double d = -PolygonalMesh::DOT(nor, p);
				double Q_tmp[10];
				MATRIX(Q_tmp, nor, d);
				
				int i0 = index[j][k+1][i];
				MAT_SUM(Q[i0], Q_tmp);
				vertex[i0][0] += p[0];
				vertex[i0][1] += p[1];
				vertex[i0][2] += p[2];
				degree[i0]++;

				i0 = index[j+1][k+1][i];
				MAT_SUM(Q[i0], Q_tmp);
				vertex[i0][0] += p[0];
				vertex[i0][1] += p[1];
				vertex[i0][2] += p[2];
				degree[i0]++;

				i0 = index[j+1][k+1][i+1];
				MAT_SUM(Q[i0], Q_tmp);
				vertex[i0][0] += p[0];
				vertex[i0][1] += p[1];
				vertex[i0][2] += p[2];
				degree[i0]++;

				i0 = index[j][k+1][i+1];
				MAT_SUM(Q[i0], Q_tmp);
				vertex[i0][0] += p[0];
				vertex[i0][1] += p[1];
				vertex[i0][2] += p[2];
				degree[i0]++;
			}

	for(int i=0; i<dimY; i++)
		for(int j=0; j<dimX; j++)
			for(int k=0; k<dimZ-1; k++){
				if(!isIn[k][i][j] || !isIn[k+1][i][j])
					continue;
				if(value[k][i][j] > 0 && value[k+1][i][j] <= 0){
					face[face_N][0] = index[k+1][i][j];
					face[face_N][1] = index[k+1][i][j+1];
					face[face_N][2] = index[k+1][i+1][j+1];
					face[face_N][3] = index[k+1][i+1][j];
					face_N++;
					flag = true;
				}
				else if(value[k][i][j] <= 0 && value[k+1][i][j] > 0){
					face[face_N][0] = index[k+1][i][j];
					face[face_N][1] = index[k+1][i+1][j];
					face[face_N][2] = index[k+1][i+1][j+1];
					face[face_N][3] = index[k+1][i][j+1];
					face_N++;
					flag = true;
				}
				if(!flag)
					continue;
				flag = false;
				float p[3], s[3], e[3];
				s[2] = originZ + k*spaceZ;
				e[2] = originZ + (k+1)*spaceZ;
				s[0] = e[0] = originX + j*spaceX;
				s[1] = e[1] = originY + i*spaceY;
				float f1 = value[k][i][j];
				float f2 = value[k+1][i][j];
				searchZero(p, s, e, f1, f2, epsilon);
				
				float g[3];
				//func->gradient(g, p[0], p[1], p[2]);
				double len = PolygonalMesh::LENGTH(g);
				
				double nor[3];
				nor[0] = g[0]/len;
				nor[1] = g[1]/len;
				nor[2] = g[2]/len;

				double d = -PolygonalMesh::DOT(nor, p);
				double Q_tmp[10];
				MATRIX(Q_tmp, nor, d);

				int i0 = index[k+1][i][j];
				MAT_SUM(Q[i0], Q_tmp);
				vertex[i0][0] += p[0];
				vertex[i0][1] += p[1];
				vertex[i0][2] += p[2];
				degree[i0]++;

				i0 = index[k+1][i][j+1];
				MAT_SUM(Q[i0], Q_tmp);
				vertex[i0][0] += p[0];
				vertex[i0][1] += p[1];
				vertex[i0][2] += p[2];
				degree[i0]++;

				i0 = index[k+1][i+1][j+1];
				MAT_SUM(Q[i0], Q_tmp);
				vertex[i0][0] += p[0];
				vertex[i0][1] += p[1];
				vertex[i0][2] += p[2];
				degree[i0]++;

				i0 = index[k+1][i+1][j];
				MAT_SUM(Q[i0], Q_tmp);
				vertex[i0][0] += p[0];
				vertex[i0][1] += p[1];
				vertex[i0][2] += p[2];
				degree[i0]++;
			}
	
	//FOR SVD
			
  vnl_matrix< float > A( 3, 3, 0. ); 
  vnl_vector<float> b( 3, 0. );

	for(int i=0; i<vertex_N; i++){
		if(degree[i] == 0)
			continue;
		vertex[i][0] /= degree[i];
		vertex[i][1] /= degree[i];
		vertex[i][2] /= degree[i];
continue;
		A[0][0] = (float)Q[i][0];
		A[1][0] = A[0][1] = (float)Q[i][1];
		A[2][0] = A[0][2] = (float)Q[i][2];
		A[1][1] = (float)Q[i][3];
		A[1][2] = A[2][1] = (float)Q[i][4];
		A[2][2] = (float)Q[i][5];

		float Av[3];
		MAT_BY_VEC(Av, Q[i], vertex[i]);
		b[1] = -(float)Q[i][6] - Av[0];
		b[2] = -(float)Q[i][7] - Av[1];
		b[3] = -(float)Q[i][8] - Av[2];

    vnl_svd<float> svd( A );
    svd.zero_out_absolute( tau );

    vnl_vector< float > x = svd.solve( b );

		if(fabs(x[1]) > spaceX || fabs(x[2]) > spaceY || fabs(x[3]) > spaceZ)
			continue;

		mesh->vertex[i][0] += x[1];
		mesh->vertex[i][1] += x[2];
		mesh->vertex[i][2] += x[3];
	}

	for(int i=0; i<dimZ+1; i++){
		for(int j=0; j<dimY+1; j++){
			delete[] isIn[i][j];
			delete[] index[i][j];
			delete[] value[i][j];
		}
		delete[] isIn[i];
		delete[] index[i];
		delete[] value[i];
	}
	delete[] isIn;
	delete[] index;
	delete[] Q;
	delete[] value;

	return mesh;
}

void Polygonizer::searchZero(float p[], float p1[], float p2[], float f1, float f2, float e)
{
	//regula falsa
	float p3[3], f3;
	if(f1 > 0){
		p3[0] = p1[0];
		p3[1] = p1[1];
		p3[2] = p1[2];

		p1[0] = p2[0];
		p1[1] = p2[1];
		p1[2] = p2[2];

		p2[0] = p3[0];
		p2[1] = p3[1];
		p2[2] = p3[2];

		float swap = f1;
		f1 = f2;
		f2 = swap;
	}

	float dt[3];
	dt[0] = (p2[0] - p1[0]);
	dt[1] = (p2[1] - p1[1]);
	dt[2] = (p2[2] - p1[2]);
	for(int j=0; j<5; j++){
		if(PolygonalMesh::LENGTH(dt) < 0.0001){
			p[0] = p3[0];
			p[1] = p3[1];
			p[2] = p3[2];
			return;
		}

		p3[0] = p1[0] + dt[0]*f1/(f1-f2);
		p3[1] = p1[1] + dt[1]*f1/(f1-f2);
		p3[2] = p1[2] + dt[2]*f1/(f1-f2);
		f3 = func->value(p3[0], p3[1], p3[2]);

		if(fabs(f3) < 0.000001){
			p[0] = p3[0];
			p[1] = p3[1];
			p[2] = p3[2];
			return;
		}
		if(f3 < 0.0){
			p1[0] = p3[0];
			p1[1] = p3[1];
			p1[2] = p3[2];

			f1 = f3;
		}
		else{
			p2[0] = p3[0];
			p2[1] = p3[1];
			p2[2] = p3[2];

			f2 = f3;
		}
		dt[0] = (p2[0] - p1[0]);
		dt[1] = (p2[1] - p1[1]);
		dt[2] = (p2[2] - p1[2]);
	}

	//Bisection method
	for(int j=0; j<5; j++){
		if(PolygonalMesh::DIST(p1, p2) < 0.000001){
			p[0] = p3[0];
			p[1] = p3[1];
			p[2] = p3[2];
			return;
		}

		p3[0] = 0.5f*(p1[0] + p2[0]);
		p3[1] = 0.5f*(p1[1] + p2[1]);
		p3[2] = 0.5f*(p1[2] + p2[2]);
		f3 = func->value(p3[0], p3[1], p3[2]);
		if(fabs(f3) < 0.000001){
			p[0] = p3[0];
			p[1] = p3[1];
			p[2] = p3[2];
			return;
		}
		else if(f1*f3 >= 0){
			p1[0] = p3[0];
			p1[1] = p3[1];
			p1[2] = p3[2];
			f1 = f3;
		}
		else{
			p2[0] = p3[0];
			p2[1] = p3[1];
			p2[2] = p3[2];
			f2 = f3;
		}
	}

	p[0] = p3[0];
	p[1] = p3[1];
	p[2] = p3[2];
}

PolygonalMesh* Polygonizer::computeSurfaceNetOctTree(float tol, int min, int max)
{
	OctTreeP oct;
	oct.func = func;
	oct.originX = originX;
	oct.originY = originY;
	oct.originZ = originZ;
	oct.sizeX = dimX*spaceX;
	oct.sizeY = dimY*spaceY;
	oct.sizeZ = dimZ*spaceZ;
	oct.constructStandard(min, max);

	int vertex_N, face_N;
	oct.countVertexAndFace(vertex_N, face_N);
	PolygonalMesh* mesh = new PolygonalMesh;
	mesh->setVertexCount(vertex_N);
	mesh->setFaceCount(face_N);
	for(int i=0; i<face_N; i++)
		mesh->setPolygonCount(i, 4);
	
	double (*Q)[10] = new double[vertex_N][10];
	for(int i=0; i<vertex_N; i++)
		MAT_INIT(Q[i]);

	oct.simplify(Q, tol);
	oct.polygonize2(mesh->vertex, mesh->face);
	//return mesh;

	//oct.simplify(Q, tol);
	//oct.polygonize2(mesh->vertex, mesh->face);

	//FOR SVD
			
  vnl_matrix< float > A( 3, 3, 0. );
  
  vnl_vector<float> b( 3, 0. ), x( 3, 0. );
	for(int i=0; i<vertex_N; i++){
		A[0][0] = Q[i][0];
		A[1][0] = A[0][1] = Q[i][1];
		A[2][0] = A[0][2] = Q[i][2];
		A[1][1] = Q[i][3];
		A[1][2] = A[2][1] = Q[i][4];
		A[2][2] = Q[i][5];

		float tmp[3];
		MAT_BY_VEC(tmp, Q[i], mesh->vertex[i]);
		
		b[1] = -Q[i][6] - tmp[0];
		b[2] = -Q[i][7] - tmp[1];
		b[3] = -Q[i][8] - tmp[2];
		
        vnl_svd< float > svd( A );
        svd.zero_out_relative( 0.025 );
		x = svd.solve( b );

		mesh->vertex[i][0] += x[1];
		mesh->vertex[i][1] += x[2];
		mesh->vertex[i][2] += x[3];
	}

	delete[] Q;
	return mesh;
}

PolygonalMesh* Polygonizer::bloomenthal(float size, float o[3], float box[6])
{
	/*
	int bounds[6];
	bounds[0] = (int)((box[0] - o[0])/size) +1;
	bounds[1] = (int)((box[1] - o[0])/size) -1;
	bounds[2] = (int)((box[2] - o[1])/size) +1;
	bounds[3] = (int)((box[3] - o[1])/size) -1;
	bounds[4] = (int)((box[4] - o[2])/size) +1;
	bounds[5] = (int)((box[5] - o[2])/size) -1;
	*/

	float s = (box[0] - box[1])*0.1f;
	box[0] += s; box[1] -= s;
	s = (box[2] - box[3])*0.1f;
	box[2] += s; box[3] -= s;
	s = (box[4] - box[5])*0.1f;
	box[4] += s; box[5] -= s;
	BLOOMENTHAL::Bloomenthal bloom;
	bloom.func = func;
	PolygonalMesh* mesh = new PolygonalMesh;
	bloom.polygonize(mesh, size, box, o[0], o[1], o[2], 1); 
	return mesh;
}
