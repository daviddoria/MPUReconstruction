// PolygonalMesh.cpp: PolygonalMesh 
//
//////////////////////////////////////////////////////////////////////

//#include "stdafx.h"
#include "PolygonalMesh.h"
#include <stdio.h>

//////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////

PolygonalMesh::PolygonalMesh()
{
	face_N = 0;
	vertex_N = 0;

	vertex = NULL;
	face = NULL;

	//vertex_link_v = NULL;
	//vertex_link_f = NULL;

	//degree_v = NULL;
	//degree_f = NULL;

	//isBound = NULL;

	normal_f = NULL;
	normal = NULL;

	center = NULL;

	value = NULL;
}

PolygonalMesh::~PolygonalMesh()
{
	if(vertex != NULL)
		delete[] vertex;
	if(face != NULL){
		for(int i=0; i<face_N; i++){
			if(poly_N != 0)
				delete[] face[i];
		}
		delete[] face;
	}
	delete[] poly_N;
	if(normal != NULL)
		delete[] normal;
	if(normal_f != NULL)
		delete[] normal_f;

	if(center != NULL)
		delete[] center;

	if(value != NULL)
		delete[] value;
}

void PolygonalMesh::setVertexCount(int vertex_N)
{
	if(vertex != NULL)
		delete[] vertex;
	this->vertex_N = vertex_N;
	vertex = new float[vertex_N][3];
}

void PolygonalMesh::setFaceCount(int face_N)
{
	if(face != NULL){
		delete[] face;
		delete[] poly_N;
	}
	this->face_N = face_N;
	face = new int*[face_N];
	poly_N = new int[face_N];
}

void PolygonalMesh::setPolygonCount(int index, int n)
{
	poly_N[index] = n;
	face[index] = new int[n];
	face[index][0] = -1;
}

void PolygonalMesh::computeFaceNormal()
{
	if(normal_f == NULL)
		normal_f = new float[face_N][3];
	for(int i=0; i<face_N; i++){
		int n = poly_N[i];
		int *f = face[i];
		if(f[0] < 0)
			continue;
		double nor[3];
		nor[0] = nor[1] = nor[2] = 0;
		for(int j=0; j<n; j++){
			int i1 = f[j];
			int i2 = f[(j+1)%n];
			nor[0] += (vertex[i1][1] - vertex[i2][1])*(vertex[i1][2] + vertex[i2][2]);
			nor[1] += (vertex[i1][2] - vertex[i2][2])*(vertex[i1][0] + vertex[i2][0]);
			nor[2] += (vertex[i1][0] - vertex[i2][0])*(vertex[i1][1] + vertex[i2][1]);
		}
		double len = LENGTH(nor);
		if((float)len != 0){
			normal_f[i][0] = (float)(nor[0]/len);
			normal_f[i][1] = (float)(nor[1]/len);
			normal_f[i][2] = (float)(nor[2]/len);
		}
		else{
			normal_f[i][0] = normal_f[i][1] = normal_f[i][2] = 0;
		}
	}
}

void PolygonalMesh::computeNormal()
{
	if(normal == NULL)
		normal = new float[vertex_N][3];
	for(int i=0; i<vertex_N; i++)
		normal[i][0] = normal[i][1] = normal[i][2] = 0;

	for(int i=0; i<face_N; i++){
		int n = poly_N[i];
		int *f = face[i];
		double nor[3];
		nor[0] = nor[1] = nor[2] = 0;
		if(f[0] < 0)
			continue;
		for(int j=0; j<n; j++){
			int i1 = f[j];
			int i2 = f[(j+1)%n];
			nor[0] += (vertex[i1][1] - vertex[i2][1])*(vertex[i1][2] + vertex[i2][2]);
			nor[1] += (vertex[i1][2] - vertex[i2][2])*(vertex[i1][0] + vertex[i2][0]);
			nor[2] += (vertex[i1][0] - vertex[i2][0])*(vertex[i1][1] + vertex[i2][1]);
		}
		for(int j=0; j<n; j++){
			int v = face[i][j];
			normal[v][0] += (float)nor[0];
			normal[v][1] += (float)nor[1];
			normal[v][2] += (float)nor[2];
		}
	}

	for(int i=0; i<vertex_N; i++){
		double len = LENGTH(normal[i]);
		if((float)len != 0){
			normal[i][0] /= (float)len;
			normal[i][1] /= (float)len;
			normal[i][2] /= (float)len;
		}
	}
}

void PolygonalMesh::computeCenter()
{
	if(center == NULL)
		center = new float[face_N][3];
	for(int i=0; i<face_N; i++){
		center[i][0] = center[i][1] = center[i][2] = 0;
		int m = poly_N[i];
		for(int j=0; j<m; j++){
			center[i][0] += vertex[face[i][j]][0];
			center[i][1] += vertex[face[i][j]][1];
			center[i][2] += vertex[face[i][j]][2];
		}
		center[i][0] /= m;
		center[i][1] /= m;
		center[i][2] /= m;
	}
}
