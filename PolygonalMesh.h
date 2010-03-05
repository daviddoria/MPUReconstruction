// PolygonalMesh.h: PolygonalMesh 
//////////////////////////////////////////////////////////////////////

#ifndef POLYGONALMESH_H
#define POLYGONALMESH_H

#include <math.h>

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class PolygonalMesh  
{
public:
	//The number of vertex and face (File data)
	int face_N, vertex_N;

	//Vertex coordinates
	float (*vertex)[3];
	//Face index
	int **face;
	int *poly_N;

	int *degree_f;

	float (*normal_f)[3];
	float (*normal)[3];
	float (*center)[3];

	//int **vertex_link_f;

	//Is the vertex specified a index boundary?
	//value = 0 if not boundary 
	//velue = index of adjacent boundary face
	//int *isBound;

	float *value;

public:
	void computeCenter();
	void computeNormal();
	void computeFaceNormal();
	void setPolygonCount(int index, int n);
	void setFaceCount(int face_N);
	void setVertexCount(int vertex_N);
	PolygonalMesh();
	virtual ~PolygonalMesh();


	static inline double LENGTH(float v[3]){
		return sqrt((double)v[0]*(double)v[0] + (double)v[1]*(double)v[1] + (double)v[2]*(double)v[2]);
	}

	static inline double LENGTH(double v[3]){
		return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
	}

	static inline void VEC(float v[3], float p1[3], float p2[3]){
		v[0] = p2[0]-p1[0];
		v[1] = p2[1]-p1[1];
		v[2] = p2[2]-p1[2];
	}

	static inline void VEC(double v[3], float p1[3], float p2[3]){
		v[0] = p2[0]-p1[0];
		v[1] = p2[1]-p1[1];
		v[2] = p2[2]-p1[2];
	}

	static inline void VEC(double v[3], double p1[3], double p2[3]){
		v[0] = p2[0]-p1[0];
		v[1] = p2[1]-p1[1];
		v[2] = p2[2]-p1[2];
	}

	static inline void VEC(double v[3], double p1[3], float p2[3]){
		v[0] = (double)p2[0]-p1[0];
		v[1] = (double)p2[1]-p1[1];
		v[2] = (double)p2[2]-p1[2];
	}

	static inline void VEC(double v[3], float p1[3], double p2[3]){
		v[0] = p2[0]-(double)p1[0];
		v[1] = p2[1]-(double)p1[1];
		v[2] = p2[2]-(double)p1[2];
	}
	
	static inline void TIMES(float kv[3], float k, float v[3]){
		kv[0] = k*v[0];
		kv[1] = k*v[1];
		kv[2] = k*v[2];
	}

	static inline double DIST(float p1[3], float p2[3]){
		float v[3];
		VEC(v,p1,p2);
		return LENGTH(v);
	}

	static inline double DIST(double p1[3], double p2[3]){
		double v[3];
		VEC(v,p1,p2);
		return LENGTH(v);
	}

	static inline double DIST2(float p1[3], float p2[3]){
		float v[3];
		VEC(v,p1,p2);
		return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
	}

	static inline double DIST2(double p1[3], double p2[3]){
		double v[3];
		VEC(v,p1,p2);
		return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
	}

	static inline double DIST2(double p1[3], float p2[3]){
		double v[3];
		VEC(v,p1,p2);
		return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
	}

	static inline double DOT(float v1[3], float v2[3]){
		return (double)v1[0]*(double)v2[0] + (double)v1[1]*(double)v2[1] + (double)v1[2]*(double)v2[2];
	}

	static inline double DOT(double v1[3], float v2[3]){
		return v1[0]*(double)v2[0] + v1[1]*(double)v2[1] + v1[2]*(double)v2[2];
	}

	static inline double DOT(float v1[3], double v2[3]){
		return (double)v1[0]*v2[0] + (double)v1[1]*v2[1] + (double)v1[2]*v2[2];
	}

	static inline double DOT(double v1[3], double v2[3]){
		return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
	}

	static inline void CROSS(float n[3], float v1[3], float v2[3]){
		n[0] = v1[1]*v2[2] - v1[2]*v2[1];
		n[1] = v1[2]*v2[0] - v1[0]*v2[2];
		n[2] = v1[0]*v2[1] - v1[1]*v2[0];
	}

	static inline void CROSS(double n[3], double v1[3], double v2[3]){
		n[0] = v1[1]*v2[2] - v1[2]*v2[1];
		n[1] = v1[2]*v2[0] - v1[0]*v2[2];
		n[2] = v1[0]*v2[1] - v1[1]*v2[0];
	}

	static inline void CROSS(double n[3], float v1[3], float v2[3]){
		n[0] = (double)v1[1]*(double)v2[2] - (double)v1[2]*(double)v2[1];
		n[1] = (double)v1[2]*(double)v2[0] - (double)v1[0]*(double)v2[2];
		n[2] = (double)v1[0]*(double)v2[1] - (double)v1[1]*(double)v2[0];
	}

	static inline void CROSS(double n[3], double v1[3], float v2[3]){
		n[0] = v1[1]*(double)v2[2] - v1[2]*(double)v2[1];
		n[1] = v1[2]*(double)v2[0] - v1[0]*(double)v2[2];
		n[2] = v1[0]*(double)v2[1] - v1[1]*(double)v2[0];
	}

	static inline void CROSS(double n[3], float v1[3], double v2[3]){
		n[0] = (double)v1[1]*v2[2] - (double)v1[2]*v2[1];
		n[1] = (double)v1[2]*v2[0] - (double)v1[0]*v2[2];
		n[2] = (double)v1[0]*v2[1] - (double)v1[1]*v2[0];
	}

	static inline void CROSS(float n[3], double v1[3], double v2[3]){
		n[0] = (float)(v1[1]*v2[2] - v1[2]*v2[1]);
		n[1] = (float)(v1[2]*v2[0] - v1[0]*v2[2]);
		n[2] = (float)(v1[0]*v2[1] - v1[1]*v2[0]);
	}

	static inline void CROSS(float n[3], double v1[3], float v2[3]){
		n[0] = (float)(v1[1]*v2[2] - v1[2]*v2[1]);
		n[1] = (float)(v1[2]*v2[0] - v1[0]*v2[2]);
		n[2] = (float)(v1[0]*v2[1] - v1[1]*v2[0]);
	}

	static inline void CROSS(float n[3], float v1[3], double v2[3]){
		n[0] = (float)(v1[1]*v2[2] - v1[2]*v2[1]);
		n[1] = (float)(v1[2]*v2[0] - v1[0]*v2[2]);
		n[2] = (float)(v1[0]*v2[1] - v1[1]*v2[0]);
	}

	static inline double AREA(float p1[3], float p2[3], float p3[3])
	{
		double n[3];
		float v1[3], v2[3];
		VEC(v1, p2, p1);
		VEC(v2, p3, p1);
		CROSS(n, v1, v2);

		return 0.5*LENGTH(n);
	}

	static inline double AREA(double p1[3], double p2[3], double p3[3])
	{
		double n[3];
		double v1[3], v2[3];
		VEC(v1, p2, p1);
		VEC(v2, p3, p1);
		CROSS(n, v1, v2);

		return 0.5*LENGTH(n);
	}

	static inline double AREA(float p1[3], float p2[3], double p3[3])
	{
		double n[3];
		float v1[3];
		double v2[3];
		VEC(v1, p2, p1);
		VEC(v2, p3, p1);
		CROSS(n, v1, v2);

		return 0.5*LENGTH(n);
	}
};

#endif