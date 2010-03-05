// OctTree.h: OctTree 
//
//////////////////////////////////////////////////////////////////////

#include "ImplicitFunction.h"

#define BD 0
#define FD 1
#define RD 2
#define LD 3
#define BU 4
#define FU 5
#define RU 6
#define LU 7
#define BL 8
#define FL 9
#define BR 10
#define FR 11

#define DOWN 0
#define UP 1
#define BACK 2
#define FRONT 3
#define RIGHT 5
#define LEFT 4

class OctTreeP  
{
	struct Cell{
		Cell  *parent;
		Cell *child[8];
		int level;
		float o[3];
		bool isLeaf;
			
		float value[8];

		int vertex_id;
	};

public:
	void edgeProcS(Cell *c1, Cell *c2, Cell *c3, Cell *c4, int edge);
	void faceProcS(Cell *c1, Cell *c2, int face);
	void cellProcS(Cell *c);
	void polygonize2(float (*vertex)[3], int **face);
	void computeOptP(float P[3], double Q[10], Cell* c);
	void simplifyCell(Cell* c, float tol);
	void edgeProcQ(Cell *c1, Cell *c2, Cell *c3, Cell *c4, int edge);
	void faceProcQ(Cell *c1, Cell *c2, int face);
	void cellProcQ(Cell *c);
	void simplify(double (*Q)[10], float tol);
	void edgeProc2(Cell *c1, Cell *c2, Cell *c3, Cell *c4, int edge);
	void faceProc2(Cell *c1, Cell *c2, int face);
	void cellProc2(Cell* c);
	void countVertexAndFace(int &vertex_N, int &face_N);
	void bisection(float p[], float start[], float end[], float f1, float f2);
	void polygonize(float (*vertex)[3], int **face, double (*Q)[10]);
	void edgeProc(Cell* c1, Cell* c2, Cell* c3, Cell* c4, int edge);
	void faceProc(Cell* c1, Cell* c2, int face);
	void cellProc(Cell* c);
	int getNeighborVertex(Cell* c, int dire);
	void incrementZeroCrossEdges(Cell* c);
	int countZeroCrossEdges();
	inline void CELL_P(float p[3], Cell* c);
	inline void FACE_P(float p[3], Cell* c, int i);
	inline void EDGE_P(float p[3], Cell *c, int i);
	inline void CONNER_P(float p[3], Cell* c, int i);
	void subdivideADF(Cell *c, float tol, int min, int max);
	void setVerticeAtCellCenter(Cell* c);
	void setVerticeAtCellCenter();
	void incrementCountAtBoundCell(Cell* c);
	int countBoundCell();
	void generateChild(Cell* c);
	void subdivideStandard(Cell* c, int min, int max);
	void constructStandard(int min, int max);
	bool checkBound(Cell* c);
	bool checkTol(Cell* c, float tol);
	void setValue(Cell* c);
	Cell* initCell(int level, float o[3]);
	void freeCell(Cell* c);
	void constructADF(float tol, int min, int max);
	OctTreeP();
	virtual ~OctTreeP();

	float sizeX;
	float sizeY;
	float sizeZ;

	float originX;
	float originY;
	float originZ;

	ImplicitFunction* func;

	Cell* root;

	float (*vertex)[3];
	int **face;
	double (*Q)[10];

private:
	int counter;
	int counter2;

	static inline double DET(double A[10]){
		return A[0]*A[3]*A[5] + 2.0*A[1]*A[4]*A[2] 
			-A[2]*A[2]*A[3] - A[1]*A[1]*A[5] - A[4]*A[4]*A[0];
	}

	static inline void MATRIX(double A[10], double n[3], double d){
		A[0] = n[0]*n[0];
		A[1] = n[0]*n[1];
		A[2] = n[0]*n[2];
		A[3] = n[1]*n[1];
		A[4] = n[1]*n[2];
		A[5] = n[2]*n[2];
		A[6] = d*n[0];
		A[7] = d*n[1];
		A[8] = d*n[2];
		A[9] = d*d;
	}

	static inline void MAT_TIMES(double A[10], double k){
		A[0] *= k;
		A[1] *= k;
		A[2] *= k;
		A[3] *= k;
		A[4] *= k;
		A[5] *= k;
		A[6] *= k;
		A[7] *= k;
		A[8] *= k;
		A[9] *= k;
	}

	static inline void MAT_SUM(double B[6], double A[6]){
		B[0] += A[0];
		B[1] += A[1];
		B[2] += A[2];
		B[3] += A[3];
		B[4] += A[4];
		B[5] += A[5];
		B[6] += A[6];
		B[7] += A[7];
		B[8] += A[8];
		B[9] += A[9];
	}

	static inline void MAT_BY_VEC(double v[3], double A[10], double b[3]){
		v[0] = A[0]*b[0] + A[1]*b[1] + A[2]*b[2];
		v[1] = A[1]*b[0] + A[3]*b[1] + A[4]*b[2];
		v[2] = A[2]*b[0] + A[4]*b[1] + A[5]*b[2];
	}

	static inline void MAT_BY_VEC(float v[3], double A[10], float b[3]){
		v[0] = (float)(A[0]*b[0] + A[1]*b[1] + A[2]*b[2]);
		v[1] = (float)(A[1]*b[0] + A[3]*b[1] + A[4]*b[2]);
		v[2] = (float)(A[2]*b[0] + A[4]*b[1] + A[5]*b[2]);
	}

	static inline void MAT_INIT(double A[10]){
		A[0] = A[1] = A[2] = A[3] = A[4] = A[5] = A[6] = A[7] = A[8] = A[9] = 0;
	}

	static inline void MAT_PLUS(double C[10], double A[10], double B[10]){
		C[0] = A[0] + B[0];
		C[1] = A[1] + B[1];
		C[2] = A[2] + B[2];
		C[3] = A[3] + B[3];
		C[4] = A[4] + B[4];
		C[5] = A[5] + B[5];
		C[6] = A[6] + B[6];
		C[7] = A[7] + B[7];
		C[8] = A[8] + B[8];
		C[9] = A[9] + B[9];
	}

	static inline void MAT_COPY(double A[10], double B[10]){
		A[0] = B[0];
		A[1] = B[1];
		A[2] = B[2];
		A[3] = B[3];
		A[4] = B[4];
		A[5] = B[5];
		A[6] = B[6];
		A[7] = B[7];
		A[8] = B[8];
		A[9] = B[9];
	}

	static inline double Q_ERR(double A[10], float v[3]){
		return v[0]*(A[0]*v[0] + A[1]*v[1] + A[2]*v[2]) +
			   v[1]*(A[1]*v[0] + A[3]*v[1] + A[4]*v[2]) +
			   v[2]*(A[2]*v[0] + A[4]*v[1] + A[5]*v[2]) +
			   2.0*(A[6]*v[0] + A[7]*v[1] + A[8]*v[2]) + 
			   A[9];
	}

	static inline double Q_ERR(double A[10], double v[3]){
		return v[0]*(A[0]*v[0] + A[1]*v[1] + A[2]*v[2]) +
			   v[1]*(A[1]*v[0] + A[3]*v[1] + A[4]*v[2]) +
			   v[2]*(A[2]*v[0] + A[4]*v[1] + A[5]*v[2]) +
			   2.0*(A[6]*v[0] + A[7]*v[1] + A[8]*v[2]) + 
			   A[9];
	}
};
