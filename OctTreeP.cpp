// OctTreeP.cpp: OctTreeP 
//
//////////////////////////////////////////////////////////////////////

//#include "stdafx.h"
#include "OctTreeP.h"
#include "math.h"
#include <vnl/algo/vnl_svd.h>
#include "PolygonalMesh.h"

static int EDGE[12][2] = {{0,1}, {2,3}, {1,3}, {0,2},
						{4,5}, {6,7}, {5,7}, {4,6},
						{0,4}, {2,6}, {1,5}, {3,7}};
static int FACE[6][4] = {{0,1,2,3}, {4,5,6,7}, {0,1,4,5}, {2,3,6,7}, {0,2,4,6}, {1,3,5,7}};

OctTreeP::OctTreeP()  
{
	root = NULL;
}

OctTreeP::~OctTreeP()
{
	freeCell(root);
}

void OctTreeP::constructADF(float tol, int min, int max)
{
	if(root != NULL)
		freeCell(root);
	float o[3];
	o[0] = originX;
	o[1] = originY;
	o[2] = originZ;
	root = initCell(0, o);
	setValue(root);
	subdivideADF(root, tol, min, max);
}

void OctTreeP::freeCell(Cell *c)
{
	if(!c->isLeaf){
		for(int i=0; i<8; i++)
			freeCell(c->child[i]);
	}
	delete c;
}

OctTreeP::Cell* OctTreeP::initCell(int level, float o[3])
{
	Cell* c = new Cell;
	c->isLeaf = true;
	c->level = level;
	c->o[0] = o[0];
	c->o[1] = o[1];
	c->o[2] = o[2];
	c->vertex_id = -1;
	return c;
}

void OctTreeP::setValue(Cell *c)
{
	int l = c->level;
	float spaceX = (float)(sizeX*pow(2.0, -l));
	float spaceY = (float)(sizeY*pow(2.0, -l));
	float spaceZ = (float)(sizeZ*pow(2.0, -l));
	for(int i=0; i<8; i++){
		float x[3];
		if(i%2 == 0)
			x[0] = c->o[0] ;
		else
			x[0] = c->o[0] + spaceX;
		if(i%4 < 2)
			x[1] = c->o[1];
		else
			x[1] = c->o[1] + spaceY;
		if(i < 4)
			x[2] = c->o[2];
		else
			x[2] = c->o[2] + spaceZ;
		c->value[i] = func->value(x[0], x[1], x[2]);
	}
/*
	x[0] = c->o[0];
	x[1] = c->o[1];
	x[2] = c->o[2];
	c->value[0] = func->value(x);
	x[0] += spaceX;
	c->value[1] = func->value(x);
	x[1] += spaceY;
	c->value[2] = func->value(x);
	x[0] = c->o[0];
	c->value[3] = func->value(x);
	x[1] = c->o[1];
	x[2] += spaceZ;
	c->value[4] = func->value(x);
	x[0] += spaceX;
	c->value[5] = func->value(x);
	x[1] += spaceY;
	c->value[6] = func->value(x);
	x[0] = c->o[0];
	c->value[7] = func->value(x);
	*/
}

bool OctTreeP::checkTol(Cell *c, float tol)
{
	//underconstruction
	return false;
}

bool OctTreeP::checkBound(Cell *c)
{
	bool flag = (c->value[0] >= 0);
	for(int i=1; i<8; i++){
		if(flag != (c->value[i] >= 0))
			return true;
	}
	return false;
}

void OctTreeP::constructStandard(int min, int max)
{
	if(root != NULL)
		freeCell(root);
	float o[3];
	o[0] = originX;
	o[1] = originY;
	o[2] = originZ;
	root = initCell(0, o);
	setValue(root);
	subdivideStandard(root, min, max);
}

void OctTreeP::subdivideStandard(Cell *c, int min, int max)
{
	if(c->level == max)
		return;
	if(c->level < min || checkBound(c)){
		generateChild(c);
		for(int i=0; i<8; i++)
			subdivideStandard(c->child[i], min, max);
	}
}

void OctTreeP::generateChild(Cell *c)
{
	c->isLeaf = false;
	int l = c->level;
	float spaceX = (float)(sizeX*pow(2.0, -l));
	float spaceY = (float)(sizeY*pow(2.0, -l));
	float spaceZ = (float)(sizeZ*pow(2.0, -l));
	for(int i=0; i<8; i++){
		float o[3];
		if(i%2 == 0)
			o[0] = c->o[0] ;
		else
			o[0] = c->o[0] + 0.5f*spaceX;
		if(i%4 < 2)
			o[1] = c->o[1];
		else
			o[1] = c->o[1] + 0.5f*spaceY;
		if(i < 4)
			o[2] = c->o[2];
		else
			o[2] = c->o[2] + 0.5f*spaceZ;
		c->child[i] = initCell(l+1, o);
		c->child[i]->parent = c;
		//setValue(c->child[i]);
	}
	for(int i=0; i<8; i++)
		c->child[i]->value[i] = c->value[i];

	float x[3], f;
	
	for(int i=0; i<12; i++){
		EDGE_P(x, c, i);
		f = func->value(x[0], x[1], x[2]);
		c->child[EDGE[i][0]]->value[EDGE[i][1]] = f;
		c->child[EDGE[i][1]]->value[EDGE[i][0]] = f;
	}
	for(int i=0; i<6; i++){
		FACE_P(x, c, i);
		f  = func->value(x[0], x[1], x[2]);
		c->child[FACE[i][0]]->value[FACE[i][3]] = f;
		c->child[FACE[i][3]]->value[FACE[i][0]] = f;
		c->child[FACE[i][1]]->value[FACE[i][2]] = f;
		c->child[FACE[i][2]]->value[FACE[i][1]] = f;
	}
	CELL_P(x, c);
	f = func->value(x[0], x[1], x[2]);
	for(int i=0; i<8; i++)
		c->child[i]->value[7-i] = f;
}

int OctTreeP::countBoundCell()
{
	counter = 0;
	incrementCountAtBoundCell(root);
	return counter;
}

void OctTreeP::incrementCountAtBoundCell(Cell *c)
{
	if(!c->isLeaf){
		for(int i=0; i<8; i++)
			incrementCountAtBoundCell(c->child[i]);
	}
	else if(checkBound(c)){
		counter++;
	}
	
}

void OctTreeP::setVerticeAtCellCenter()
{
	setVerticeAtCellCenter(root);
}

void OctTreeP::setVerticeAtCellCenter(Cell *c)
{
	if(!c->isLeaf){
		for(int i=0; i<8; i++)
			setVerticeAtCellCenter(c->child[i]);
	}
	if(c->vertex_id >= 0)
		CELL_P(vertex[c->vertex_id], c);
}

void OctTreeP::subdivideADF(Cell *c, float tol, int min, int max)
{
	if(c->level == max)
		return;
	if(c->level < min){
		generateChild(c);
		for(int i=0; i<8; i++)
			subdivideADF(c->child[i], tol, min, max);
		return;
	}
	if(!checkBound(c))
		return;

	bool flag = false;
	float *vc = c->value;
	float vf[12];
	for(int i=0; i<12; i++){
		float l = 0.5f*(vc[EDGE[i][0]] + vc[EDGE[i][1]]);
		float x[3];
		EDGE_P(x, c, i);
		vf[i] = func->value(x[0], x[1], x[2]);
		if(!flag){
			flag = (fabs(l - vf[i]) > tol);
		}
	}
	for(int i=0; i<6; i++){
		float l = 0.25f*(vc[FACE[i][0]] + vc[FACE[i][1]] + vc[FACE[i][2]] + vc[FACE[i][3]]);
		float x[3];
		FACE_P(x, c, i);
		vf[i] = func->value(x[0], x[1], x[2]);
		if(!flag){
			flag = (fabs(l - vf[i]) > tol);
		}
	}
	float l = 0.125f*(vc[0] + vc[1] + vc[2] + vc[3] + vc[4] + vc[5] + vc[6] + vc[7]);
	float x[3];
	CELL_P(x, c);
	vf[0] = func->value(x[0], x[1], x[2]);
	if(!flag){
		flag = (fabs(l - vf[0]) > tol);
	}
	if(flag){
		generateChild(c);
		for(int i=0; i<8; i++)
			subdivideADF(c->child[i], tol, min, max);
	}
}

inline void OctTreeP::CONNER_P(float p[3], Cell *c, int i)
{
	p[0] = c->o[0];
	p[1] = c->o[1];
	p[2] = c->o[2];
	if(i > 3)
		p[2] += (float)(sizeZ*pow(2.0, -c->level));
	if(i%4 > 1)
		p[1] += (float)(sizeY*pow(2.0, -c->level));
	if(i%2 == 1)
		p[0] += (float)(sizeX*pow(2.0, -c->level));
}

inline void OctTreeP::EDGE_P(float p[], Cell *c, int i)
{
	float p1[3], p2[3];
	CONNER_P(p1, c, EDGE[i][0]);
	CONNER_P(p2, c, EDGE[i][1]);
	p[0] = 0.5f*(p1[0] + p2[0]);
	p[1] = 0.5f*(p1[1] + p2[1]);
	p[2] = 0.5f*(p1[2] + p2[2]);
}

inline void OctTreeP::FACE_P(float p[], Cell *c, int i)
{
	float q[3];
	CONNER_P(p, c, FACE[i][0]);
	for(int j=1; j<4; j++){
		CONNER_P(q, c, FACE[i][j]);
		p[0] += q[0];
		p[1] += q[1];
		p[2] += q[2];
	}
	p[0] *= 0.25f;
	p[1] *= 0.25f;
	p[2] *= 0.25f;
}

inline void OctTreeP::CELL_P(float p[], Cell *c)
{
	p[0] = c->o[0] + (float)(sizeX*pow(2.0, -c->level-1));
	p[1] = c->o[1] + (float)(sizeY*pow(2.0, -c->level-1));
	p[2] = c->o[2] + (float)(sizeZ*pow(2.0, -c->level-1));
}

int OctTreeP::countZeroCrossEdges()
{
	counter = 0;
	incrementZeroCrossEdges(root);
	return counter;
}

void OctTreeP::incrementZeroCrossEdges(Cell *c)
{
	if(c->value[EDGE[RU][0]] * c->value[EDGE[RU][1]] <= 0)
		counter++;
	if(c->value[EDGE[LD][0]] * c->value[EDGE[LD][1]] <= 0)
		counter++;
	if(c->value[EDGE[FU][0]] * c->value[EDGE[FU][1]] <= 0)
		counter++;
	if(c->value[EDGE[BD][0]] * c->value[EDGE[BD][1]] <= 0)
		counter++;
	if(c->value[EDGE[FR][0]] * c->value[EDGE[FR][1]] <= 0)
		counter++;
	if(c->value[EDGE[BL][0]] * c->value[EDGE[BL][1]] <= 0)
		counter++;
	if(c->isLeaf)
		return;
	for(int i=0; i<8; i++)
		incrementZeroCrossEdges(c->child[i]);
}

int OctTreeP::getNeighborVertex(Cell *c, int dire)
{
	/*
	Cell* q = c;
	for(Cell* p = c->parent; p != root; p = p->parent){
		int index;
		for(int i=0; i<8; i++){
			if(p->child[i] == q)
				break;
		}
		index = i;
		bool isOK = false;
		switch(index){
		case 0:
			if(dire == UP)
				q = p->child[4];
				
				|| dire == FRONT || RIGHT)
				q = 
		}
		q = p;
	}
	*/
	return -1;
}

void OctTreeP::cellProc(Cell *c)
{
	if(c->isLeaf)
		return;
	
	for(int i=0; i<8; i++)
		cellProc(c->child[i]);

	faceProc(c->child[EDGE[0][0]], c->child[EDGE[0][1]], RIGHT);
	faceProc(c->child[EDGE[1][0]], c->child[EDGE[1][1]], RIGHT);
	faceProc(c->child[EDGE[2][0]], c->child[EDGE[2][1]], FRONT);
	faceProc(c->child[EDGE[3][0]], c->child[EDGE[3][1]], FRONT);
	faceProc(c->child[EDGE[4][0]], c->child[EDGE[4][1]], RIGHT);
	faceProc(c->child[EDGE[5][0]], c->child[EDGE[5][1]], RIGHT);
	faceProc(c->child[EDGE[6][0]], c->child[EDGE[6][1]], FRONT);
	faceProc(c->child[EDGE[7][0]], c->child[EDGE[7][1]], FRONT);
	faceProc(c->child[EDGE[8][0]], c->child[EDGE[8][1]], UP);
	faceProc(c->child[EDGE[9][0]], c->child[EDGE[9][1]], UP);
	faceProc(c->child[EDGE[10][0]], c->child[EDGE[10][1]], UP);
	faceProc(c->child[EDGE[11][0]], c->child[EDGE[11][1]], UP);

	edgeProc(c->child[FACE[0][0]], c->child[FACE[0][1]], c->child[FACE[0][2]], c->child[FACE[0][3]], FR);
	edgeProc(c->child[FACE[1][0]], c->child[FACE[1][1]], c->child[FACE[1][2]], c->child[FACE[1][3]], FR);
	edgeProc(c->child[FACE[2][0]], c->child[FACE[2][1]], c->child[FACE[2][2]], c->child[FACE[2][3]], RU);
	edgeProc(c->child[FACE[3][0]], c->child[FACE[3][1]], c->child[FACE[3][2]], c->child[FACE[3][3]], RU);
	edgeProc(c->child[FACE[4][0]], c->child[FACE[4][1]], c->child[FACE[4][2]], c->child[FACE[4][3]], FU);
	edgeProc(c->child[FACE[5][0]], c->child[FACE[5][1]], c->child[FACE[5][2]], c->child[FACE[5][3]], FU);
}

void OctTreeP::faceProc(Cell *c1, Cell *c2, int face)
{
	if(c1->isLeaf && c2->isLeaf)
		return;

	Cell *tmp1[4];
	if(c1->isLeaf){
		for(int i=0; i<4; i++)
			tmp1[i] = c1;
	}
	else{
		for(int i=0; i<4; i++)
			tmp1[i] = c1->child[FACE[face][i]];
	}

	Cell *tmp2[4];
	if(c2->isLeaf){
		for(int i=0; i<4; i++)
			tmp2[i] = c2;
	}
	else{
		for(int i=0; i<4; i++)
			tmp2[i] = c2->child[FACE[face-1][i]];
	}

	for(int i=0; i<4; i++)
		faceProc(tmp1[i], tmp2[i], face);

	if(face == RIGHT){
		edgeProc(tmp1[0], tmp2[0], tmp1[2], tmp2[2], RU);
		edgeProc(tmp1[1], tmp2[1], tmp1[3], tmp2[3], RU);
		edgeProc(tmp1[0], tmp2[0], tmp1[1], tmp2[1], FR);
		edgeProc(tmp1[2], tmp2[2], tmp1[3], tmp2[3], FR);
	}
	else if(face == UP){
		edgeProc(tmp1[0], tmp1[2], tmp2[0], tmp2[2], FU);
		//edgeProc(tmp1[0], tmp2[0], tmp1[2], tmp2[2], FU);
		//edgeProc(tmp1[1], tmp2[1], tmp1[3], tmp2[3], FU);
		edgeProc(tmp1[1], tmp1[3], tmp2[1], tmp2[3], FU);
		edgeProc(tmp1[0], tmp1[1], tmp2[0], tmp2[1], RU);
		edgeProc(tmp1[2], tmp1[3], tmp2[2], tmp2[3], RU);
	}
	else{
		edgeProc(tmp1[0], tmp2[0], tmp1[2], tmp2[2], FU);
		edgeProc(tmp1[1], tmp2[1], tmp1[3], tmp2[3], FU);
		edgeProc(tmp1[0], tmp1[1], tmp2[0], tmp2[1], FR);
		edgeProc(tmp1[2], tmp1[3], tmp2[2], tmp2[3], FR);
	}
}

void OctTreeP::edgeProc(Cell *c1, Cell *c2, Cell *c3, Cell *c4, int edge)
{
	if(c1->isLeaf && c2->isLeaf && c3->isLeaf && c4->isLeaf){
		int i1 = c1->vertex_id;
		int i2 = c2->vertex_id;
		int i3 = c3->vertex_id;
		int i4 = c4->vertex_id;
		if(i1 == -1 || i2 == -1 || i3 == -1 || i4== -1)
			return;
		int max = c1->level;
		int index = 1;
		if(max < c2->level){
			max = c2->level;
			index = 2;
		}
		if(max < c3->level){
			max = c3->level;
			index = 3;
		}
		if(max < c4->level){
			max = c4->level;
			index = 4;
		}

		float p1[3], p2[3], f1, f2;
		if(index == 1){
			if(edge == RU){
				f1 = c1->value[7];
				f2 = c1->value[5];
				CONNER_P(p1, c1, 7);
				CONNER_P(p2, c1, 5);
			}
			else if(edge == FU){
				f1 = c1->value[6];
				f2 = c1->value[7];
				CONNER_P(p1, c1, 6);
				CONNER_P(p2, c1, 7);
			}
			else{
				f1 = c1->value[3];
				f2 = c1->value[7];
				CONNER_P(p1, c1, 3);
				CONNER_P(p2, c1, 7);
			}
		}
		else if(index == 2){
			if(edge == RU){
				f1 = c2->value[6];
				f2 = c2->value[4];
				CONNER_P(p1, c2, 6);
				CONNER_P(p2, c2, 4);
			}
			else if(edge == FU){
				f1 = c2->value[4];
				f2 = c2->value[5];
				CONNER_P(p1, c2, 4);
				CONNER_P(p2, c2, 5);
			}
			else{
				f1 = c2->value[2];
				f2 = c2->value[6];
				CONNER_P(p1, c2, 2);
				CONNER_P(p2, c2, 6);
			}
		}
		else if(index == 3){
			if(edge == RU){
				f1 = c3->value[3];
				f2 = c3->value[1];
				CONNER_P(p1, c3, 3);
				CONNER_P(p2, c3, 1);
			}
			else if(edge == FU){
				f1 = c3->value[2];
				f2 = c3->value[3];
				CONNER_P(p1, c3, 2);
				CONNER_P(p2, c3, 3);
			}
			else{
				f1 = c3->value[1];
				f2 = c3->value[5];
				CONNER_P(p1, c3, 1);
				CONNER_P(p2, c3, 5);
			}
		}
		else{
			if(edge == RU){
				f1 = c4->value[2];
				f2 = c4->value[0];
				CONNER_P(p1, c4, 2);
				CONNER_P(p2, c4, 0);
			}
			else if(edge == FU){
				f1 = c4->value[0];
				f2 = c4->value[1];
				CONNER_P(p1, c4, 0);
				CONNER_P(p2, c4, 1);
			}
			else{
				f1 = c4->value[0];
				f2 = c4->value[4];
				CONNER_P(p1, c4, 0);
				CONNER_P(p2, c4, 4);
			}
		}
		
		if((f1 < 0 || f2 >= 0) && (f2 < 0 || f1 >= 0))
			return;

		//bisection
		float p[3];
		bisection(p, p1, p2, f1, f2);
		float g[3];
		func->gradient(g, p[0], p[1], p[2]);
		double len = PolygonalMesh::LENGTH(g);
		double nor[3];
		if((float)len != 0){
			nor[0] = g[0]/len;
			nor[1] = g[1]/len;
			nor[2] = g[2]/len;
		}
		else
			nor[0] = nor[1] = nor[2] = 0;

		float v[3];
		float center[3];
		CELL_P(center, c1);
		PolygonalMesh::VEC(v, center, p);
		double d = -PolygonalMesh::DOT(v, nor);
		double Q_tmp[10];
		MATRIX(Q_tmp, nor, d);
		MAT_SUM(Q[c1->vertex_id], Q_tmp);

		if(c2->vertex_id != c1->vertex_id){
			CELL_P(center, c2);
			PolygonalMesh::VEC(v, center, p);
			d = -PolygonalMesh::DOT(v, nor);
			MATRIX(Q_tmp, nor, d);
			MAT_SUM(Q[c2->vertex_id], Q_tmp);
		}

		if(c3->vertex_id != c1->vertex_id 
			&& c3->vertex_id != c2->vertex_id){
			CELL_P(center, c3);
			PolygonalMesh::VEC(v, center, p);
			d = -PolygonalMesh::DOT(v, nor);
			MATRIX(Q_tmp, nor, d);
			MAT_SUM(Q[c3->vertex_id], Q_tmp);
		}

		if(c4->vertex_id != c1->vertex_id 
			&& c4->vertex_id != c2->vertex_id
			&& c4->vertex_id != c3->vertex_id){
			CELL_P(center, c4);
			PolygonalMesh::VEC(v, center, p);
			d = -PolygonalMesh::DOT(v, nor);
			MATRIX(Q_tmp, nor, d);
			MAT_SUM(Q[c4->vertex_id], Q_tmp);
		}

		if(f1 > 0){
			face[counter][0] = c1->vertex_id;
			face[counter][1] = c2->vertex_id;
			face[counter][2] = c4->vertex_id;
			face[counter][3] = c3->vertex_id;
		}
		else{
			face[counter][0] = c3->vertex_id;
			face[counter][1] = c4->vertex_id;
			face[counter][2] = c2->vertex_id;
			face[counter][3] = c1->vertex_id;
		}
		counter++;

		return;
	}

	Cell *tmp1[4], *tmp2[4];
	if(edge == FU){
		if(c1->isLeaf){
			tmp1[0] = c1;
			tmp2[0] = c1;
		}
		else{
			tmp1[0] = c1->child[6];
			tmp2[0] = c1->child[7];
		}

		if(c2->isLeaf){
			tmp1[1] = c2;
			tmp2[1] = c2;
		}
		else{
			tmp1[1] = c2->child[4];
			tmp2[1] = c2->child[5];
		}

		if(c3->isLeaf){
			tmp1[2] = c3;
			tmp2[2] = c3;
		}
		else{
			tmp1[2] = c3->child[2];
			tmp2[2] = c3->child[3];
		}

		if(c4->isLeaf){
			tmp1[3] = c4;
			tmp2[3] = c4;
		}
		else{
			tmp1[3] = c4->child[0];
			tmp2[3] = c4->child[1];
		}
	}
	else if(edge == RU){
		if(c1->isLeaf){
			tmp1[0] = c1;
			tmp2[0] = c1;
		}
		else{
			tmp1[0] = c1->child[5];
			tmp2[0] = c1->child[7];
		}

		if(c2->isLeaf){
			tmp1[1] = c2;
			tmp2[1] = c2;
		}
		else{
			tmp1[1] = c2->child[4];
			tmp2[1] = c2->child[6];
		}

		if(c3->isLeaf){
			tmp1[2] = c3;
			tmp2[2] = c3;
		}
		else{
			tmp1[2] = c3->child[1];
			tmp2[2] = c3->child[3];
		}

		if(c4->isLeaf){
			tmp1[3] = c4;
			tmp2[3] = c4;
		}
		else{
			tmp1[3] = c4->child[0];
			tmp2[3] = c4->child[2];
		}
	}
	else{
		if(c1->isLeaf){
			tmp1[0] = c1;
			tmp2[0] = c1;
		}
		else{
			tmp1[0] = c1->child[3];
			tmp2[0] = c1->child[7];
		}

		if(c2->isLeaf){
			tmp1[1] = c2;
			tmp2[1] = c2;
		}
		else{
			tmp1[1] = c2->child[2];
			tmp2[1] = c2->child[6];
		}

		if(c3->isLeaf){
			tmp1[2] = c3;
			tmp2[2] = c3;
		}
		else{
			tmp1[2] = c3->child[1];
			tmp2[2] = c3->child[5];
		}

		if(c4->isLeaf){
			tmp1[3] = c4;
			tmp2[3] = c4;
		}
		else{
			tmp1[3] = c4->child[0];
			tmp2[3] = c4->child[4];
		}
	}
	edgeProc(tmp1[0], tmp1[1], tmp1[2], tmp1[3], edge);
	edgeProc(tmp2[0], tmp2[1], tmp2[2], tmp2[3], edge);
}

void OctTreeP::polygonize(float (*vertex)[3], int **face, double (*Q)[10])
{
	this->face = face;
	this->vertex = vertex;
	this->Q = Q;
	setVerticeAtCellCenter();
	counter = 0;
	cellProc(root);
}

void OctTreeP::bisection(float p[], float p1[], float p2[], float f1, float f2)
{
	float p3[3];
	float f3;

	//regula falsa
	if(f1 > 0.0){
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

		bool flag;
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

	//bisection
	for(int j=0; j<5; j++){
		if(PolygonalMesh::DIST(p1, p2) < 0.0001){
			p[0] = p3[0];
			p[1] = p3[1];
			p[2] = p3[2];
		}

		p3[0] = 0.5f*(p1[0] + p2[0]);
		p3[1] = 0.5f*(p1[1] + p2[1]);
		p3[2] = 0.5f*(p1[2] + p2[2]);
		f3 = func->value(p3[0], p3[1], p3[2]);
		if(fabs(f3) < 0.000001){
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
	p[0] = p3[0];
	p[1] = p3[1];
	p[2] = p3[2];
}

void OctTreeP::countVertexAndFace(int &vertex_N, int &face_N)
{
	counter = 0;
	counter2 = 0;
	cellProc2(root);
	vertex_N = counter;
	face_N = counter2;
}

void OctTreeP::cellProc2(Cell *c)
{
	if(c->isLeaf)
		return;
	
	for(int i=0; i<8; i++)
		cellProc2(c->child[i]);

	faceProc2(c->child[EDGE[0][0]], c->child[EDGE[0][1]], RIGHT);
	faceProc2(c->child[EDGE[1][0]], c->child[EDGE[1][1]], RIGHT);
	faceProc2(c->child[EDGE[2][0]], c->child[EDGE[2][1]], FRONT);
	faceProc2(c->child[EDGE[3][0]], c->child[EDGE[3][1]], FRONT);
	faceProc2(c->child[EDGE[4][0]], c->child[EDGE[4][1]], RIGHT);
	faceProc2(c->child[EDGE[5][0]], c->child[EDGE[5][1]], RIGHT);
	faceProc2(c->child[EDGE[6][0]], c->child[EDGE[6][1]], FRONT);
	faceProc2(c->child[EDGE[7][0]], c->child[EDGE[7][1]], FRONT);
	faceProc2(c->child[EDGE[8][0]], c->child[EDGE[8][1]], UP);
	faceProc2(c->child[EDGE[9][0]], c->child[EDGE[9][1]], UP);
	faceProc2(c->child[EDGE[10][0]], c->child[EDGE[10][1]], UP);
	faceProc2(c->child[EDGE[11][0]], c->child[EDGE[11][1]], UP);

	edgeProc2(c->child[FACE[0][0]], c->child[FACE[0][1]], c->child[FACE[0][2]], c->child[FACE[0][3]], FR);
	edgeProc2(c->child[FACE[1][0]], c->child[FACE[1][1]], c->child[FACE[1][2]], c->child[FACE[1][3]], FR);
	edgeProc2(c->child[FACE[2][0]], c->child[FACE[2][1]], c->child[FACE[2][2]], c->child[FACE[2][3]], RU);
	edgeProc2(c->child[FACE[3][0]], c->child[FACE[3][1]], c->child[FACE[3][2]], c->child[FACE[3][3]], RU);
	edgeProc2(c->child[FACE[4][0]], c->child[FACE[4][1]], c->child[FACE[4][2]], c->child[FACE[4][3]], FU);
	edgeProc2(c->child[FACE[5][0]], c->child[FACE[5][1]], c->child[FACE[5][2]], c->child[FACE[5][3]], FU);
}

void OctTreeP::faceProc2(Cell *c1, Cell *c2, int face)
{
	if(c1->isLeaf && c2->isLeaf)
		return;

	Cell *tmp1[4];
	if(c1->isLeaf){
		for(int i=0; i<4; i++)
			tmp1[i] = c1;
	}
	else{
		for(int i=0; i<4; i++)
			tmp1[i] = c1->child[FACE[face][i]];
	}

	Cell *tmp2[4];
	if(c2->isLeaf){
		for(int i=0; i<4; i++)
			tmp2[i] = c2;
	}
	else{
		for(int i=0; i<4; i++)
			tmp2[i] = c2->child[FACE[face-1][i]];
	}

	for(int i=0; i<4; i++)
		faceProc2(tmp1[i], tmp2[i], face);

	if(face == RIGHT){
		edgeProc2(tmp1[0], tmp2[0], tmp1[2], tmp2[2], RU);
		edgeProc2(tmp1[1], tmp2[1], tmp1[3], tmp2[3], RU);
		edgeProc2(tmp1[0], tmp2[0], tmp1[1], tmp2[1], FR);
		edgeProc2(tmp1[2], tmp2[2], tmp1[3], tmp2[3], FR);
	}
	else if(face == UP){
		edgeProc2(tmp1[0], tmp1[2], tmp2[0], tmp2[2], FU);
		//edgeProc(tmp1[0], tmp2[0], tmp1[2], tmp2[2], FU);
		//edgeProc(tmp1[1], tmp2[1], tmp1[3], tmp2[3], FU);
		edgeProc2(tmp1[1], tmp1[3], tmp2[1], tmp2[3], FU);
		edgeProc2(tmp1[0], tmp1[1], tmp2[0], tmp2[1], RU);
		edgeProc2(tmp1[2], tmp1[3], tmp2[2], tmp2[3], RU);
	}
	else{
		edgeProc2(tmp1[0], tmp2[0], tmp1[2], tmp2[2], FU);
		edgeProc2(tmp1[1], tmp2[1], tmp1[3], tmp2[3], FU);
		edgeProc2(tmp1[0], tmp1[1], tmp2[0], tmp2[1], FR);
		edgeProc2(tmp1[2], tmp1[3], tmp2[2], tmp2[3], FR);
	}
}

void OctTreeP::edgeProc2(Cell *c1, Cell *c2, Cell *c3, Cell *c4, int edge)
{
	if(c1->isLeaf && c2->isLeaf && c3->isLeaf && c4->isLeaf){
		int max = c1->level;
		int index = 1;
		if(max < c2->level){
			max = c2->level;
			index = 2;
		}
		if(max < c3->level){
			max = c3->level;
			index = 3;
		}
		if(max < c4->level){
			max = c4->level;
			index = 4;
		}

		float f1, f2;
		if(index == 1){
			if(edge == RU){
				f1 = c1->value[7];
				f2 = c1->value[5];
			}
			else if(edge == FU){
				f1 = c1->value[6];
				f2 = c1->value[7];
			}
			else{
				f1 = c1->value[3];
				f2 = c1->value[7];
			}
		}
		else if(index == 2){
			if(edge == RU){
				f1 = c2->value[6];
				f2 = c2->value[4];
			}
			else if(edge == FU){
				f1 = c2->value[4];
				f2 = c2->value[5];
			}
			else{
				f1 = c2->value[2];
				f2 = c2->value[6];
			}
		}
		else if(index == 3){
			if(edge == RU){
				f1 = c3->value[3];
				f2 = c3->value[1];
			}
			else if(edge == FU){
				f1 = c3->value[2];
				f2 = c3->value[3];
			}
			else{
				f1 = c3->value[1];
				f2 = c3->value[5];
			}
		}
		else{
			if(edge == RU){
				f1 = c4->value[2];
				f2 = c4->value[0];
			}
			else if(edge == FU){
				f1 = c4->value[0];
				f2 = c4->value[1];
			}
			else{
				f1 = c4->value[0];
				f2 = c4->value[4];
			}
		}
	
		if((f1 < 0 || f2 >= 0) && (f2 < 0 || f1 >= 0))
			return;

		counter2++;

		if(c1->vertex_id < 0){
			c1->vertex_id = counter;
			counter++;
		}
		if(c2->vertex_id < 0){
			c2->vertex_id = counter;
			counter++;
		}
		if(c3->vertex_id < 0){
			c3->vertex_id = counter;
			counter++;
		}
		if(c4->vertex_id < 0){
			c4->vertex_id = counter;
			counter++;
		}
		return;
	}

	Cell *tmp1[4], *tmp2[4];
	if(edge == FU){
		if(c1->isLeaf){
			tmp1[0] = c1;
			tmp2[0] = c1;
		}
		else{
			tmp1[0] = c1->child[6];
			tmp2[0] = c1->child[7];
		}

		if(c2->isLeaf){
			tmp1[1] = c2;
			tmp2[1] = c2;
		}
		else{
			tmp1[1] = c2->child[4];
			tmp2[1] = c2->child[5];
		}

		if(c3->isLeaf){
			tmp1[2] = c3;
			tmp2[2] = c3;
		}
		else{
			tmp1[2] = c3->child[2];
			tmp2[2] = c3->child[3];
		}

		if(c4->isLeaf){
			tmp1[3] = c4;
			tmp2[3] = c4;
		}
		else{
			tmp1[3] = c4->child[0];
			tmp2[3] = c4->child[1];
		}
	}
	else if(edge == RU){
		if(c1->isLeaf){
			tmp1[0] = c1;
			tmp2[0] = c1;
		}
		else{
			tmp1[0] = c1->child[5];
			tmp2[0] = c1->child[7];
		}

		if(c2->isLeaf){
			tmp1[1] = c2;
			tmp2[1] = c2;
		}
		else{
			tmp1[1] = c2->child[4];
			tmp2[1] = c2->child[6];
		}

		if(c3->isLeaf){
			tmp1[2] = c3;
			tmp2[2] = c3;
		}
		else{
			tmp1[2] = c3->child[1];
			tmp2[2] = c3->child[3];
		}

		if(c4->isLeaf){
			tmp1[3] = c4;
			tmp2[3] = c4;
		}
		else{
			tmp1[3] = c4->child[0];
			tmp2[3] = c4->child[2];
		}
	}
	else{
		if(c1->isLeaf){
			tmp1[0] = c1;
			tmp2[0] = c1;
		}
		else{
			tmp1[0] = c1->child[3];
			tmp2[0] = c1->child[7];
		}

		if(c2->isLeaf){
			tmp1[1] = c2;
			tmp2[1] = c2;
		}
		else{
			tmp1[1] = c2->child[2];
			tmp2[1] = c2->child[6];
		}

		if(c3->isLeaf){
			tmp1[2] = c3;
			tmp2[2] = c3;
		}
		else{
			tmp1[2] = c3->child[1];
			tmp2[2] = c3->child[5];
		}

		if(c4->isLeaf){
			tmp1[3] = c4;
			tmp2[3] = c4;
		}
		else{
			tmp1[3] = c4->child[0];
			tmp2[3] = c4->child[4];
		}
	}
	edgeProc2(tmp1[0], tmp1[1], tmp1[2], tmp1[3], edge);
	edgeProc2(tmp2[0], tmp2[1], tmp2[2], tmp2[3], edge);
}

void OctTreeP::simplify(double (*Q)[10], float tol)
{
	this->Q = Q;
	cellProcQ(root);
	if(tol != 0)
		simplifyCell(root, tol);
}

void OctTreeP::cellProcQ(Cell *c)
{
	if(c->isLeaf)
		return;
	
	for(int i=0; i<8; i++)
		cellProcQ(c->child[i]);

	faceProcQ(c->child[EDGE[0][0]], c->child[EDGE[0][1]], RIGHT);
	faceProcQ(c->child[EDGE[1][0]], c->child[EDGE[1][1]], RIGHT);
	faceProcQ(c->child[EDGE[2][0]], c->child[EDGE[2][1]], FRONT);
	faceProcQ(c->child[EDGE[3][0]], c->child[EDGE[3][1]], FRONT);
	faceProcQ(c->child[EDGE[4][0]], c->child[EDGE[4][1]], RIGHT);
	faceProcQ(c->child[EDGE[5][0]], c->child[EDGE[5][1]], RIGHT);
	faceProcQ(c->child[EDGE[6][0]], c->child[EDGE[6][1]], FRONT);
	faceProcQ(c->child[EDGE[7][0]], c->child[EDGE[7][1]], FRONT);
	faceProcQ(c->child[EDGE[8][0]], c->child[EDGE[8][1]], UP);
	faceProcQ(c->child[EDGE[9][0]], c->child[EDGE[9][1]], UP);
	faceProcQ(c->child[EDGE[10][0]], c->child[EDGE[10][1]], UP);
	faceProcQ(c->child[EDGE[11][0]], c->child[EDGE[11][1]], UP);

	edgeProcQ(c->child[FACE[0][0]], c->child[FACE[0][1]], c->child[FACE[0][2]], c->child[FACE[0][3]], FR);
	edgeProcQ(c->child[FACE[1][0]], c->child[FACE[1][1]], c->child[FACE[1][2]], c->child[FACE[1][3]], FR);
	edgeProcQ(c->child[FACE[2][0]], c->child[FACE[2][1]], c->child[FACE[2][2]], c->child[FACE[2][3]], RU);
	edgeProcQ(c->child[FACE[3][0]], c->child[FACE[3][1]], c->child[FACE[3][2]], c->child[FACE[3][3]], RU);
	edgeProcQ(c->child[FACE[4][0]], c->child[FACE[4][1]], c->child[FACE[4][2]], c->child[FACE[4][3]], FU);
	edgeProcQ(c->child[FACE[5][0]], c->child[FACE[5][1]], c->child[FACE[5][2]], c->child[FACE[5][3]], FU);
}

void OctTreeP::faceProcQ(Cell *c1, Cell *c2, int face)
{
	if(c1->isLeaf && c2->isLeaf)
		return;

	Cell *tmp1[4];
	if(c1->isLeaf){
		for(int i=0; i<4; i++)
			tmp1[i] = c1;
	}
	else{
		for(int i=0; i<4; i++)
			tmp1[i] = c1->child[FACE[face][i]];
	}

	Cell *tmp2[4];
	if(c2->isLeaf){
		for(int i=0; i<4; i++)
			tmp2[i] = c2;
	}
	else{
		for(int i=0; i<4; i++)
			tmp2[i] = c2->child[FACE[face-1][i]];
	}

	for(int i=0; i<4; i++)
		faceProcQ(tmp1[i], tmp2[i], face);

	if(face == RIGHT){
		edgeProcQ(tmp1[0], tmp2[0], tmp1[2], tmp2[2], RU);
		edgeProcQ(tmp1[1], tmp2[1], tmp1[3], tmp2[3], RU);
		edgeProcQ(tmp1[0], tmp2[0], tmp1[1], tmp2[1], FR);
		edgeProcQ(tmp1[2], tmp2[2], tmp1[3], tmp2[3], FR);
	}
	else if(face == UP){
		edgeProcQ(tmp1[0], tmp1[2], tmp2[0], tmp2[2], FU);
		//edgeProc(tmp1[0], tmp2[0], tmp1[2], tmp2[2], FU);
		//edgeProc(tmp1[1], tmp2[1], tmp1[3], tmp2[3], FU);
		edgeProcQ(tmp1[1], tmp1[3], tmp2[1], tmp2[3], FU);
		edgeProcQ(tmp1[0], tmp1[1], tmp2[0], tmp2[1], RU);
		edgeProcQ(tmp1[2], tmp1[3], tmp2[2], tmp2[3], RU);
	}
	else{
		edgeProcQ(tmp1[0], tmp2[0], tmp1[2], tmp2[2], FU);
		edgeProcQ(tmp1[1], tmp2[1], tmp1[3], tmp2[3], FU);
		edgeProcQ(tmp1[0], tmp1[1], tmp2[0], tmp2[1], FR);
		edgeProcQ(tmp1[2], tmp1[3], tmp2[2], tmp2[3], FR);
	}
}

void OctTreeP::edgeProcQ(Cell *c1, Cell *c2, Cell *c3, Cell *c4, int edge)
{
	if(c1->isLeaf && c2->isLeaf && c3->isLeaf && c4->isLeaf){
		int i1 = c1->vertex_id;
		int i2 = c2->vertex_id;
		int i3 = c3->vertex_id;
		int i4 = c4->vertex_id;
		if(i1 == -1 || i2 == -1 || i3 == -1 || i4== -1)
			return;
		int max = c1->level;
		int index = 1;
		if(max < c2->level){
			max = c2->level;
			index = 2;
		}
		if(max < c3->level){
			max = c3->level;
			index = 3;
		}
		if(max < c4->level){
			max = c4->level;
			index = 4;
		}

		float p1[3], p2[3], f1, f2;
		if(index == 1){
			if(edge == RU){
				f1 = c1->value[7];
				f2 = c1->value[5];
				CONNER_P(p1, c1, 7);
				CONNER_P(p2, c1, 5);
			}
			else if(edge == FU){
				f1 = c1->value[6];
				f2 = c1->value[7];
				CONNER_P(p1, c1, 6);
				CONNER_P(p2, c1, 7);
			}
			else{
				f1 = c1->value[3];
				f2 = c1->value[7];
				CONNER_P(p1, c1, 3);
				CONNER_P(p2, c1, 7);
			}
		}
		else if(index == 2){
			if(edge == RU){
				f1 = c2->value[6];
				f2 = c2->value[4];
				CONNER_P(p1, c2, 6);
				CONNER_P(p2, c2, 4);
			}
			else if(edge == FU){
				f1 = c2->value[4];
				f2 = c2->value[5];
				CONNER_P(p1, c2, 4);
				CONNER_P(p2, c2, 5);
			}
			else{
				f1 = c2->value[2];
				f2 = c2->value[6];
				CONNER_P(p1, c2, 2);
				CONNER_P(p2, c2, 6);
			}
		}
		else if(index == 3){
			if(edge == RU){
				f1 = c3->value[3];
				f2 = c3->value[1];
				CONNER_P(p1, c3, 3);
				CONNER_P(p2, c3, 1);
			}
			else if(edge == FU){
				f1 = c3->value[2];
				f2 = c3->value[3];
				CONNER_P(p1, c3, 2);
				CONNER_P(p2, c3, 3);
			}
			else{
				f1 = c3->value[1];
				f2 = c3->value[5];
				CONNER_P(p1, c3, 1);
				CONNER_P(p2, c3, 5);
			}
		}
		else{
			if(edge == RU){
				f1 = c4->value[2];
				f2 = c4->value[0];
				CONNER_P(p1, c4, 2);
				CONNER_P(p2, c4, 0);
			}
			else if(edge == FU){
				f1 = c4->value[0];
				f2 = c4->value[1];
				CONNER_P(p1, c4, 0);
				CONNER_P(p2, c4, 1);
			}
			else{
				f1 = c4->value[0];
				f2 = c4->value[4];
				CONNER_P(p1, c4, 0);
				CONNER_P(p2, c4, 4);
			}
		}
		
		if((f1 < 0 || f2 >= 0) && (f2 < 0 || f1 >= 0))
			return;

		//bisection
		float p[3];
		bisection(p, p1, p2, f1, f2);
		float g[3];
		func->gradient(g, p[0], p[1], p[2]);
		double len = PolygonalMesh::LENGTH(g);
		double nor[3];
		if((float)len != 0){
			nor[0] = g[0]/len;
			nor[1] = g[1]/len;
			nor[2] = g[2]/len;
		}
		else
			nor[0] = nor[1] = nor[2] = 0;

		float v[3];
		float center[3];
		//CELL_P(center, c1);
		//PolygonalMesh::VEC(v, center, p);
		double d = -PolygonalMesh::DOT(p, nor);
		double Q_tmp[10];
		MATRIX(Q_tmp, nor, d);
		MAT_SUM(Q[c1->vertex_id], Q_tmp);

		if(c2->vertex_id != c1->vertex_id){
			CELL_P(center, c2);
			//PolygonalMesh::VEC(v, center, p);
			//d = -PolygonalMesh::DOT(p, nor);
			MATRIX(Q_tmp, nor, d);
			MAT_SUM(Q[c2->vertex_id], Q_tmp);
		}

		if(c3->vertex_id != c1->vertex_id 
			&& c3->vertex_id != c2->vertex_id){
			//CELL_P(center, c3);
			//PolygonalMesh::VEC(v, center, p);
			d = -PolygonalMesh::DOT(p, nor);
			MATRIX(Q_tmp, nor, d);
			MAT_SUM(Q[c3->vertex_id], Q_tmp);
		}

		if(c4->vertex_id != c1->vertex_id 
			&& c4->vertex_id != c2->vertex_id
			&& c4->vertex_id != c3->vertex_id){
			//CELL_P(center, c4);
			//PolygonalMesh::VEC(v, center, p);
			d = -PolygonalMesh::DOT(p, nor);
			MATRIX(Q_tmp, nor, d);
			MAT_SUM(Q[c4->vertex_id], Q_tmp);
		}
		return;
	}

	Cell *tmp1[4], *tmp2[4];
	if(edge == FU){
		if(c1->isLeaf){
			tmp1[0] = c1;
			tmp2[0] = c1;
		}
		else{
			tmp1[0] = c1->child[6];
			tmp2[0] = c1->child[7];
		}

		if(c2->isLeaf){
			tmp1[1] = c2;
			tmp2[1] = c2;
		}
		else{
			tmp1[1] = c2->child[4];
			tmp2[1] = c2->child[5];
		}

		if(c3->isLeaf){
			tmp1[2] = c3;
			tmp2[2] = c3;
		}
		else{
			tmp1[2] = c3->child[2];
			tmp2[2] = c3->child[3];
		}

		if(c4->isLeaf){
			tmp1[3] = c4;
			tmp2[3] = c4;
		}
		else{
			tmp1[3] = c4->child[0];
			tmp2[3] = c4->child[1];
		}
	}
	else if(edge == RU){
		if(c1->isLeaf){
			tmp1[0] = c1;
			tmp2[0] = c1;
		}
		else{
			tmp1[0] = c1->child[5];
			tmp2[0] = c1->child[7];
		}

		if(c2->isLeaf){
			tmp1[1] = c2;
			tmp2[1] = c2;
		}
		else{
			tmp1[1] = c2->child[4];
			tmp2[1] = c2->child[6];
		}

		if(c3->isLeaf){
			tmp1[2] = c3;
			tmp2[2] = c3;
		}
		else{
			tmp1[2] = c3->child[1];
			tmp2[2] = c3->child[3];
		}

		if(c4->isLeaf){
			tmp1[3] = c4;
			tmp2[3] = c4;
		}
		else{
			tmp1[3] = c4->child[0];
			tmp2[3] = c4->child[2];
		}
	}
	else{
		if(c1->isLeaf){
			tmp1[0] = c1;
			tmp2[0] = c1;
		}
		else{
			tmp1[0] = c1->child[3];
			tmp2[0] = c1->child[7];
		}

		if(c2->isLeaf){
			tmp1[1] = c2;
			tmp2[1] = c2;
		}
		else{
			tmp1[1] = c2->child[2];
			tmp2[1] = c2->child[6];
		}

		if(c3->isLeaf){
			tmp1[2] = c3;
			tmp2[2] = c3;
		}
		else{
			tmp1[2] = c3->child[1];
			tmp2[2] = c3->child[5];
		}

		if(c4->isLeaf){
			tmp1[3] = c4;
			tmp2[3] = c4;
		}
		else{
			tmp1[3] = c4->child[0];
			tmp2[3] = c4->child[4];
		}
	}
	edgeProcQ(tmp1[0], tmp1[1], tmp1[2], tmp1[3], edge);
	edgeProcQ(tmp2[0], tmp2[1], tmp2[2], tmp2[3], edge);
}

void OctTreeP::simplifyCell(Cell *c, float tol)
{
	if(c->isLeaf)
		return;

	for(int i=0; i<8; i++)
		simplifyCell(c->child[i], tol);		
	
	for(int i=0; i<8; i++){
		if(!c->child[i]->isLeaf)
			return;
	}

	int index = -1;
	double Q[10];
	MAT_INIT(Q);
	for(int i=0; i<8; i++){
		int id = c->child[i]->vertex_id;
		if(id < 0)
			continue;
		if(index < 0)
			index = id;
		MAT_SUM(Q, this->Q[id]);
	}
	if(index < 0)
		return;
	float x[3];
	computeOptP(x, Q, c);
	float err = Q_ERR(Q, x);
	if(err < tol*tol){
		c->vertex_id = index;
		MAT_COPY(this->Q[index], Q);
		for(int i=0; i<8; i++){
			delete[] c->child[i];
		}
		c->isLeaf = true;
	}
}

void OctTreeP::computeOptP(float P[], double Q[], Cell* c)
{
  vnl_matrix< float > A( 3, 3 ); 
  vnl_vector<float> b(3), x(3);

	A[1][1] = Q[0];
	A[2][1] = A[1][2] = Q[1];
	A[3][1] = A[1][3] = Q[2];
	A[2][2] = Q[3];
	A[2][3] = A[3][2] = Q[4];
	A[3][3] = Q[5];

	float cen[3];
	CELL_P(cen, c);
	float tmp[3];
	MAT_BY_VEC(tmp, Q, cen);
		
	b[1] = -Q[6] - tmp[0];
	b[2] = -Q[7] - tmp[1];
	b[3] = -Q[8] - tmp[2];

	//b[1] = -Q[6];
	//b[2] = -Q[7];
	//b[3] = -Q[8];

  vnl_svd<float> svd( A );
  svd.zero_out_absolute( 0.0000001 );

  x = svd.solve( b );

	P[0] = cen[0] + x[1];
	P[1] = cen[1] + x[2];
	P[2] = cen[2] + x[3];
}

void OctTreeP::polygonize2(float (*vertex)[3], int **face)
{
	this->face = face;
	this->vertex = vertex;
	setVerticeAtCellCenter();
	counter = 0;
	cellProcS(root);
}

void OctTreeP::cellProcS(Cell *c)
{
	if(c->isLeaf)
		return;
	
	for(int i=0; i<8; i++)
		cellProcS(c->child[i]);

	faceProcS(c->child[EDGE[0][0]], c->child[EDGE[0][1]], RIGHT);
	faceProcS(c->child[EDGE[1][0]], c->child[EDGE[1][1]], RIGHT);
	faceProcS(c->child[EDGE[2][0]], c->child[EDGE[2][1]], FRONT);
	faceProcS(c->child[EDGE[3][0]], c->child[EDGE[3][1]], FRONT);
	faceProcS(c->child[EDGE[4][0]], c->child[EDGE[4][1]], RIGHT);
	faceProcS(c->child[EDGE[5][0]], c->child[EDGE[5][1]], RIGHT);
	faceProcS(c->child[EDGE[6][0]], c->child[EDGE[6][1]], FRONT);
	faceProcS(c->child[EDGE[7][0]], c->child[EDGE[7][1]], FRONT);
	faceProcS(c->child[EDGE[8][0]], c->child[EDGE[8][1]], UP);
	faceProcS(c->child[EDGE[9][0]], c->child[EDGE[9][1]], UP);
	faceProcS(c->child[EDGE[10][0]], c->child[EDGE[10][1]], UP);
	faceProcS(c->child[EDGE[11][0]], c->child[EDGE[11][1]], UP);

	edgeProcS(c->child[FACE[0][0]], c->child[FACE[0][1]], c->child[FACE[0][2]], c->child[FACE[0][3]], FR);
	edgeProcS(c->child[FACE[1][0]], c->child[FACE[1][1]], c->child[FACE[1][2]], c->child[FACE[1][3]], FR);
	edgeProcS(c->child[FACE[2][0]], c->child[FACE[2][1]], c->child[FACE[2][2]], c->child[FACE[2][3]], RU);
	edgeProcS(c->child[FACE[3][0]], c->child[FACE[3][1]], c->child[FACE[3][2]], c->child[FACE[3][3]], RU);
	edgeProcS(c->child[FACE[4][0]], c->child[FACE[4][1]], c->child[FACE[4][2]], c->child[FACE[4][3]], FU);
	edgeProcS(c->child[FACE[5][0]], c->child[FACE[5][1]], c->child[FACE[5][2]], c->child[FACE[5][3]], FU);
}

void OctTreeP::faceProcS(Cell *c1, Cell *c2, int face)
{
	if(c1->isLeaf && c2->isLeaf)
		return;

	Cell *tmp1[4];
	if(c1->isLeaf){
		for(int i=0; i<4; i++)
			tmp1[i] = c1;
	}
	else{
		for(int i=0; i<4; i++)
			tmp1[i] = c1->child[FACE[face][i]];
	}

	Cell *tmp2[4];
	if(c2->isLeaf){
		for(int i=0; i<4; i++)
			tmp2[i] = c2;
	}
	else{
		for(int i=0; i<4; i++)
			tmp2[i] = c2->child[FACE[face-1][i]];
	}

	for(int i=0; i<4; i++)
		faceProcS(tmp1[i], tmp2[i], face);

	if(face == RIGHT){
		edgeProcS(tmp1[0], tmp2[0], tmp1[2], tmp2[2], RU);
		edgeProcS(tmp1[1], tmp2[1], tmp1[3], tmp2[3], RU);
		edgeProcS(tmp1[0], tmp2[0], tmp1[1], tmp2[1], FR);
		edgeProcS(tmp1[2], tmp2[2], tmp1[3], tmp2[3], FR);
	}
	else if(face == UP){
		edgeProcS(tmp1[0], tmp1[2], tmp2[0], tmp2[2], FU);
		//edgeProc(tmp1[0], tmp2[0], tmp1[2], tmp2[2], FU);
		//edgeProc(tmp1[1], tmp2[1], tmp1[3], tmp2[3], FU);
		edgeProcS(tmp1[1], tmp1[3], tmp2[1], tmp2[3], FU);
		edgeProcS(tmp1[0], tmp1[1], tmp2[0], tmp2[1], RU);
		edgeProcS(tmp1[2], tmp1[3], tmp2[2], tmp2[3], RU);
	}
	else{
		edgeProcS(tmp1[0], tmp2[0], tmp1[2], tmp2[2], FU);
		edgeProcS(tmp1[1], tmp2[1], tmp1[3], tmp2[3], FU);
		edgeProcS(tmp1[0], tmp1[1], tmp2[0], tmp2[1], FR);
		edgeProcS(tmp1[2], tmp1[3], tmp2[2], tmp2[3], FR);
	}
}

void OctTreeP::edgeProcS(Cell *c1, Cell *c2, Cell *c3, Cell *c4, int edge)
{
	if(c1->isLeaf && c2->isLeaf && c3->isLeaf && c4->isLeaf){
		int i1 = c1->vertex_id;
		int i2 = c2->vertex_id;
		int i3 = c3->vertex_id;
		int i4 = c4->vertex_id;
		if(i1 == -1 || i2 == -1 || i3 == -1 || i4== -1)
			return;
		int max = c1->level;
		int index = 1;
		if(max < c2->level){
			max = c2->level;
			index = 2;
		}
		if(max < c3->level){
			max = c3->level;
			index = 3;
		}
		if(max < c4->level){
			max = c4->level;
			index = 4;
		}

		float f1, f2;
		if(index == 1){
			if(edge == RU){
				f1 = c1->value[7];
				f2 = c1->value[5];
			}
			else if(edge == FU){
				f1 = c1->value[6];
				f2 = c1->value[7];
			}
			else{
				f1 = c1->value[3];
				f2 = c1->value[7];
			}
		}
		else if(index == 2){
			if(edge == RU){
				f1 = c2->value[6];
				f2 = c2->value[4];
			}
			else if(edge == FU){
				f1 = c2->value[4];
				f2 = c2->value[5];
			}
			else{
				f1 = c2->value[2];
				f2 = c2->value[6];
			}
		}
		else if(index == 3){
			if(edge == RU){
				f1 = c3->value[3];
				f2 = c3->value[1];
			}
			else if(edge == FU){
				f1 = c3->value[2];
				f2 = c3->value[3];
			}
			else{
				f1 = c3->value[1];
				f2 = c3->value[5];
			}
		}
		else{
			if(edge == RU){
				f1 = c4->value[2];
				f2 = c4->value[0];
			}
			else if(edge == FU){
				f1 = c4->value[0];
				f2 = c4->value[1];
			}
			else{
				f1 = c4->value[0];
				f2 = c4->value[4];
			}
		}
		
		if((f1 < 0 || f2 >= 0) && (f2 < 0 || f1 >= 0))
			return;

		if(f1 > 0){
			face[counter][0] = c1->vertex_id;
			face[counter][1] = c2->vertex_id;
			face[counter][2] = c4->vertex_id;
			face[counter][3] = c3->vertex_id;
		}
		else{
			face[counter][0] = c3->vertex_id;
			face[counter][1] = c4->vertex_id;
			face[counter][2] = c2->vertex_id;
			face[counter][3] = c1->vertex_id;
		}
		counter++;

		return;
	}

	Cell *tmp1[4], *tmp2[4];
	if(edge == FU){
		if(c1->isLeaf){
			tmp1[0] = c1;
			tmp2[0] = c1;
		}
		else{
			tmp1[0] = c1->child[6];
			tmp2[0] = c1->child[7];
		}

		if(c2->isLeaf){
			tmp1[1] = c2;
			tmp2[1] = c2;
		}
		else{
			tmp1[1] = c2->child[4];
			tmp2[1] = c2->child[5];
		}

		if(c3->isLeaf){
			tmp1[2] = c3;
			tmp2[2] = c3;
		}
		else{
			tmp1[2] = c3->child[2];
			tmp2[2] = c3->child[3];
		}

		if(c4->isLeaf){
			tmp1[3] = c4;
			tmp2[3] = c4;
		}
		else{
			tmp1[3] = c4->child[0];
			tmp2[3] = c4->child[1];
		}
	}
	else if(edge == RU){
		if(c1->isLeaf){
			tmp1[0] = c1;
			tmp2[0] = c1;
		}
		else{
			tmp1[0] = c1->child[5];
			tmp2[0] = c1->child[7];
		}

		if(c2->isLeaf){
			tmp1[1] = c2;
			tmp2[1] = c2;
		}
		else{
			tmp1[1] = c2->child[4];
			tmp2[1] = c2->child[6];
		}

		if(c3->isLeaf){
			tmp1[2] = c3;
			tmp2[2] = c3;
		}
		else{
			tmp1[2] = c3->child[1];
			tmp2[2] = c3->child[3];
		}

		if(c4->isLeaf){
			tmp1[3] = c4;
			tmp2[3] = c4;
		}
		else{
			tmp1[3] = c4->child[0];
			tmp2[3] = c4->child[2];
		}
	}
	else{
		if(c1->isLeaf){
			tmp1[0] = c1;
			tmp2[0] = c1;
		}
		else{
			tmp1[0] = c1->child[3];
			tmp2[0] = c1->child[7];
		}

		if(c2->isLeaf){
			tmp1[1] = c2;
			tmp2[1] = c2;
		}
		else{
			tmp1[1] = c2->child[2];
			tmp2[1] = c2->child[6];
		}

		if(c3->isLeaf){
			tmp1[2] = c3;
			tmp2[2] = c3;
		}
		else{
			tmp1[2] = c3->child[1];
			tmp2[2] = c3->child[5];
		}

		if(c4->isLeaf){
			tmp1[3] = c4;
			tmp2[3] = c4;
		}
		else{
			tmp1[3] = c4->child[0];
			tmp2[3] = c4->child[4];
		}
	}
	edgeProcS(tmp1[0], tmp1[1], tmp1[2], tmp1[3], edge);
	edgeProcS(tmp2[0], tmp2[1], tmp2[2], tmp2[3], edge);
}
