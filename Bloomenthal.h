// Bloomenthal.h: Bloomenthal's polygonizer 
//
//////////////////////////////////////////////////////////////////////

/*
#if !defined(AFX_BLOOMENTHAL_H__AAB677A8_023A_4833_83BC_FB97144ECC3B__INCLUDED_)
#define AFX_BLOOMENTHAL_H__AAB677A8_023A_4833_83BC_FB97144ECC3B__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
  */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include <sys/types.h>

#include "ImplicitFunction.h"
#include "PolygonalMesh.h"

/*
 * C code from the article
 * "An Implicit Surface Polygonizer"
 * by Jules Bloomenthal, jbloom@beauty.gmu.edu
 * in "Graphics Gems IV", Academic Press, 1994
 */

/* implicit.c
 *     an implicit surface polygonizer, translated from Mesa
 *     applications should call polygonize()
 *
 * to compile a test program for ASCII output:
 *     cc implicit.c -o implicit -lm
 *
 * to compile a test program for display on an SGI workstation:
 *     cc -DSGIGFX implicit.c -o implicit -lgl_s -lm
 *
 * Authored by Jules Bloomenthal, Xerox PARC.
 * Copyright (c) Xerox Corporation, 1991.  All rights reserved.
 * Permission is granted to reproduce, use and distribute this code for
 * any and all purposes, provided that this notice appears in all copies. */

namespace BLOOMENTHAL{

#define TET	0  /* use tetrahedral decomposition */
#define NOTET	1  /* no tetrahedral decomposition  */

#define RES	30 /* # converge iterations    */

#define L	0  /* left direction:	-x, -i */
#define R	1  /* right direction:	+x, +i */
#define B	2  /* bottom direction: -y, -j */
#define T	3  /* top direction:	+y, +j */
#define N	4  /* near direction:	-z, -k */
#define F	5  /* far direction:	+z, +k */
#define LBN	0  /* left bottom near corner  */
#define LBF	1  /* left bottom far corner   */
#define LTN	2  /* left top near corner     */
#define LTF	3  /* left top far corner      */
#define RBN	4  /* right bottom near corner */
#define RBF	5  /* right bottom far corner  */
#define RTN	6  /* right top near corner    */
#define RTF	7  /* right top far corner     */

/* the LBN corner of cube (i, j, k), corresponds with location
 * (start.x+(i-.5)*size, start.y+(j-.5)*size, start.z+(k-.5)*size) */

#define RAND()	    ((rand()&32767)/32767.)    /* random number between 0 and 1 */
#define HASHBIT	    (5)
#define HASHSIZE    (size_t)(1<<(3*HASHBIT))   /* hash table size (32768) */
#define MASK	    ((1<<HASHBIT)-1)
#define HASH(i,j,k) ((((((i)&MASK)<<HASHBIT)|((j)&MASK))<<HASHBIT)|((k)&MASK))
#define BIT(i, bit) (((i)>>(bit))&1)
#define FLIP(i,bit) ((i)^1<<(bit)) /* flip the given bit of i */

#define LB	0  /* left bottom edge	*/
#define LT	1  /* left top edge	*/
#define LN	2  /* left near edge	*/
#define LF	3  /* left far edge	*/
#define RB	4  /* right bottom edge */
#define RT	5  /* right top edge	*/
#define RN	6  /* right near edge	*/
#define RF	7  /* right far edge	*/
#define BN	8  /* bottom near edge	*/
#define BF	9  /* bottom far edge	*/
#define TN	10 /* top near edge	*/
#define TF	11 /* top far edge	*/

//void *calloc();
//char *mycalloc();

#define FILEOUT false
#define TRIFILE "tri.tmp"
#define VERFILE "ver.tmp"
#define MESHFILE "mesh.ply2"

class Bloomenthal  
{
	public:
typedef struct point {		   /* a three-dimensional point */
    double x, y, z;		   /* its coordinates */
} POINT1;

typedef struct test {		   /* test the function for a signed value */
    POINT1 p;			   /* location of test */
    double value;		   /* function value at p */
    int ok;			   /* if value is of correct sign */
} TEST;

typedef struct vertex {		   /* surface vertex */
    POINT1 position;//, normal;	   /* position and surface normal */
} VERTEX;

typedef struct vertices {	   /* list of vertices in polygonization */
    int count, max;		   /* # vertices, max # allowed */
    VERTEX *ptr;		   /* dynamically allocated */
} VERTICES;

typedef struct corner {		   /* corner of a cube */
    int i, j, k;		   /* (i, j, k) is index within lattice */
    double x, y, z, value;	   /* location and function value */
} CORNER;

typedef struct cube {		   /* partitioning cell (cube) */
    int i, j, k;		   /* lattice location of cube */
    CORNER *corners[8];		   /* eight corners */
} CUBE;

typedef struct cubes {		   /* linked list of cubes acting as stack */
    CUBE cube;			   /* a single cube */
    struct cubes *next;		   /* remaining elements */
} CUBES;

typedef struct centerlist {	   /* list of cube locations */
    int i, j, k;		   /* cube location */
    struct centerlist *next;	   /* remaining elements */
} CENTERLIST;

typedef struct cornerlist {	   /* list of corners */
    int i, j, k;		   /* corner id */
    double value;		   /* corner value */
    struct cornerlist *next;	   /* remaining elements */
} CORNERLIST;

typedef struct edgelist {	   /* list of edges */
    int i1, j1, k1, i2, j2, k2;	   /* edge corner ids */
    int vid;			   /* vertex id */
    struct edgelist *next;	   /* remaining elements */
} EDGELIST;

typedef struct intlist {	   /* list of integers */
    int i;			   /* an integer */
    struct intlist *next;	   /* remaining elements */
} INTLIST;

typedef struct trilist {
	int i1, i2, i3;
	struct trilist *next;
} TRILIST;

typedef struct intlists {	   /* list of list of integers */
    INTLIST *list;		   /* a list of integers */
    struct intlists *next;	   /* remaining elements */
} INTLISTS;

typedef struct process {	   /* parameters, function, storage */
    //double (*function)(double x, double y, double z);	   /* implicit surface function */
    //int (*triproc)(int i1, int i2, int i3, VERTICES vertices);		   /* triangle output function */
    double size, delta;		   /* cube size, normal delta */
    //int bounds;			   /* cube range within lattice */
	float boundsX1, boundsX2;
	float boundsY1, boundsY2;
	float boundsZ1, boundsZ2;
	POINT1 start;		   /* start point on surface */
    CUBES *cubes;		   /* active cubes */
    VERTICES vertices;		   /* surface vertices */
    CENTERLIST **centers;	   /* cube center hash table */
    CORNERLIST **corners;	   /* corner value hash table */
    EDGELIST **edges;		   /* edge and vertex id hash table */
} PROCESS;

public:
	ImplicitFunction *func;
	//INTLISTS *cubetable[256];

	int gntris;	     /* global needed by application */
	VERTICES gvertices;  /* global needed by application */
	TRILIST *tris;

	FILE* tri_file;
	FILE* ver_file;

public:
	double function(double x, double y, double z);
	void converge (POINT1 *p1, POINT1 *p2, double v, POINT1 *p);
	void vnormal (POINT1 *point, PROCESS *p, POINT1 *v);
	void addtovertices (VERTICES *vertices, VERTEX v);
	int vertid (CORNER *c1, CORNER *c2, PROCESS *p);
	int getedge (EDGELIST *table[], int i1, int j1, int k1, int i2, int j2, int k2);
	void setedge (EDGELIST *table[], int i1, int j1, int k1, int i2, int j2, int k2, int vid);
	int setcenter(CENTERLIST *table[], int i, int j, int k);
	char * mycalloc (int nitems, int nbytes);
	void makecubetable ();
	int otherface (int edge, int face);
	int nextcwedge (int edge, int face);
	int docube  (CUBE* cube, PROCESS* p);
	int dotet (CUBE* cube, int c1, int c2, int c3, int c4, PROCESS* p);
	TEST find1(int sign, PROCESS* p, double x, double y, double z);
	CORNER * setcorner1(PROCESS* p, int i, int j, int k);
	void testface (int i, int j, int k, CUBE* old, int face, int c1, int c2, int c3, int c4, PROCESS* p);
	char* polygonize (PolygonalMesh *mesh, double size, float bounds[6], double x, double y, double z, int mode);
	int triproc (int i1, int i2, int i3, VERTICES vertices);
	double blob (double x, double y, double z);
	double sphere (double x, double y, double z);
	double torus (double x, double y, double z);
	Bloomenthal();
	virtual ~Bloomenthal();
};
}

//#endif // !defined(AFX_BLOOMENTHAL_H__AAB677A8_023A_4833_83BC_FB97144ECC3B__INCLUDED_)
