// Bloomenthal.cpp: Bloomenthal 
//
//////////////////////////////////////////////////////////////////////

//#include "stdafx.h"
//#include "RBF3D.h"
#include "Bloomenthal.h"

/*
#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif
*/

using namespace BLOOMENTHAL;

//static int INTLISTS *cubetable[256];
static struct Bloomenthal::intlists *cubetable[256];

/*			edge: LB, LT, LN, LF, RB, RT, RN, RF, BN, BF, TN, TF */
static int corner1[12]	   = {LBN,LTN,LBN,LBF,RBN,RTN,RBN,RBF,LBN,LBF,LTN,LTF};
static int corner2[12]	   = {LBF,LTF,LTN,LTF,RBF,RTF,RTN,RTF,RBN,RBF,RTN,RTF};
static int leftface[12]	   = {B,  L,  L,  F,  R,  T,  N,  R,  N,  B,  T,  F};
			     /* face on left when going corner1 to corner2 */
static int rightface[12]   = {L,  T,  N,  L,  B,  R,  R,  F,  B,  F,  N,  T};
			     /* face on right when going corner1 to corner2 */

//////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////

Bloomenthal::Bloomenthal()
{

}

Bloomenthal::~Bloomenthal()
{

}


double Bloomenthal::torus(double x, double y, double z)
{
	double x2 = x*x, y2 = y*y, z2 = z*z;
    double a = x2+y2+z2+(0.5*0.5)-(0.1*0.1);
    return a*a-4.0*(0.5*0.5)*(y2+z2);
}

double Bloomenthal::sphere(double x, double y, double z)
{
	double rsq = x*x+y*y+z*z;
    return 1.0/(rsq < 0.00001? 0.00001 : rsq);
}


double Bloomenthal::blob(double x, double y, double z)
{
	return 4.0-sphere(x+1.0,y,z)-sphere(x,y+1.0,z)-sphere(x,y,z+1.0);
}


int Bloomenthal::triproc(int i1, int i2, int i3, VERTICES vertices)
{
	gvertices = vertices;
	gntris++;

	if(FILEOUT){
		fprintf(tri_file, "%d %d %d\n", i1, i2, i3);
	}
	else{
		TRILIST* n_tris = new TRILIST;
		n_tris->next = tris;
		n_tris->i1 = i1;
		n_tris->i2 = i2;
		n_tris->i3 = i3;
		tris = n_tris;
	}
    return 1;
}	

/**** An Implicit Surface Polygonizer ****/


/* polygonize: polygonize the implicit surface function
 *   arguments are:
 *	 double function (x, y, z)
 *		 double x, y, z (an arbitrary 3D point)
 *	     the implicit surface function
 *	     return negative for inside, positive for outside
 *	 double size
 *	     width of the partitioning cube
 *	 int bounds
 *	     max. range of cubes (+/- on the three axes) from first cube
 *	 double x, y, z
 *	     coordinates of a starting point on or near the surface
 *	     may be defaulted to 0., 0., 0.
 *	 int triproc (i1, i2, i3, vertices)
 *		 int i1, i2, i3 (indices into the vertex array)
 *		 VERTICES vertices (the vertex array, indexed from 0)
 *	     called for each triangle
 *	     the triangle coordinates are (for i = i1, i2, i3):
 *		 vertices.ptr[i].position.x, .y, and .z
 *	     vertices are ccw when viewed from the out (positive) side
 *		 in a left-handed coordinate system
 *	     vertex normals point outwards
 *	     return 1 to continue, 0 to abort
 *	 int mode
 *	     TET: decompose cube and polygonize six tetrahedra
 *	     NOTET: polygonize cube directly
 *   returns error or NULL
 */

int counter = 0;
char* Bloomenthal::polygonize(PolygonalMesh* mesh, double size, float bounds[6], double x, double y, double z, int mode)
{
	gntris = 0;
	if(FILEOUT){
		ver_file = fopen(VERFILE, "w");
		tri_file = fopen(TRIFILE, "w");
	}
	else{
		tris = new TRILIST;
		tris->i1 = -1;
		tris->i2 = -1;
		tris->i3 = -1;
		tris->next = NULL;
	}

	PROCESS p;
    int n, noabort;
    //CORNER *setcorner();
    TEST in, out;// find();

    //p.function = this->function;//function;
    //p.triproc = triproc;
    p.size = size;
    //p.bounds = bounds;
	p.boundsX1 = bounds[0];
	p.boundsX2 = bounds[1];
	p.boundsY1 = bounds[2];
	p.boundsY2 = bounds[3];
	p.boundsZ1 = bounds[4];
	p.boundsZ2 = bounds[5];
    p.delta = size/(double)(RES*RES);

    /* allocate hash tables and build cube polygon table: */
    p.centers = (CENTERLIST **) mycalloc(HASHSIZE,sizeof(CENTERLIST *));
    p.corners = (CORNERLIST **) mycalloc(HASHSIZE,sizeof(CORNERLIST *));
    p.edges =	(EDGELIST   **) mycalloc(2*HASHSIZE,sizeof(EDGELIST *));
    makecubetable();

    /* find point on surface, beginning search at (x, y, z): */
    srand(1);
    in = find1(1, &p, x, y, z);
    out = find1(0, &p, x, y, z);
    if (!in.ok || !out.ok) return (char*)"can't find starting point";
    converge(&in.p, &out.p, in.value, &p.start);

    /* push initial cube on stack: */
    p.cubes = (CUBES *) mycalloc(1, sizeof(CUBES)); /* list of 1 */
    p.cubes->cube.i = p.cubes->cube.j = p.cubes->cube.k = 0;
    p.cubes->next = NULL;

    /* set corners of initial cube: */
    for (n = 0; n < 8; n++)
	p.cubes->cube.corners[n] = setcorner1(&p, BIT(n,2), BIT(n,1), BIT(n,0));

    p.vertices.count = p.vertices.max = 0; /* no vertices yet */
    p.vertices.ptr = NULL;

    setcenter(p.centers, 0, 0, 0);

    while (p.cubes != NULL) { /* process active cubes till none left */
	CUBE c;
	CUBES *temp = p.cubes;
	c = p.cubes->cube;

	noabort = mode == TET?
	       /* either decompose into tetrahedra and polygonize: */
	       dotet(&c, LBN, LTN, RBN, LBF, &p) &&
	       dotet(&c, RTN, LTN, LBF, RBN, &p) &&
	       dotet(&c, RTN, LTN, LTF, LBF, &p) &&
	       dotet(&c, RTN, RBN, LBF, RBF, &p) &&
	       dotet(&c, RTN, LBF, LTF, RBF, &p) &&
	       dotet(&c, RTN, LTF, RTF, RBF, &p)
	       :
	       /* or polygonize the cube directly: */
	       docube(&c, &p);
	if (! noabort) return (char*)"aborted";

	/* pop current cube from stack */
	p.cubes = p.cubes->next;
	free((char *) temp);

	/* test six face directions, maybe add to stack: */
	testface(c.i-1, c.j, c.k, &c, L, LBN, LBF, LTN, LTF, &p);
	testface(c.i+1, c.j, c.k, &c, R, RBN, RBF, RTN, RTF, &p);
	testface(c.i, c.j-1, c.k, &c, B, LBN, LBF, RBN, RBF, &p);
	testface(c.i, c.j+1, c.k, &c, T, LTN, LTF, RTN, RTF, &p);
	testface(c.i, c.j, c.k-1, &c, N, LBN, LTN, RBN, RTN, &p);
	testface(c.i, c.j, c.k+1, &c, F, LBF, LTF, RBF, RTF, &p);

	for(n=0; n<8; n++)
		free((char *)c.corners[n]);
    }

	//free
	for(unsigned int i=0; i<HASHSIZE; i++){
		while(p.centers[i] != NULL){
			CENTERLIST* current = p.centers[i];
			p.centers[i] = current->next;
			free((char*)current);
		}
	}
	free((char*)p.centers);

	for(unsigned int i=0; i<HASHSIZE; i++){
		while(p.corners[i] != NULL){
			CORNERLIST* current = p.corners[i];
			p.corners[i] = current->next;
			free((char*)current);
		}
	}
	free((char*)p.corners);

	for(unsigned int i=0; i<2*HASHSIZE; i++){
		while(p.edges[i] != NULL){
			EDGELIST* current = p.edges[i];
			p.edges[i] = current->next;
			free((char*)current);
		}
	}
    free((char*)p.edges);

	//out
	int vertex_N = gvertices.count;
	int face_N = gntris;
	if(FILEOUT){
		fclose(ver_file);
		fclose(tri_file);
		ver_file = fopen(VERFILE, "r");
		tri_file = fopen(TRIFILE, "r");
		
		FILE* mesh_file = fopen(MESHFILE, "w");
		fprintf(mesh_file, "%d\n%d\n", vertex_N, face_N);
		for(int i=0; i<vertex_N; i++){
			//float *v = mesh->vertex[i];
			//fscanf(ver_file, "%f %f %f", &v[0], &v[1], &v[2]);

			float vx, vy, vz;
			fscanf(ver_file, "%f %f %f", &vx, &vy, &vz);
			fprintf(mesh_file, "%f %f %f\n", vx, vy, vz);
		}
		for(int i=0; i<face_N; i++){
			//mesh->setPolygonCount(i, 3);
			//int* f = mesh->face[i];
			//fscanf(tri_file, "%d %d %d", &f[0], &f[1], &f[2]);

			int i0, i1, i2;
			fscanf(tri_file, "%d %d %d", &i0, &i1, &i2);
			fprintf(mesh_file, "3 %d %d %d\n", i0, i1, i2);
		}
		fclose(ver_file);
		fclose(tri_file);

		fclose(mesh_file);
	}
	else{
		mesh->setFaceCount(face_N);
		mesh->setVertexCount(vertex_N);
		for(int i=0; i<vertex_N; i++){
			VERTEX v;
			v = gvertices.ptr[i];
			mesh->vertex[i][0] = (float)v.position.x;
			mesh->vertex[i][1] = (float)v.position.y;
			mesh->vertex[i][2] = (float)v.position.z;
		}
		if (gvertices.ptr != NULL) 
    {
      //free((char *)gvertices.ptr); !!! fix this - segfaults
    }

		TRILIST* current = tris;
		for(int i=0; i<face_N; i++){
			mesh->setPolygonCount(i, 3);
			mesh->face[i][0] = current->i1;
			mesh->face[i][1] = current->i2;
			mesh->face[i][2] = current->i3;
			TRILIST* previous = current;
			current = current->next;
			delete previous;
		}
	}

    return NULL;
}

/* testface: given cube at lattice (i, j, k), and four corners of face,
 * if surface crosses face, compute other four corners of adjacent cube
 * and add new cube to cube stack */

void Bloomenthal::testface(int i, int j, int k, CUBE *old, 
						   int face, int c1, int c2, int c3, int c4, PROCESS *p)
{
	CUBE new1;
    CUBES *oldcubes = p->cubes;
    //CORNER *setcorner();
//    static int facebit[6] = {2, 2, 1, 1, 0, 0};
    int n, pos = old->corners[c1]->value > 0.0 ? 1 : 0;//, bit = facebit[face];

    /* test if no surface crossing, cube out of bounds, or already visited: */
    if ((old->corners[c2]->value > 0) == pos &&
	(old->corners[c3]->value > 0) == pos &&
	(old->corners[c4]->value > 0) == pos) return;
    //if (abs(i) > p->bounds || abs(j) > p->bounds || abs(k) > p->bounds) return;
	float size = p->size;
	if(i*size + p->start.x > p->boundsX1 
		|| i*size + p->start.x < p->boundsX2 
		|| j*size + p->start.y > p->boundsY1 
		|| j*size + p->start.y < p->boundsY2 
		|| k*size + p->start.z > p->boundsZ1 
		|| k*size + p->start.z < p->boundsZ2)
		return;
    if (setcenter(p->centers, i, j, k)) return;

    /* create new cube: */
    new1.i = i;
    new1.j = j;
    new1.k = k;
	/*
    for (n = 0; n < 8; n++) new1.corners[n] = NULL;
    new1.corners[FLIP(c1, bit)] = old->corners[c1];
    new1.corners[FLIP(c2, bit)] = old->corners[c2];
    new1.corners[FLIP(c3, bit)] = old->corners[c3];
    new1.corners[FLIP(c4, bit)] = old->corners[c4];*/
    for (n = 0; n < 8; n++)
	//if (new1.corners[n] == NULL)
	    new1.corners[n] = setcorner1(p, i+BIT(n,2), j+BIT(n,1), k+BIT(n,0));

    /*add cube to top of stack: */
    p->cubes = (CUBES *) mycalloc(1, sizeof(CUBES));
    p->cubes->cube = new1;
    p->cubes->next = oldcubes;
}

/* setcorner: return corner with the given lattice location
   set (and cache) its function value */
Bloomenthal::CORNER * Bloomenthal::setcorner1(PROCESS *p, int i, int j, int k)
{
	/* for speed, do corner value caching here */
    CORNER *c = (CORNER *) mycalloc(1, sizeof(CORNER));
    int index = HASH(i, j, k);
    CORNERLIST *l = p->corners[index];
    c->i = i; c->x = p->start.x+((double)i-.5)*p->size;
    c->j = j; c->y = p->start.y+((double)j-.5)*p->size;
    c->k = k; c->z = p->start.z+((double)k-.5)*p->size;
    for (; l != NULL; l = l->next)
	if (l->i == i && l->j == j && l->k == k) {
	    c->value = l->value;
	    return c;
	    }
    l = (CORNERLIST *) mycalloc(1, sizeof(CORNERLIST));
    l->i = i; l->j = j; l->k = k;
    l->value = c->value = function(c->x, c->y, c->z);
    l->next = p->corners[index];
    p->corners[index] = l;
    return c;
}

/* find: search for point with value of given sign (0: neg, 1: pos) */


Bloomenthal::TEST Bloomenthal::find1(int sign, PROCESS *p, double x, double y, double z)
{
	int i;
    TEST test;
    double range = p->size;
    test.ok = 1;
    for (i = 0; i < 100000; i++) {
	test.p.x = x+range*(RAND()-0.5);
	test.p.y = y+range*(RAND()-0.5);
	test.p.z = z+range*(RAND()-0.5);
	test.value = function(test.p.x, test.p.y, test.p.z);
	if (sign == (test.value > 0.0)) return test;
	range = range*1.0005; /* slowly expand search outwards */
    }
    test.ok = 0;
    return test;
}

/**** Tetrahedral Polygonization ****/

/* dotet: triangulate the tetrahedron
 * b, c, d should appear clockwise when viewed from a
 * return 0 if client aborts, 1 otherwise */

int Bloomenthal::dotet(CUBE *cube, int c1, int c2, int c3, int c4, PROCESS *p)
{
	CORNER *a = cube->corners[c1];
    CORNER *b = cube->corners[c2];
    CORNER *c = cube->corners[c3];
    CORNER *d = cube->corners[c4];
    int index = 0, e1, e2, e3, e4, e5, e6;
    bool apos, bpos, cpos, dpos;
    
    apos = a->value > 0.0;
    if (apos) index += 8;
    
    bpos = (b->value > 0.0);
    if (bpos) index += 4;
    
    cpos = (c->value > 0.0);
    if (cpos) index += 2;
    
    dpos = (d->value > 0.0);
    if (dpos) index += 1;

    /* index is now 4-bit number representing one of the 16 possible cases */
    if (apos != bpos) e1 = vertid(a, b, p);
    if (apos != cpos) e2 = vertid(a, c, p);
    if (apos != dpos) e3 = vertid(a, d, p);
    if (bpos != cpos) e4 = vertid(b, c, p);
    if (bpos != dpos) e5 = vertid(b, d, p);
    if (cpos != dpos) e6 = vertid(c, d, p);
    /* 14 productive tetrahedral cases (0000 and 1111 do not yield polygons */
    switch (index) {
	case 1:	 return triproc(e5, e6, e3, p->vertices);
	case 2:	 return triproc(e2, e6, e4, p->vertices);
	case 3:	 return triproc(e3, e5, e4, p->vertices) &&
			triproc(e3, e4, e2, p->vertices);
	case 4:	 return triproc(e1, e4, e5, p->vertices);
	case 5:	 return triproc(e3, e1, e4, p->vertices) &&
			triproc(e3, e4, e6, p->vertices);
	case 6:	 return triproc(e1, e2, e6, p->vertices) &&
			triproc(e1, e6, e5, p->vertices);
	case 7:	 return triproc(e1, e2, e3, p->vertices);
	case 8:	 return triproc(e1, e3, e2, p->vertices);
	case 9:	 return triproc(e1, e5, e6, p->vertices) &&
			triproc(e1, e6, e2, p->vertices);
	case 10: return triproc(e1, e3, e6, p->vertices) &&
			triproc(e1, e6, e4, p->vertices);
	case 11: return triproc(e1, e5, e4, p->vertices);
	case 12: return triproc(e3, e2, e4, p->vertices) &&
			triproc(e3, e4, e5, p->vertices);
	case 13: return triproc(e6, e2, e4, p->vertices);
	case 14: return triproc(e5, e3, e6, p->vertices);
    }
    return 1;
}

/**** Cubical Polygonization (optional) ****/

/* docube: triangulate the cube directly, without decomposition */

int Bloomenthal::docube(CUBE *cube, PROCESS *p)
{
	INTLISTS *polys;
    int i, index = 0;
    for (i = 0; i < 8; i++) if (cube->corners[i]->value > 0.0) index += (1<<i);
    for (polys = cubetable[index]; polys; polys = polys->next) {
	INTLIST *edges;
	int a = -1, b = -1, count = 0;
	for (edges = polys->list; edges; edges = edges->next) {
	    CORNER *c1 = cube->corners[corner1[edges->i]];
	    CORNER *c2 = cube->corners[corner2[edges->i]];
	    int c = vertid(c1, c2, p);
	    if (++count > 2 && ! triproc(a, b, c, p->vertices)) return 0;
	    if (count < 3) a = b;
	    b = c;
	}
    }
    return 1;
}

/* nextcwedge: return next clockwise edge from given edge around given face */

 int Bloomenthal::nextcwedge(int edge, int face)
{
	switch (edge) {
	case LB: return (face == L)? LF : BN;
	case LT: return (face == L)? LN : TF;
	case LN: return (face == L)? LB : TN;
	case LF: return (face == L)? LT : BF;
	case RB: return (face == R)? RN : BF;
	case RT: return (face == R)? RF : TN;
	case RN: return (face == R)? RT : BN;
	case RF: return (face == R)? RB : TF;
	case BN: return (face == B)? RB : LN;
	case BF: return (face == B)? LB : RF;
	case TN: return (face == T)? LT : RN;
	case TF: return (face == T)? RT : LF;
	default: return 0;
    }
}

/* otherface: return face adjoining edge that is not the given face */

int Bloomenthal::otherface(int edge, int face)
{
	int other = leftface[edge];
    return face == other? rightface[edge] : other;
}

/* makecubetable: create the 256 entry table for cubical polygonization */

void Bloomenthal::makecubetable()
{
	static int visited = 0;

	if (visited)
		return;

	visited = 1;

	int i, e, c, done[12], pos[8];
    for (i = 0; i < 256; i++) {
	for (e = 0; e < 12; e++) done[e] = 0;
	for (c = 0; c < 8; c++) pos[c] = BIT(i, c);
	for (e = 0; e < 12; e++)
	    if (!done[e] && (pos[corner1[e]] != pos[corner2[e]])) {
		INTLIST *ints = 0;
		INTLISTS *lists = (INTLISTS *) mycalloc(1, sizeof(INTLISTS));
		int start = e, edge = e;
		/* get face that is to right of edge from pos to neg corner: */
		int face = pos[corner1[e]]? rightface[e] : leftface[e];
		while (1) {
		    edge = nextcwedge(edge, face);
		    done[edge] = 1;
		    if (pos[corner1[edge]] != pos[corner2[edge]]) {
			INTLIST *tmp = ints;
			ints = (INTLIST *) mycalloc(1, sizeof(INTLIST));
			ints->i = edge;
			ints->next = tmp; /* add edge to head of list */
			if (edge == start) break;
			face = otherface(edge, face);
		    }
		}
		lists->list = ints; /* add ints to head of table entry */
		lists->next = cubetable[i];
		cubetable[i] = lists;
	    }
    }
}

/**** Storage ****/

/* mycalloc: return successful calloc or exit program */

char * Bloomenthal::mycalloc(int nitems, int nbytes)
{
	char *ptr = (char*)calloc(nitems, nbytes);
	if (ptr != NULL) return ptr;
	fprintf(stderr, "can't calloc %d bytes\n", nitems*nbytes);
	exit(1);
}

/* setcenter: set (i,j,k) entry of table[]
 * return 1 if already set; otherwise, set and return 0 */

int Bloomenthal::setcenter(CENTERLIST *table[], int i, int j, int k)
{
	int index = HASH(i, j, k);
    CENTERLIST *new1, *l, *q = table[index];
    for (l = q; l != NULL; l = l->next)
	if (l->i == i && l->j == j && l->k == k) return 1;
    new1 = (CENTERLIST *) mycalloc(1, sizeof(CENTERLIST));
    new1->i = i; new1->j = j; new1->k = k; new1->next = q;
    table[index] = new1;
    return 0;
}

/* setedge: set vertex id for edge */

void Bloomenthal::setedge(EDGELIST *table[], int i1, int j1, int k1, int i2, int j2, int k2, int vid)
{
	unsigned int index;
    EDGELIST *new1;
    if (i1>i2 || (i1==i2 && (j1>j2 || (j1==j2 && k1>k2)))) {
	int t=i1; i1=i2; i2=t; t=j1; j1=j2; j2=t; t=k1; k1=k2; k2=t;
    }
    index = HASH(i1, j1, k1) + HASH(i2, j2, k2);
    new1 = (EDGELIST *) mycalloc(1, sizeof(EDGELIST));
    new1->i1 = i1; new1->j1 = j1; new1->k1 = k1;
    new1->i2 = i2; new1->j2 = j2; new1->k2 = k2;
    new1->vid = vid;
    new1->next = table[index];
    table[index] = new1;	
}

/* getedge: return vertex id for edge; return -1 if not set */

int Bloomenthal::getedge(EDGELIST *table[], int i1, int j1, int k1, int i2, int j2, int k2)
{
	EDGELIST *q;
    if (i1>i2 || (i1==i2 && (j1>j2 || (j1==j2 && k1>k2)))) {
	int t=i1; i1=i2; i2=t; t=j1; j1=j2; j2=t; t=k1; k1=k2; k2=t;
    };
    q = table[HASH(i1, j1, k1)+HASH(i2, j2, k2)];
    for (; q != NULL; q = q->next)
	if (q->i1 == i1 && q->j1 == j1 && q->k1 == k1 &&
	    q->i2 == i2 && q->j2 == j2 && q->k2 == k2)
	    return q->vid;
    return -1;
}

/**** Vertices ****/


/* vertid: return index for vertex on edge:
 * c1->value and c2->value are presumed of different sign
 * return saved index if any; else compute vertex and save */

int Bloomenthal::vertid(CORNER *c1, CORNER *c2, PROCESS *p)
{
	VERTEX v;
    POINT1 a, b;
    int vid = getedge(p->edges, c1->i, c1->j, c1->k, c2->i, c2->j, c2->k);
    if (vid != -1) return vid;			     /* previously computed */
    a.x = c1->x; a.y = c1->y; a.z = c1->z;
    b.x = c2->x; b.y = c2->y; b.z = c2->z;

	double x1, y1, z1, x2, y2, z2, f1, f2;
	if(c1->value > 0){
		f1 = c1->value;
		x1 = c1->x; y1 = c1->y; z1 = c1->z;
		f2 = c2->value;
		x2 = c2->x; y2 = c2->y; z2 = c2->z;
	}
	else{
		f1 = c2->value;
		x1 = c2->x; y1 = c2->y; z1 = c2->z;
		f2 = c1->value;
		x2 = c1->x; y2 = c1->y; z2 = c1->z;
	}

	//Repeated Linear Interpolation
	double x, y, z, f;
	for(int i=0; i<10; i++){
		double w1 = fabs(f1);
		double w2 = fabs(f2);
		x = (w2*x1 + w1*x2)/(w1+w2);
		y = (w2*y1 + w1*y2)/(w1+w2);
		z = (w2*z1 + w1*z2)/(w1+w2);
		break;

		f = function(x,y,z);
		if(fabs(f) < 0.00001)
			break;
		if(f > 0){
			f1 = f;
			x1 = x;  y1 = y;  z1 = z;
		}
		else{
			f2 = f;
			x2 = x;  y2 = y;  z2 = z;
		}
	}
	v.position.x = x;
	v.position.y = y;
	v.position.z = z;

	

    //converge(&a, &b, c1->value, &v.position); /* position */
    //vnormal(&v.position, p, &v.normal);			   /* normal */
    addtovertices(&p->vertices, v);			   /* save vertex */
    vid = p->vertices.count-1;
    setedge(p->edges, c1->i, c1->j, c1->k, c2->i, c2->j, c2->k, vid);
    return vid;
}

/* addtovertices: add v to sequence of vertices */

void Bloomenthal::addtovertices(VERTICES *vertices, VERTEX v)
{
	if(FILEOUT){
		fprintf(ver_file, "%f %f %f\n", v.position.x, v.position.y, v.position.z);
		vertices->count++;
	}
	else{
		if(vertices->count == vertices->max) {
			int i;
			VERTEX *new1;
			vertices->max = vertices->count == 0 ? 10 : 2*vertices->count;
			new1 = (VERTEX *) mycalloc((unsigned) vertices->max, sizeof(VERTEX));
			for (i = 0; i < vertices->count; i++) new1[i] = vertices->ptr[i];
			if (vertices->ptr != NULL) free((char *)vertices->ptr);
			vertices->ptr = new1;
		}
		vertices->ptr[vertices->count++] = v;
	}
}

/* vnormal: compute unit length surface normal at point */

void Bloomenthal::vnormal(POINT1 *point, PROCESS *p, POINT1 *v)
{
	double f = function(point->x, point->y, point->z);
    v->x = function(point->x+p->delta, point->y, point->z)-f;
    v->y = function(point->x, point->y+p->delta, point->z)-f;
    v->z = function(point->x, point->y, point->z+p->delta)-f;
    f = sqrt(v->x*v->x + v->y*v->y + v->z*v->z);
    if (f != 0.0) {v->x /= f; v->y /= f; v->z /= f;}
}

/* converge: from two points of differing sign, converge to zero crossing */
void Bloomenthal::converge(POINT1 *p1, POINT1 *p2, double v, POINT1 *p)
{
	int i = 0;
    POINT1 pos, neg;
    if (v < 0) {
	pos.x = p2->x; pos.y = p2->y; pos.z = p2->z;
	neg.x = p1->x; neg.y = p1->y; neg.z = p1->z;
    }
    else {
	pos.x = p1->x; pos.y = p1->y; pos.z = p1->z;
	neg.x = p2->x; neg.y = p2->y; neg.z = p2->z;
    }
    while (1) {
	p->x = 0.5*(pos.x + neg.x);
	p->y = 0.5*(pos.y + neg.y);
	p->z = 0.5*(pos.z + neg.z);
	if (i++ == RES) return;
	if ((function(p->x, p->y, p->z)) > 0.0)
	     {pos.x = p->x; pos.y = p->y; pos.z = p->z;}
	else {neg.x = p->x; neg.y = p->y; neg.z = p->z;}
    }
}

double Bloomenthal::function(double x, double y, double z)
{
	return func->value(x,y,z);
}
