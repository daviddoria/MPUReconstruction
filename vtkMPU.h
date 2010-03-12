#ifndef __vtkMPU_h
#define __vtkMPU_h

#include "vtkPolyDataAlgorithm.h"

class vtkImageData;
class ImplicitPOU;


/** 
\class vtkMPU
\brief vtk wrapper class for the MPU reconstruction method (see reference below).
Default values should work fine! But if you are willing to tweal the parameters, have a look first at the original papers.



\ref Bibtex:
@article{MPU,
 author = {Ohtake, Yutaka and Belyaev, Alexander and Alexa, Marc and Turk, Greg and Seidel, Hans-Peter},
 title = {Multi-level partition of unity implicits},
 journal = {ACM Trans. Graph.},
 volume = {22},
 number = {3},
 year = {2003},
 issn = {0730-0301},
 pages = {463--470},
 doi = {http://doi.acm.org/10.1145/882262.882293},
 publisher = {ACM},
 address = {New York, NY, USA},
 }
*/
class vtkMPU : public vtkPolyDataAlgorithm 
{
public:
  vtkTypeRevisionMacro(vtkMPU,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  static vtkMPU *New();
	  
  vtkSetMacro( Support, float );
  vtkGetMacro( Support, float );
  
  vtkSetMacro( Box, float );
  vtkGetMacro( Box, float );
  
  vtkSetMacro( Lambda, float );
  vtkGetMacro( Lambda, float );
  
  vtkSetMacro( Nmin, float );
  vtkGetMacro( Nmin, float );
  
  vtkSetMacro( Edge, float );
  vtkGetMacro( Edge, float );
  
  vtkSetMacro( Corner, float );
  vtkGetMacro( Corner, float );
  
  vtkSetMacro( Grid, float );
  vtkGetMacro( Grid, float );
  
  vtkSetMacro( Level, float );
  vtkGetMacro( Level, float );
  
  vtkSetMacro( Error, float );
  vtkGetMacro( Error, float );
  
  vtkSetMacro( Iso, float );
  vtkGetMacro( Iso, float );
  
  vtkSetMacro( Sharp, bool );
  vtkGetMacro( Sharp, bool );
  vtkBooleanMacro( Sharp, bool );

  vtkSetMacro( Bloomenthal, bool );
  vtkGetMacro( Bloomenthal, bool );
  vtkBooleanMacro( Bloomenthal, bool );

  vtkImageData* GetImplicitImage();
  
protected:
  vtkMPU();
  ~vtkMPU();
  
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  void ComputeImplicitImage();

  ImplicitPOU* _func;
  vtkImageData* ImplicitImage;

  float min[3];
  float max[3];
  float Support;
  float Box;
  float Lambda;
  float Nmin;
  float Edge;
  float Corner;
  float Grid;
  float Level;
  float Error;
  float Iso;
  bool Sharp;
  bool Bloomenthal;


private:
  vtkMPU(const vtkMPU&);  // Not implemented.
  void operator=(const vtkMPU&);  // Not implemented.

};

#endif
