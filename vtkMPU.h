#ifndef __vtkMPU_h
#define __vtkMPU_h

#include "vtkPolyDataAlgorithm.h"

class vtkMPU : public vtkPolyDataAlgorithm 
{
public:
  vtkTypeRevisionMacro(vtkMPU,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  static vtkMPU *New();
	  
  
protected:
  vtkMPU();
  
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
  vtkMPU(const vtkMPU&);  // Not implemented.
  void operator=(const vtkMPU&);  // Not implemented.

};

#endif
