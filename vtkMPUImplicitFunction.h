#ifndef __vtkMPUImplicitFunction_h
#define __vtkMPUImplicitFunction_h

#include "ImplicitOctTree.h"
#include "vtkImplicitFunction.h"

class vtkMPUImplicitFunction : public vtkImplicitFunction
{
  public:
    static vtkMPUImplicitFunction *New();

    void PrintSelf(ostream& os, vtkIndent indent);
    vtkTypeRevisionMacro(vtkMPUImplicitFunction,vtkImplicitFunction);

    virtual void SetMPUFunction( ImplicitFunction* iF );

    virtual double EvaluateFunction( double x[3] );
    virtual void EvaluateGradient( double x[3], double g[3] );

  protected:
    vtkMPUImplicitFunction();
    virtual ~vtkMPUImplicitFunction();

    ImplicitFunction* function;

  private:
    vtkMPUImplicitFunction( const vtkMPUImplicitFunction& );
    void operator = ( const vtkMPUImplicitFunction& );
};
#endif

