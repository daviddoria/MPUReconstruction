#include "vtkMPUImplicitFunction.h"
#include "vtkObjectFactory.h"

vtkCxxRevisionMacro(vtkMPUImplicitFunction, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkMPUImplicitFunction);

vtkMPUImplicitFunction::vtkMPUImplicitFunction() : function( 0 )
{
}

vtkMPUImplicitFunction::~vtkMPUImplicitFunction()
{
}

void vtkMPUImplicitFunction::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

void vtkMPUImplicitFunction::SetMPUFunction( ImplicitFunction* iF )
{
  this->function = iF;
}

double vtkMPUImplicitFunction::EvaluateFunction( double iX[3] )
{
  if( this->function )
    {
    float x = static_cast< float >( iX[0] );
    float y = static_cast< float >( iX[1] );
    float z = static_cast< float >( iX[2] );
    return static_cast< double >( this->function->value( x, y, z ) );
    }
  else
    {
    std::cerr <<"this->function is NULL" <<std::endl;
    return 0.;
    }
}

void vtkMPUImplicitFunction::EvaluateGradient( double iX[3], double oG[3] )
{
  if( this->function )
    {
    float x = static_cast< float >( iX[0] );
    float y = static_cast< float >( iX[1] );
    float z = static_cast< float >( iX[2] );

    float g[3];
    this->function->gradient( g, x, y, z);

    oG[0] = static_cast< double >( g[0] );
    oG[1] = static_cast< double >( g[1] );
    oG[2] = static_cast< double >( g[2] );
    }
  else
    {
    std::cerr <<"this->function is NULL" <<std::endl;
    }
}