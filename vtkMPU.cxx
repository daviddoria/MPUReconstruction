#include "vtkMPU.h"

#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkDataObject.h"
#include "vtkSmartPointer.h"
#include "vtkTriangle.h"
#include "vtkCellArray.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"
#include "vtkCellData.h"

#include "PointSet.h"
#include "ImplicitOctTree.h"
#include "Polygonizer.h"

#include "vtkMPUImplicitFunction.h"
#include "vtkSampleFunction.h"
#include "vtkImageData.h"
#include "vtkMarchingCubes.h"

vtkCxxRevisionMacro(vtkMPU, "$Revision: 1.70 $");
vtkStandardNewMacro(vtkMPU);

vtkMPU::vtkMPU()
{
  // parameter setting
  Support = 0.75;
  Box = 1.2;
  Lambda = 0.1;
  Nmin = 15;
  Edge = 0.9;
  Corner = 0.7;
  Grid = 0.01;
  Level = 20;
  Error = 0.005;
  Iso = 0;
  Sharp=true;
  Bloomenthal = true;

  min[0] = min[1] = min[2] = 0.;
  max[0] = max[1] = max[2] = 0.;

  ImplicitImage = vtkImageData::New();
  _func = new ImplicitPOU;
}

vtkMPU::~vtkMPU()
{
  ImplicitImage->Delete();
  delete _func;
}

int vtkMPU::RequestData(vtkInformation *vtkNotUsed(request),
                                             vtkInformationVector **inputVector,
                                             vtkInformationVector *outputVector)
{
  
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
    
  // get the input and ouptut
  vtkPolyData *input = vtkPolyData::SafeDownCast(
      inInfo->Get(vtkDataObject::DATA_OBJECT()));
  
  vtkPolyData *output = vtkPolyData::SafeDownCast(
      outInfo->Get(vtkDataObject::DATA_OBJECT()));
    
  PointSet *_ps = new PointSet;
  _ps->setPointSize(input->GetNumberOfPoints());
  
  vtkSmartPointer<vtkFloatArray> normals = 
      vtkFloatArray::SafeDownCast(input->GetPointData()->GetNormals());
  
  cout << "There are " << input->GetNumberOfPoints() << " input points." << endl;
  for(vtkIdType i = 0; i < input->GetNumberOfPoints(); i++)
    {
    double p[3];
    input->GetPoint(i, p);
    //cout << "Point " << i << " : " << p[0] << " " << p[1] << " " << p[2] << endl;
    _ps->setPoint(i, p[0], p[1], p[2]);
    
    double n[3];
    normals->GetTuple(i, n);
    //cout << "Normal " << i << " : " << n[0] << " " << n[1] << " " << n[2] << endl;
    _ps->setNormal(i, n[0], n[1], n[2]);
    }
    
  _ps->bound(this->min, this->max);
  float mid[3];
  mid[0] = 0.5f*(this->min[0] + this->max[0]);
  mid[1] = 0.5f*(this->min[1] + this->max[1]);
  mid[2] = 0.5f*(this->min[2] + this->max[2]);
  float sizeX = this->max[0] - this->min[0];
  float sizeY = this->max[1] - this->min[1];
  float sizeZ = this->max[2] - this->min[2];
  //  float size = MAX(sizeX, MAX(sizeY, sizeZ))*0.5f*m_box;
  float size = std::max(sizeX, std::max(sizeY, sizeZ))*0.5f*Box;
  float maxC[3], minC[3];
  minC[0] = mid[0] - size;
  minC[1] = mid[1] - size;
  minC[2] = mid[2] - size;
  maxC[0] = mid[0] + size;
  maxC[1] = mid[1] + size;
  maxC[2] = mid[2] + size;

  _func->imp = new ImplicitOctTree(_ps, minC, maxC);
  _func->setMaxCellN(100000000);
  _func->imp->_support = Support;
  _func->imp->_Nmin = Nmin;
  _func->imp->_lambda = Lambda;
  _func->imp->_sharp = Sharp;
  _func->imp->_edgeT = Edge;
  _func->imp->_cornerT = Corner;
  _func->_max_level = Level;
  _func->_tol = Error*(float)sqrt(sizeX*sizeX + sizeY*sizeY + sizeZ*sizeZ);
  _func->_iso = Iso*(float)sqrt(sizeX*sizeX + sizeY*sizeY + sizeZ*sizeZ);
  // polygonization
  int N = _ps->_pointN;
  Polygonizer poly;
  poly.func = _func;
  float box[] = {this->max[0], this->min[0], this->max[1], this->min[1], this->max[2], this->min[2]};
  
  float boxC[] = {0.5f*(box[0]+box[1]), 0.5f*(box[2]+box[3]), 0.5f*(box[4]+box[5])};
  for(int i=0; i<N; i++)
    {
    float *p = _ps->_point[i];
    if(p[0] < box[0] && p[0] > box[1] &&
       p[1] < box[2] && p[1] > box[3] &&
       p[2] < box[4] && p[2] > box[5])
      {
      boxC[0] = p[0];
      boxC[1] = p[1];
      boxC[2] = p[2];

      break;
      }
    }
  _func->count = 0;

  if( this->Bloomenthal )
    {
    PolygonalMesh* meshMPU = poly.bloomenthal(Grid*(float)sqrt(sizeX*sizeX + sizeY*sizeY + sizeZ*sizeZ),boxC,box);

    // Visualize
    vtkSmartPointer<vtkPolyData> outputMesh = 
        vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkCellArray> outputTriangles = 
        vtkSmartPointer<vtkCellArray>::New();
    
    vtkSmartPointer<vtkPoints> outputPoints = 
        vtkSmartPointer<vtkPoints>::New();
    
    vtkIdType pts[3];
    
    vtkSmartPointer<vtkFloatArray> outputNormals = 
        vtkSmartPointer<vtkFloatArray>::New();
    // set vertex
    int numVertex = meshMPU->vertex_N;
    cout << "There are " << numVertex << " vertices." << endl;
    outputPoints->SetNumberOfPoints(numVertex);
    float (*vertexMPU)[3] = meshMPU->vertex;
    
    for (int i = 0; i < numVertex; i++)
      {
      outputPoints->SetPoint(i, vertexMPU[i]);
      }
      
    outputMesh->SetPoints(outputPoints);
    
    // set triangle
    int numTriangle = meshMPU->face_N;
    cout << "There are " << numTriangle << " triangles." << endl;
    int **face = meshMPU->face;
//    int *poly_N = meshMPU->poly_N; // store num of vertex for ith face, shouldn't be a constant like 3?
    
    for(int i = 0; i < numTriangle; i++)
      {
      //pts->SetComponent(0,0,face[i][0]);
      //pts->SetComponent(1,0,face[i][1]);
      //pts->SetComponent(2,0,face[i][2]);
      pts[0]=face[i][0];  pts[1]=face[i][1];  pts[2]=face[i][2];
      outputTriangles->InsertNextCell(3,pts);
      }
      
    outputMesh->SetPolys(outputTriangles);
    // set normal
    meshMPU->computeFaceNormal();
    float (*normalMPU)[3] = meshMPU->normal_f;
    outputNormals->SetNumberOfComponents(3);
    outputNormals->SetNumberOfTuples(numTriangle);
    
    for (int i=0; i<numTriangle; i++)
      {
      outputNormals->SetComponent(i,0,normalMPU[i][0]);
      outputNormals->SetComponent(i,1,normalMPU[i][1]);
      outputNormals->SetComponent(i,2,normalMPU[i][2]);
      }
      
    outputMesh->GetCellData()->SetNormals(outputNormals);
      
    output->ShallowCopy(outputMesh);
    }
  else // Marching Cubes
    {
    ComputeImplicitImage( );

    vtkSmartPointer< vtkMarchingCubes > mc =
      vtkSmartPointer< vtkMarchingCubes >::New();
    mc->SetInput( this->ImplicitImage );
    mc->SetNumberOfContours(1);
    mc->ComputeScalarsOn();
    mc->ComputeGradientsOn();
    mc->ComputeNormalsOn();
    mc->SetValue(0, 0.);
    mc->Update();

    output->ShallowCopy( mc->GetOutput() );
    }

  delete _ps;
  
  return 1;
}

void vtkMPU::ComputeImplicitImage( )
{
  vtkSmartPointer< vtkMPUImplicitFunction > implicit =
    vtkSmartPointer< vtkMPUImplicitFunction >::New();
  implicit->SetMPUFunction( this->_func );

  int dimension = static_cast< int >( 1. / this->Grid );
  vtkSmartPointer< vtkSampleFunction > sample =
    vtkSmartPointer< vtkSampleFunction >::New();
  sample->SetImplicitFunction( implicit );
  sample->SetOutputScalarTypeToFloat();
  sample->SetSampleDimensions( dimension, dimension, dimension );
  sample->SetModelBounds( this->min[0], this->max[0], this->min[1], this->max[1], this->min[2], this->max[2] );
  sample->Update();

  this->ImplicitImage->ShallowCopy( sample->GetOutput() );
}

vtkImageData* vtkMPU::GetImplicitImage()
{
  if( this->Bloomenthal )
    {
    ComputeImplicitImage( );
    }

  return this->ImplicitImage;
}

//----------------------------------------------------------------------------
void vtkMPU::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

