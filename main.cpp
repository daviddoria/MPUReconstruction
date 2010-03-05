#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>

#include "vtkMPU.h"

int main( int argc, char *argv[])
{
	vtkSmartPointer<vtkXMLPolyDataReader> reader = 
      vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(argv[1]);
  reader->Update();
  
  vtkSmartPointer<vtkMPU> mpu = 
      vtkSmartPointer<vtkMPU>::New();
  mpu->SetInputConnection(reader->GetOutputPort());
  mpu->Update();
  
  vtkSmartPointer<vtkXMLPolyDataWriter> writer = 
      vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(argv[2]);
  writer->SetInputConnection(mpu->GetOutputPort());
  writer->Write();
  

}
