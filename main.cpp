#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>

#include "vtkMPU.h"

int main( int argc, char *argv[])
{
  if( argc != 4 )
    {
    std::cerr <<"MPU(.exe) takes 3 arguments" <<std::endl;
    std::cerr <<"1-Input filename" <<std::endl;
    std::cerr <<"2-Level" <<std::endl;
    std::cerr <<"3-Output filename" <<std::endl;
    return EXIT_FAILURE;
    }

	vtkSmartPointer<vtkXMLPolyDataReader> reader = 
      vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(argv[1]);
  reader->Update();
  
  vtkSmartPointer<vtkMPU> mpu = 
      vtkSmartPointer<vtkMPU>::New();
  mpu->SetInputConnection(reader->GetOutputPort());
  mpu->SetLevel( atoi( argv[2] ) );
  mpu->Update();
  
  vtkSmartPointer<vtkXMLPolyDataWriter> writer = 
      vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName( argv[3] );
  writer->SetInputConnection(mpu->GetOutputPort());
  writer->Write();
  

}
