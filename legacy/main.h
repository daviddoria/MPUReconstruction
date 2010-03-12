
// Standard library
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <cstring>
/// GDCM library
//#include "gdcmFile.h"
//#include "gdcmFileHelper.h"
//#include "gdcmCommon.h"
//#include "gdcmDebug.h"
//#include "gdcmDocEntry.h"
//#include "gdcmSeqEntry.h"
//#include "gdcmContentEntry.h"
//#include "gdcmValEntry.h"
//#include "gdcmBinEntry.h"
//#include "gdcmSQItem.h"
/// VTK library

#include "vtkPoints.h"
#include "vtkIntArray.h"
#include "vtkFloatArray.h"
#include "vtkPolygon.h"
#include "vtkUnsignedCharArray.h"
#include "vtkPointData.h"
#include "vtkCellArray.h"
#include "vtkIdTypeArray.h"
#include "vtkDataArray.h"
#include "vtkImageData.h"
#include "vtkPolyDataNormals.h"
#include "vtkCellData.h"

#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkLookupTable.h"
#include "vtkImageMapToColors.h"
#include "vtkImageActor.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkScalarBarActor.h"
#include "vtkCamera.h"
#include "vtkProperty.h"
#include "vtkImageMapper.h"
#include "vtkVolumeMapper.h"
#include "vtkVolume.h"
#include "vtkVolumeRayCastCompositeFunction.h"
#include "vtkVolumeRayCastMapper.h"
#include "vtkXMLImageDataWriter.h"
#include "vtkXMLPolyDataWriter.h"

#include "vtkMath.h"
#include "vtkContourFilter.h"
#include "vtkImageGaussianSmooth.h"
#include "vtkImageGradient.h"

#include "ImplicitOctTree.h"
#include "PolygonalMesh.h"
#include "Polygonizer.h"
#include "PointSet.h"

#ifndef MIN
#define MIN(a,b)    ((a)>(b)?(b):(a))
#endif

#ifndef MAX
#define MAX(a,b)    ((a)>=(b)?(a):(b))
#endif
