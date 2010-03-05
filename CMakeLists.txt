cmake_minimum_required(VERSION 2.6)

PROJECT(MPU)

FIND_PACKAGE (VTK)
IF ( VTK_FOUND)
	INCLUDE( ${USE_VTK_FILE} )
ENDIF( VTK_FOUND)

INCLUDE_DIRECTORIES(  
  ${Braude_SOURCE_DIR}
)
 
ADD_EXECUTABLE(MPU
main.cpp 
vtkMPU.cxx
ImplicitFunction.cpp
Polygonizer.cpp
PolygonalMesh.cpp
OctTreeP.cpp
ImplicitOctTree.cpp
SVD.cpp
Bloomenthal.cpp
QuadricOEdge.cpp
Quadric.cpp
QuadricO.cpp
QuadricOCorner.cpp
AxisKdTree.cpp
)

TARGET_LINK_LIBRARIES(MPU
vtkRendering 
vtkGraphics 
vtkHybrid 
vtkImaging 
vtkIO 
vtkFiltering 
vtkCommon
vtkVolumeRendering
)  


