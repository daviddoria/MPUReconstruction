cmake_minimum_required(VERSION 2.6)

PROJECT(MPU)

FIND_PACKAGE (VTK)
IF ( VTK_FOUND)
	INCLUDE( ${USE_VTK_FILE} )
ENDIF( VTK_FOUND)

# only to get the svd with right license!
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)
########################################

SET(BUILD_PARAVIEW_PLUGIN OFF CACHE BOOL "Build Paraview plugin?")

if(BUILD_PARAVIEW_PLUGIN)
FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES} ../)

ADD_PARAVIEW_PLUGIN(MPUReconstructionPlugin "1.0"
  SERVER_MANAGER_XML vtkMPU.xml
SERVER_MANAGER_SOURCES vtkMPU.cxx vtkMPUImplicitFunction.cxx
SERVER_SOURCES vtkMPU.cxx ImplicitFunction.cpp Polygonizer.cpp PolygonalMesh.cpp OctTreeP.cpp ImplicitOctTree.cpp SVD.cpp Bloomenthal.cpp QuadricOEdge.cpp Quadric.cpp QuadricO.cpp QuadricOCorner.cpp AxisKdTree.cpp
)

endif(BUILD_PARAVIEW_PLUGIN)

INCLUDE_DIRECTORIES(  
  ${Braude_SOURCE_DIR}
  ${ITK_INCLUDE_DIRS}
)
 
ADD_EXECUTABLE(MPU
main.cpp 
vtkMPU.cxx
ImplicitFunction.cpp
Polygonizer.cpp
PolygonalMesh.cpp
OctTreeP.cpp
ImplicitOctTree.cpp
Bloomenthal.cpp
QuadricOEdge.cpp
Quadric.cpp
QuadricO.cpp
QuadricOCorner.cpp
AxisKdTree.cpp
vtkMPUImplicitFunction.cxx
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
ITKNumerics
)  


