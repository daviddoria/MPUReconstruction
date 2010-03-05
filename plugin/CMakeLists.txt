FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES} ../)

ADD_PARAVIEW_PLUGIN(MPUReconstructionPlugin "1.0"
  SERVER_MANAGER_XML ../vtkMPUReconstruction.xml
SERVER_MANAGER_SOURCES ../vtkMPUReconstruction.cxx  
../main.cpp 
../vtkMPU.cxx
../ImplicitFunction.cpp
../Polygonizer.cpp
../PolygonalMesh.cpp
../OctTreeP.cpp
../ImplicitOctTree.cpp
../SVD.cpp
../Bloomenthal.cpp
../QuadricOEdge.cpp
../Quadric.cpp
../QuadricO.cpp
../QuadricOCorner.cpp
../AxisKdTree.cpp
)

