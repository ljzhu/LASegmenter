add_library(PoissonSurface vtkPoissonReconstruction.cxx Factor.cpp Geometry.cpp MarchingCubes.cpp MultiGridOctest.cpp MultiGridOctreeData.cpp)
TARGET_LINK_LIBRARIES(PoissonSurface ${VTK_LIBRARIES})
