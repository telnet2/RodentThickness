CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(MarinaTools) 

FIND_PACKAGE(SlicerExecutionModel REQUIRED)
IF(SlicerExecutionModel_FOUND)
	INCLUDE(${SlicerExecutionModel_USE_FILE})
	INCLUDE(${GenerateCLP_USE_FILE})
ELSE(SlicerExecutionModel_FOUND)
	MESSAGE(FATAL_ERROR "SlicerExecutionModel not found. Please set SlicerExecutionModel_DIR")
ENDIF(SlicerExecutionModel_FOUND)


FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})



GENERATECLP(MeshPointsIntensitySampling.cxx MeshPointsIntensitySampling.xml)
ADD_EXECUTABLE(MeshPointsIntensitySampling MeshPointsIntensitySampling.cxx MeshPointsIntensitySamplingCLP.h vtkPolyDataToitkMesh.cxx vtkPolyDataToitkMesh.h vtkio.cpp vtkio.h)
TARGET_LINK_LIBRARIES(MeshPointsIntensitySampling ${ITK_LIBRARIES} vtkFiltering vtkIO vtkGraphics)

#GENERATECLP(vtkTools.cxx vtkTools.xml)
ADD_EXECUTABLE(vtktools vtkTools.cxx piOptions.cpp)
TARGET_LINK_LIBRARIES(vtktools vtkFiltering vtkIO vtkGraphics vtkHybrid)


GENERATECLP(BinaryToDist.cxx BinToDistPreprocessing.xml)
ADD_EXECUTABLE(BinaryToDistanceMap BinaryToDist.cxx BinToDistPreprocessingCLP.h)
TARGET_LINK_LIBRARIES(BinaryToDistanceMap ${ITK_LIBRARIES})

    
ADD_EXECUTABLE(kmesh kmesh.h kmesh.cpp kimage.h kimage.cpp kstreamtracer.h kstreamtracer.cpp piOptions.h piOptions.cpp vtkio.h vtkio.cpp)
TARGET_LINK_LIBRARIES(kmesh ${VTK_LIBRARIES} ${ITK_LIBRARIES})

INCLUDE_DIRECTORIES(muparser)
FILE(GLOB muparser_SRC "muparser/*.cpp")
ADD_EXECUTABLE(kcalc kcalc.cpp PixelMathImageFilter.h PixelMathImageFilter.hxx piOptions.h piOptions.cpp piImageIO.h ${muparser_SRC})
TARGET_LINK_LIBRARIES(kcalc ${VTK_LIBRARIES} ${ITK_LIBRARIES})
