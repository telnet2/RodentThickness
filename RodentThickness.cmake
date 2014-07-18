CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

set(MODULE_NAME ${EXTENSION_NAME}) # Do not use 'project()'
set(MODULE_TITLE ${MODULE_NAME})

include(${CMAKE_CURRENT_SOURCE_DIR}/Common.cmake)


FIND_PACKAGE(SlicerExecutionModel NO_MODULE REQUIRED)
INCLUDE(${SlicerExecutionModel_USE_FILE})
INCLUDE(${GenerateCLP_USE_FILE})


FIND_PACKAGE(BatchMake NO_MODULE REQUIRED)
IF(BatchMake_FOUND)
	INCLUDE(${BatchMake_USE_FILE})
ELSE(BatchMake_FOUND)
	MESSAGE(FATAL_ERROR "BatchMake not found. Please set BatchMake_DIR")
ENDIF(BatchMake_FOUND)

#-----------------------------------------------------------------------------
find_package(ITK ${ITK_VERSION_MAJOR} NO_MODULE REQUIRED)

include(${ITK_USE_FILE})

FIND_PACKAGE(VTK NO_MODULE REQUIRED)
INCLUDE (${VTK_USE_FILE})

FIND_PACKAGE(Qt4 REQUIRED) #/tools/devel/linux/Qt4/Qt4.8.1/Qt4.8.1_linux64/bin/qmake
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${QT_INCLUDE_DIR})
INCLUDE(${QT_USE_FILE})
ADD_DEFINITIONS(-DQT_GUI_LIBS -DQT_CORE_LIB -DQT3_SUPPORT)


#======================================================================================
include(ExternalProject)
include(SlicerMacroEmptyExternalProject)
include(SlicerMacroCheckExternalProjectDependency)


set( proj ${LOCAL_PROJECT_NAME} )
option(USE_SPHARM-PDM "Build SPHARM-PDM" ON)
if(USE_SPHARM-PDM)

 set(${LOCAL_PROJECT_NAME}_DEPENDENCIES spharm-pdm)
endif()



include(${CMAKE_CURRENT_SOURCE_DIR}/SetExternalProjectOptions.cmake)



add_subdirectory(Applications)





