CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(Shape2014Validation)

FIND_PACKAGE(STATISMO REQUIRED)

FIND_PACKAGE(ITK REQUIRED)
include (${ITK_USE_FILE})
	
FILE(GLOB INCLUDE_FILES *.h)
FILE(GLOB SOURCE_FILES  *.txx *.cpp)

SET (SOURCE_FILES 
	src/shapemodelvalidation.cpp 
	src/generalization.cpp 
	src/compactness.cpp 
	src/specificity.cpp 
	src/common.cpp 
	src/config.cpp
)

include_directories(${STATISMO_INCLUDE_DIRS} ${ITK_INCLUDE_DIRS})
link_directories( ${ITK_LIBRARY_DIRS})

add_executable (shapemodelValidation ${SOURCE_FILES})
target_link_libraries (shapemodelValidation ${ITK_LIBRARIES})
