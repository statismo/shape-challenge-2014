CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(Shape2014Validation)

FIND_PACKAGE(STATISMO REQUIRED)

FIND_PACKAGE(ITK REQUIRED)
include (${ITK_USE_FILE})
	
FILE(GLOB INCLUDE_FILES *.h)
FILE(GLOB SOURCE_FILES  *.txx *.cpp)

SET (SOURCE_FILES_COMMON
        src/generalization.cpp
        src/compactness.cpp
        src/specificity.cpp
        src/utils.cpp
        src/config.cpp
        src/correspondence.cpp
        src/distances.cpp
)
SET (SOURCE_FILES_CONSOLE
	src/shapemodelvalidation.cpp 
        ${SOURCE_FILES_COMMON}
)

SET (SOURCE_FILES_VSD
        src/shapemodelvalidation_vsd.cpp
        ${SOURCE_FILES_COMMON}
)

include_directories(${STATISMO_INCLUDE_DIRS} ${ITK_INCLUDE_DIRS}/itkhdf5 ${ITK_INCLUDE_DIRS}/itkhdf5/cpp ${ITK_INCLUDE_DIRS})

link_directories( ${ITK_LIBRARY_DIRS})

add_executable (shapemodelValidation ${SOURCE_FILES_CONSOLE})
target_link_libraries (shapemodelValidation ${ITK_LIBRARIES})

add_executable (shapemodelValidationVSD ${SOURCE_FILES_VSD})
target_link_libraries (shapemodelValidationVSD ${ITK_LIBRARIES})
