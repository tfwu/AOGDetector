# required cmake version
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(RGM)


# The build options
option(RGM_USE_DOUBLE "Use doubles instead of floats (slower, uses twice more memory, and the increase in precision is not necessarily useful)." OFF)
option(RGM_USE_FELZENSZWALB_HOG_FEATURES "Use the implementation adapted from VOC-Realease 4 and 5." ON)
option(RGM_USE_EXTRA_FEATURES "Use extra features (LBP + color) in addition to HOG." OFF)

# Select a default build configuration if none was chosen
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Debug"
		CACHE STRING "Choose the type of build, options are: None (CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

# set the directory where the binaries will be stored
set(EXECUTABLE_OUTPUT_PATH  ${PROJECT_BINARY_DIR}/bin
  CACHE PATH  "Directory where all executables will be stored"
)

# set the directory where the libraries will be stored
set(LIBRARY_OUTPUT_PATH  ${PROJECT_BINARY_DIR}/lib
  CACHE PATH  "Directory where all the libraries will be stored"
)

# Define the options
if(RGM_USE_DOUBLE)
  message("Use doubles instead of floats.")
  add_definitions(-DRGM_USE_DOUBLE)
endif()

if(RGM_USE_FELZENSZWALB_HOG_FEATURES)
  message("Use the implementation of computing HOG adapted from VOC-Realease 4 and 5.")
  add_definitions(-DRGM_USE_FELZENSZWALB_HOG_FEATURES)
endif()

if(RGM_USE_EXTRA_FEATURES)
  message("Use extra features (LBP + color) in addition to HOG.")
  add_definitions(-DRGM_USE_EXTRA_FEATURES)
endif()


# find OpenCV
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
#message("${OpenCV_LIBS}")


# find C++ Boost
set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME     OFF)
find_package(Boost 1.55.0 COMPONENTS system filesystem log serialization thread REQUIRED) 
#add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
include_directories(${Boost_INCLUDE_DIRS})
#message("${Boost_LIBRARIES}")


# find Eigen package
#file(GLOB EIGEN_ARCHIVE "eigen3" PATH /usr/include)
set(EIGEN_ARCHIVE /usr/include/eigen3)
find_path(EIGEN_INCLUDE_DIR Eigen ${EIGEN_ARCHIVE} .)
if(NOT EIGEN_INCLUDE_DIR)
  message(FATAL_ERROR "Could not find Eigen matrix library.")
endif()
include_directories(${EIGEN_INCLUDE_DIR})
add_definitions(-DEIGEN_USE_NEW_STDVECTOR
                -DEIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET)


# find FFTW package
find_path(FFTW3_INCLUDE_DIR fftw3.h)
if(FFLD_HOGPYRAMID_DOUBLE)
  find_library(FFTW3_LIBRARIES fftw3)
else()
  find_library(FFTW3_LIBRARIES fftw3f)
endif()
if(NOT FFTW3_INCLUDE_DIR OR NOT FFTW3_LIBRARIES)
  MESSAGE(FATAL_ERROR "Could not find fftw3.")
endif()
include_directories(${FFTW3_INCLUDE_DIR})

# matio lib
option(USE_MATLAB_IO "Use MatIO for interfacing Matlab mat files." ON)
if(USE_MATLAB_IO)
        add_definitions(-DRGM_USE_MATLAB)
        find_library(MATIO_LIB NAMES   matio)
        find_library(LIB_HDF5_HL NAMES hdf5_hl)
        find_library(LIB_HDF5    NAMES hdf5)
	if(NOT MATIO_LIB)
  		message(FATAL_ERROR "Could not find matio. Please download (http://sourceforge.net/projects/matio/) and install it.")
	endif()
endif()


# all external libs
set(ALL_EXTERNAL_LIBS	
	${OpenCV_LIBS}
	${Boost_LIBRARIES}
	${BLAS_LIBRARIES}
	${FFTW3_LIBRARIES}        
    ${LIB_HDF5_HL}
    ${LIB_HDF5}
    ${MATIO_LIB}
)
#message("${ALL_EXTERNAL_LIBS}")

# RGM
set(RGM_INCLUDE_DIRS
	${PROJECT_SOURCE_DIR}/src/core
	${PROJECT_SOURCE_DIR}/external/
)
include_directories(${RGM_INCLUDE_DIRS})

# CXX Flags
set(COMMON_FLAGS "${COMMON_FLAGS} -Wall -std=c++11 -fPIC") 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMMON_FLAGS}")
set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   ${COMMON_FLAGS}")

# find OpenMP, Not required, but strongly recommended on multi-core systems
find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# get all source file
file(GLOB_RECURSE	CXX_SRCS	${PROJECT_SOURCE_DIR}/src/core/*.cpp)
file(GLOB_RECURSE	HXX_SRCS	${PROJECT_SOURCE_DIR}/src/core/*.hpp)

# RGM  library
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
	add_library(RGM_Cored ${CXX_SRCS} ${HXX_SRCS})
	set(RGM_LIBRARIES RGM_Cored)
else()
	add_library(RGM_Core ${CXX_SRCS} ${HXX_SRCS})
	set(RGM_LIBRARIES RGM_Core)
endif()

# excutable
file(GLOB_RECURSE	MAIN_SRCS	${PROJECT_SOURCE_DIR}/src/main/*.cpp)
foreach (m ${MAIN_SRCS})
	get_filename_component(f ${m} NAME_WE)
	if (CMAKE_BUILD_TYPE STREQUAL "Debug")
            add_executable(${f}d ${m})
            target_link_libraries(${f}d	${RGM_LIBRARIES}	${ALL_EXTERNAL_LIBS})            
	else()
            add_executable(${f} ${m})
            target_link_libraries(${f}	${RGM_LIBRARIES}	${ALL_EXTERNAL_LIBS})            
	endif()
endforeach()


