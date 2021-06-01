cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(Analysis)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

include_directories(/usr/local/include)

#Find ROOT
find_package(ROOT REQUIRED)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")

set(SCRIPTS
config.json
EnergyLoss.dat
ProtonInSi.dat
ProtonInC2D4.dat
Geometry/QQQ5ExcitationAngle.dat
Geometry/QQQ5Geometry.dat
Calibration/UpQQQ5EnCal.dat
Calibration/UpSX3pedestals.dat
Calibration/UpSX3gains.dat
Calibration/UpSX3EnCal.dat
Calibration/UpSX3PosCal.dat
Calibration/BB10EnCal.dat
Calibration/DownSX3pedestals.dat
Calibration/DownSX3gains.dat
Calibration/DownSX3EnCal.dat
Calibration/DownSX3PosCal.dat
    )

foreach(_script ${SCRIPTS})
    configure_file(
        ${PROJECT_SOURCE_DIR}/${_script}
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${_script}
        COPYONLY
    )
endforeach()

add_library(QQQ5GammaAnalysisLibrary SHARED src/Q5Gamma.cpp)
target_link_libraries(QQQ5GammaAnalysisLibrary ${ROOT_LIBRARIES} Spectrum)
add_executable(Q5Gamma src/AnalysisRun.cpp src/jsoncpp.cpp src/EffectiveThickness.cpp src/Utilities.cpp)
target_link_libraries(Q5Gamma QQQ5GammaAnalysisLibrary)

add_library(SX3GammaAnalysisLibrary SHARED src/SX3Gamma.cpp)
target_link_libraries(SX3GammaAnalysisLibrary ${ROOT_LIBRARIES} Spectrum)
add_executable(SX3Gamma src/AnalysisRun.cpp src/jsoncpp.cpp src/EffectiveThickness.cpp src/Utilities.cpp)
target_link_libraries(SX3Gamma SX3GammaAnalysisLibrary)

add_library(DownstreamAnalysisLibrary SHARED src/DownstreamAnalysis.cpp)
target_link_libraries(DownstreamAnalysisLibrary ${ROOT_LIBRARIES} Spectrum)
add_executable(Downstream src/AnalysisRun.cpp src/jsoncpp.cpp src/EffectiveThickness.cpp src/Utilities.cpp)
target_link_libraries(Downstream DownstreamAnalysisLibrary)






