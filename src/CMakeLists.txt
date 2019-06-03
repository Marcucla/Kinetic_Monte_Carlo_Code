CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
PROJECT(DislocationKMC)

OPTION(ENABLE_DEBUG_CHECKS "Enable sanity checks that slow down the code a bit." "ON")

# Find the Boost library.
SET(Boost_NO_BOOST_CMAKE ON)
FIND_PACKAGE(Boost COMPONENTS program_options system REQUIRED)
IF(NOT Boost_FOUND)
	MESSAGE(FATAL_ERROR "Boost library not found. Reason: ${Boost_ERROR_REASON}")
ENDIF(NOT Boost_FOUND)
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

SET(SourceFiles
	util/linalg/LinAlg.cpp
	util/linalg/Tensor.cpp
	util/linalg/Matrix3.cpp
	util/Debug.cpp
	context/Context.cpp
	simulation/Output.cpp
	simulation/Initialization.cpp
	simulation/MonteCarlo.cpp
	simulation/KinkPairNucleation.cpp
	dislocations/Dislocations.cpp
	dislocations/Cleanup.cpp
	dislocations/KinkCreation.cpp
	dislocations/KinkMotion.cpp
	dislocations/CrossKinks.cpp
	dislocations/Output.cpp
	dislocations/ElasticInteraction.cpp
	pointdefects/PointDefects.cpp
	pointdefects/Output.cpp
	main/Main.cpp
)

IF(ENABLE_DEBUG_CHECKS)
	SET(BASE_COMPILE_DEFINITIONS "DEBUG_SIMULATION")
ENDIF()

SET(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -std=c++0x")

###################################################################
# Build executable.
###################################################################
ADD_EXECUTABLE(DislocationKMC ${SourceFiles})
TARGET_LINK_LIBRARIES(DislocationKMC ${Boost_LIBRARIES})
SET_TARGET_PROPERTIES(DislocationKMC PROPERTIES COMPILE_DEFINITIONS "${BASE_COMPILE_DEFINITIONS}")
