# This file is part of the RBVMS application.
#
# RBVMSis free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.
#----------------------------------------------------------------------
cmake_minimum_required(VERSION 3.8)
project(RBVMS LANGUAGES CXX C)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../config)

#-------------------------------------------------------------------------------
# Add MPI
#-------------------------------------------------------------------------------
find_package(MPI REQUIRED)

set(MPICXX mpicxx)
set(MPI_CXX_INCLUDE_DIRS ${MPI_CXX_INCLUDE_PATH})
set(CMAKE_CXX_COMPILER ${MPICXX})

#-------------------------------------------------------------------------------
# Add mfem
#-------------------------------------------------------------------------------
find_package(mfem REQUIRED NAMES MFEM HINTS "${MFEM_DIR}")

message(STATUS "Found mfem: in ${MFEM_DIR} (version ${MFEM_VERSION})")
message(STATUS "MFEM_CXX_FLAGS = ${MFEM_CXX_FLAGS}")
message(STATUS "MFEM_INCLUDE_DIRS = ${MFEM_INCLUDE_DIRS}")
message(STATUS "MFEM_LIBRARIES = ${MFEM_LIBRARIES}")

include_directories(${MFEM_INCLUDE_DIRS})

#-------------------------------------------------------------------------------
# Create header
#-------------------------------------------------------------------------------
add_custom_command(OUTPUT src/buildInfo.hpp
                   PRE_BUILD
                   COMMAND /bin/sh "${PROJECT_SOURCE_DIR}/config/buildInfo.sh" "${PROJECT_SOURCE_DIR}" "${EXT_DIR}"
                   COMMENT "Getting git info")

#-------------------------------------------------------------------------------
# Build executable
#-------------------------------------------------------------------------------
add_executable(rbvms rbvms.cpp
               monitor.cpp tau.cpp weakform.cpp
               monitor.hpp tau.hpp weakform.hpp)

target_link_libraries(rbvms PRIVATE ${MFEM_LIBRARIES})
install(TARGETS rbvms DESTINATION bin)

#-------------------------------------------------------------------------------
# Add tests
#-------------------------------------------------------------------------------
#add_test(NAME navsto_np=4
#         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${RBVMS_MPI_NP}
#         ${MPIEXEC_PREFLAGS} $<TARGET_FILE:navsto> -no-vis
#         ${MPIEXEC_POSTFLAGS}
#        -m ${PROJECT_SOURCE_DIR}/data/square-nurbs.mesh -r 1 -o 2)

