#-------------------------------------------------------------------------------
# This file is part of the RBVMS application.
#
# RBVMSis free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.
#----------------------------------------------------------------------
cmake_minimum_required(VERSION 3.8)
include(ExternalProject)
project(RBVMS} LANGUAGES CXX C)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)

#-------------------------------------------------------------------------------
# Add external
#-------------------------------------------------------------------------------
set(EXT_DIR ${PROJECT_BINARY_DIR}/external)
include_directories(${EXT_DIR}/include)

#-------------------------------------------------------------------------------
# Add MPI
#-------------------------------------------------------------------------------
find_package(MPI REQUIRED)

set(MPICXX mpicxx)
set(MPI_CXX_INCLUDE_DIRS ${MPI_CXX_INCLUDE_PATH})
set(CMAKE_CXX_COMPILER ${MPICXX})

#-------------------------------------------------------------------------------
# Add METIS
#-------------------------------------------------------------------------------
find_package(METIS REQUIRED)

#-------------------------------------------------------------------------------
# Add hypre
#-------------------------------------------------------------------------------
set(HYPRE_REPO "git@github.com:hypre-space/hypre.git" 
    CACHE STRING "Location of the hypre repository")

ExternalProject_Add(hypre
                    PREFIX external
                    CONFIGURE_COMMAND ${CMAKE_COMMAND} "${EXT_DIR}/src/hypre/src" 
                                      "-DHYPRE_INSTALL_PREFIX=${EXT_DIR}"
                    GIT_REPOSITORY ${HYPRE_REPO})

set(HYPRE_LIB ${EXT_DIR}/lib/libHYPRE.a)

#-------------------------------------------------------------------------------
# Add mfem
#-------------------------------------------------------------------------------s
set(MFEM_REPO "git@github.com:mfem/mfem.git"
    CACHE STRING "Location of the mfem repository")

ExternalProject_Add(mfem
                    PREFIX external
                    DEPENDS hypre
                    CMAKE_ARGS -DMFEM_USE_MPI=Yes -DHYPRE_DIR=${EXT_DIR} -DCMAKE_INSTALL_PREFIX=${EXT_DIR}
                    GIT_REPOSITORY ${MFEM_REPO})

set(MFEM_LIB ${EXT_DIR}/lib/libmfem.a)

#-------------------------------------------------------------------------------
# Add rbvms
#-------------------------------------------------------------------------------
add_subdirectory(src)
