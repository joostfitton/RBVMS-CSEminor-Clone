#-------------------------------------------------------------------------------
# This file is part of the RBVMS application.
#
# RBVMSis free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.
#----------------------------------------------------------------------
cmake_minimum_required(VERSION 3.8)

set(PROJECT_NAME RBVMS)
project(${PROJECT_NAME}
  LANGUAGES CXX C)

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/install")
endif ()

#-------------------------------------------------------------------------------
# Get git submodules
#-------------------------------------------------------------------------------
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

#-------------------------------------------------------------------------------
# Add hypre
#-------------------------------------------------------------------------------
message(STATUS "Configure Hypre")
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/hypre/src/CMakeLists.txt")
    message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()
add_subdirectory(hypre/src)

#-------------------------------------------------------------------------------
# Add mfem
#-------------------------------------------------------------------------------
message(STATUS "Configure MFEM")
option(MFEM_USE_MPI "Enable MPI parallel build" ON)
set(HYPRE_DIR "${PROJECT_SOURCE_DIR}/install" CACHE PATH "")

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/mfem/CMakeLists.txt")
    message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()
add_subdirectory(mfem)

#-------------------------------------------------------------------------------
# Add rbvms
#-------------------------------------------------------------------------------
message(STATUS "Configure RBVMS")
add_subdirectory(src)