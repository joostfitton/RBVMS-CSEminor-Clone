#-------------------------------------------------------------------------------
#
#
#
#-------------------------------------------------------------------------------
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
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/hypre/src/CMakeLists.txt")
    message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()
add_subdirectory(hypre/src)

#-------------------------------------------------------------------------------
# Add mfem
#-------------------------------------------------------------------------------
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/mfem/CMakeLists.txt")
    message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()
add_subdirectory(mfem)