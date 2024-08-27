# This cmake script is inspired on the add_mfem_miniapp in MfemCmakeUtilities

#-------------------------------------------------------------------------------
# Set variables
#-------------------------------------------------------------------------------
set(RBVMS_EXE rbvms)
set(RBVMS_MAIN rbvms.cpp)
set(RBVMS_SOURCES monitor.cpp tau.cpp weakform.cpp)
set(RBVMS_HEADERS monitor.hpp tau.hpp weakform.hpp)
set(LIBRARIES_LIST mfem)

#-------------------------------------------------------------------------------
# If CUDA is enabled, tag source files to be compiled with nvcc.
#-------------------------------------------------------------------------------
if (MFEM_USE_CUDA)
  set_source_files_properties(${RBVMS_MAIN} PROPERTIES LANGUAGE CUDA)
  if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.12.0)
    list(TRANSFORM EXTRA_OPTIONS_LIST PREPEND "-Xcompiler=")
  else()
    set(LIST_)
    foreach(item IN LISTS EXTRA_OPTIONS_LIST)
      list(APPEND LIST_ "-Xcompiler=${item}")
    endforeach()
    set(EXTRA_OPTIONS_LIST ${LIST_})
  endif()
endif()

#-------------------------------------------------------------------------------
# Actually add the executable
#-------------------------------------------------------------------------------
add_executable(${RBVMS_EXE} ${RBVMS_MAIN}
               ${RBVMS_SOURCES}
               ${RBVMS_HEADERS})
if (MFEM_USE_CUDA)
  set_target_properties(${RBVMS_EXE} PROPERTIES
    CUDA_RESOLVE_DEVICE_SYMBOLS ON)
endif()

#-------------------------------------------------------------------------------
# Append the additional libraries and options
#-------------------------------------------------------------------------------
if (LIBRARIES_LIST)
  target_link_libraries(${RBVMS_EXE} PRIVATE ${LIBRARIES_LIST})
endif()
if (EXTRA_OPTIONS_LIST)
  string(REPLACE ";" " " EXTRA_OPTIONS_STRING "${EXTRA_OPTIONS_LIST}")
  message(STATUS "${RBVMS_EXE}: add flags \"${EXTRA_OPTIONS_STRING}\"")
  target_compile_options(${RBVMS_EXE} PRIVATE ${EXTRA_OPTIONS_LIST})
endif()
if (EXTRA_DEFINES_LIST)
  target_compile_definitions(${RBVMS_EXE} PRIVATE ${EXTRA_DEFINES_LIST})
endif()

#-------------------------------------------------------------------------------
# Add tests
#-------------------------------------------------------------------------------
#if (MFEM_ENABLE_TESTING)
#  add_test(NAME ex_condif
#    COMMAND $<TARGET_FILE:condif> -no-vis
#    -m ${PROJECT_SOURCE_DIR}/data/square-nurbs.mesh -r 1 -o 2)
#endif()

#if (MFEM_ENABLE_TESTING)
#    add_test(NAME navsto_np=4
#      COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MFEM_MPI_NP}
#      ${MPIEXEC_PREFLAGS} $<TARGET_FILE:navsto> -no-vis
#      ${MPIEXEC_POSTFLAGS}
#      -m ${PROJECT_SOURCE_DIR}/data/square-nurbs.mesh -r 1 -o 2)
#endif()
