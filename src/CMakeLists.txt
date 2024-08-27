#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
add_executable(rbvms rbvms.cpp
               monitor.cpp tau.cpp weakform.cpp
               monitor.hpp tau.hpp weakform.hpp)

add_dependencies(rbvms mfem)
target_link_libraries(rbvms PRIVATE ${MFEM_LIB})
target_link_libraries(rbvms PRIVATE ${HYPRE_LIB})
target_link_libraries(rbvms PRIVATE ${METIS_LIB})

install(TARGETS rbvms RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
#install(PROGRAMS DESTINATION ${INSTALL_BIN_DIR})

# Install the executable
#install(TARGETS delfi-exe RUNTIME DESTINATION bin)


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
