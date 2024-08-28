#-------------------------------------------------------------------------------
# Add tests
#-------------------------------------------------------------------------------
#add_test(NAME navsto_np=4
#         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${RBVMS_MPI_NP}
#         ${MPIEXEC_PREFLAGS} $<TARGET_FILE:navsto> -no-vis
#         ${MPIEXEC_POSTFLAGS}
#        -m ${PROJECT_SOURCE_DIR}/data/square-nurbs.mesh -r 1 -o 2)

