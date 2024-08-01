list(APPEND src_rll2esg
  mod_rtma_regrid.F90
  rtma_regrid_rll2esg.F90)

add_executable(rtma_regrid_rll2esg.exe ${src_rll2esg})
target_link_libraries(rtma_regrid_rll2esg.exe PRIVATE ${PROJECT_NAME}::pesglib2
                                              PRIVATE bacio::bacio_4
                                              PRIVATE g2::g2_4
                                              PRIVATE g2tmpl::g2tmpl
                                              PRIVATE ip::ip_d
                                              PRIVATE sp::sp_d
                                              PRIVATE w3emc::w3emc_4
                                              PRIVATE NetCDF::NetCDF_Fortran
                                              PRIVATE MPI::MPI_Fortran)
if(OpenMP_Fortran_FOUND)
  target_link_libraries(rtma_regrid_rll2esg.exe PRIVATE OpenMP::OpenMP_Fortran)
endif()
list(APPEND ESG_CONVERSION_Targets rtma_regrid_rll2esg.exe)

list(APPEND src_esg2rll
  mod_rtma_regrid.F90
  rtma_regrid_esg2rll.F90)

add_executable(rtma_regrid_esg2rll.exe ${src_esg2rll})
target_link_libraries(rtma_regrid_esg2rll.exe PRIVATE ${PROJECT_NAME}::pesglib2
                                              PRIVATE bacio::bacio_4
                                              PRIVATE g2::g2_4
                                              PRIVATE g2tmpl::g2tmpl
                                              PRIVATE ip::ip_d
                                              PRIVATE sp::sp_d
                                              PRIVATE w3emc::w3emc_4
                                              PRIVATE NetCDF::NetCDF_Fortran
                                              PRIVATE MPI::MPI_Fortran)
if(OpenMP_Fortran_FOUND)
  target_link_libraries(rtma_regrid_esg2rll.exe PRIVATE OpenMP::OpenMP_Fortran)
endif()
list(APPEND ESG_CONVERSION_Targets rtma_regrid_esg2rll.exe)

# if(ip_VERSION GREATER_EQUAL 4.0.0)
#   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DIP_V4")
# endif()
if(ip_VERSION LESS 4.0.0)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DIP_V3")
endif()

install(
  TARGETS ${ESG_CONVERSION_Targets}
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
