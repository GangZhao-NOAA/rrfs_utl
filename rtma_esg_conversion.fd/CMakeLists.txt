list(APPEND src_rll2esg
  mod_rtma_regrid.f90
  rtma_regrid_rll2esg.f90)

add_executable(rtma_regrid_rll2esg.exe ${src_rll2esg})
target_link_libraries(rtma_regrid_rll2esg.exe PUBLIC ${PROJECT_NAME}::pesglib2
                                              PUBLIC bacio::bacio_4
                                              PUBLIC g2::g2_4
                                              PUBLIC g2tmpl::g2tmpl
                                              PUBLIC ip::ip_d
                                              PUBLIC sp::sp_d
                                              PUBLIC NetCDF::NetCDF_Fortran
                                              PUBLIC MPI::MPI_Fortran)
if(OpenMP_Fortran_FOUND)
  target_link_libraries(rtma_regrid_rll2esg.exe PRIVATE OpenMP::OpenMP_Fortran)
endif()
list(APPEND ESG_CONVERSION_Targets rtma_regrid_rll2esg.exe)

list(APPEND src_esg2rll
  mod_rtma_regrid.f90
  rtma_regrid_esg2rll.f90)

add_executable(rtma_regrid_esg2rll.exe ${src_esg2rll})
target_link_libraries(rtma_regrid_esg2rll.exe PUBLIC ${PROJECT_NAME}::pesglib2
                                              PUBLIC bacio::bacio_4
                                              PUBLIC g2::g2_4
                                              PUBLIC g2tmpl::g2tmpl
                                              PUBLIC ip::ip_d
                                              PUBLIC sp::sp_d
                                              PUBLIC NetCDF::NetCDF_Fortran
                                              PUBLIC MPI::MPI_Fortran)
if(OpenMP_Fortran_FOUND)
  target_link_libraries(rtma_regrid_esg2rll.exe PRIVATE OpenMP::OpenMP_Fortran)
endif()
list(APPEND ESG_CONVERSION_Targets rtma_regrid_esg2rll.exe)


install(
  TARGETS ${ESG_CONVERSION_Targets}
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
