list(APPEND src_rll2esg
regrid_rll2esg_geoll.f90)
add_executable(regrid_rll2esg_geoll.exe ${src_rll2esg})
target_link_libraries(regrid_rll2esg_geoll.exe PUBLIC ${PROJECT_NAME}::pesglib2
                                               PUBLIC bacio::bacio_4
                                               PUBLIC g2::g2_4
                                               PUBLIC g2tmpl::g2tmpl
                                               PUBLIC ip::ip_4
                                               PUBLIC sp::sp_d
                                               PUBLIC NetCDF::NetCDF_Fortran
                                               PUBLIC MPI::MPI_Fortran)
if(OpenMP_Fortran_FOUND)
  target_link_libraries(regrid_rll2esg_geoll.exe PRIVATE OpenMP::OpenMP_Fortran)
endif()
list(APPEND ESG_CONVERSION_Targets regrid_rll2esg_geoll.exe)

list(APPEND src_esg2rll
regrid_esg2rll_geoll.f90)
add_executable(regrid_esg2rll_geoll.exe ${src_esg2rll})
target_link_libraries(regrid_esg2rll_geoll.exe PUBLIC ${PROJECT_NAME}::pesglib2
                                               PUBLIC bacio::bacio_4
                                               PUBLIC g2::g2_4
                                               PUBLIC g2tmpl::g2tmpl
                                               PUBLIC ip::ip_4
                                               PUBLIC sp::sp_d
                                               PUBLIC NetCDF::NetCDF_Fortran
                                               PUBLIC MPI::MPI_Fortran)
if(OpenMP_Fortran_FOUND)
  target_link_libraries(regrid_esg2rll_geoll.exe PRIVATE OpenMP::OpenMP_Fortran)
endif()
list(APPEND ESG_CONVERSION_Targets regrid_esg2rll_geoll.exe)

install(
  TARGETS ${ESG_CONVERSION_Targets}
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
