list(APPEND src_rll2esg
regrid_rll2esg_geoll.f90)

add_executable(regrid_rll2esg_geoll.exe ${src_rll2esg})
target_link_libraries(regrid_rll2esg_geoll.exe PUBLIC bacio::bacio_4)
target_link_libraries(regrid_rll2esg_geoll.exe PUBLIC ncio::ncio)
target_link_libraries(regrid_rll2esg_geoll.exe PUBLIC g2::g2_4)
target_link_libraries(regrid_rll2esg_geoll.exe PUBLIC ip::ip_4)
target_link_libraries(regrid_rll2esg_geoll.exe PUBLIC sp::sp_d)
target_link_libraries(regrid_rll2esg_geoll.exe PUBLIC NetCDF::NetCDF_Fortran)
target_link_libraries(regrid_rll2esg_geoll.exe PUBLIC pesglib2::pesglib2)
target_link_libraries(regrid_rll2esg_geoll.exe PUBLIC MPI::MPI_Fortran)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(regrid_rll2esg_geoll.exe PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(
  TARGETS regrid_rll2esg_geoll.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
