list(APPEND src_esg_lib
  pkind.f90
  pietc.f90
  pietc_s.f90
  pmat.f90
  pmat2.f90
  pmat4.f90
  pmat5.f90
  pfun.f90
  psym2.f90
  pesg.f90
  pbswi.f90)

set(module_dir "${CMAKE_CURRENT_BINARY_DIR}/include")
add_library(pesglib2 STATIC ${src_esg_lib})
add_library(${PROJECT_NAME}::pesglib2 ALIAS pesglib2)
# target_link_libraries(pesglib2 PUBLIC ${PROJECT_NAME}::ncio)
set_target_properties(pesglib2 PROPERTIES Fortran_MODULE_DIRECTORY "${module_dir}")
target_include_directories(pesglib2 PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
                                          $<INSTALL_INTERFACE:include>)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include DESTINATION ${CMAKE_INSTALL_PREFIX})

install(
  TARGETS pesglib2
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
