cmake_minimum_required(VERSION 3.18)
project(my_proj LANGUAGES C Fortran) # C is needed to discover Homebrew OpenBLAS

# Find BLAS and LAPACK
include(CMakeFindDependencyMacro)
find_dependency(BLAS)
find_dependency(LAPACK)
find_package(BLAS)
find_package(LAPACK)

# Find OpenMP
if(WITH_OpenMP)
  find_package(OpenMP REQUIRED)
endif()

# Find X11. Can be added for GUI, if DISLIN is also available
 if(WITH_GUI)
  find_package(X11 REQUIRED)
 endif()

# Turn on preprocessor
set(CMAKE_Fortran_PREPROCESS ON)

# Add source files
file (GLOB srcs CONFIGURE_DEPENDS
    ${PROJECT_SOURCE_DIR}/*.f90
    ${PROJECT_SOURCE_DIR}/*.F
    ${PROJECT_SOURCE_DIR}/noGUI/*.f90
    ${PROJECT_SOURCE_DIR}/libreta_hybrid/*.f90
)

# Add executable
add_executable(multiwfn ${srcs})
target_link_libraries(
  multiwfn
  PRIVATE
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  $<$<BOOL:${OpenMP_Fortran_FOUND}>:OpenMP::OpenMP_Fortran>
  $<$<BOOL:${X11_FOUND}>:${X11_LIBRARIES}>
)
target_compile_options(multiwfn PRIVATE -std=legacy -ffree-line-length-none)

# Install
install(TARGETS multiwfn)
install(FILES settings.ini DESTINATION ${CMAKE_INSTALL_PREFIX}/etc/multiwfn)
