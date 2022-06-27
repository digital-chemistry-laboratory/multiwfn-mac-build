cmake_minimum_required(VERSION 3.18)
project(my_proj LANGUAGES C Fortran) # C is needed to discover Homebrew OpenBLAS

# Options
option(WITH_OpenMP "Build with OpenMP")
option(WITH_GUI "Build with GUI")

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

# Find graphics libraries
if(WITH_GUI)
  find_library(
    dislin
    NAMES dislin_d
    PATHS ${DISLIN_DIR}
    REQUIRED
  )
  find_library(
    Xm
    NAMES Xm
    REQUIRED
  )
  message("-- Found DISLIN: ${dislin}")
  message("-- Found OpenMotif: ${Xm}")  
endif()

# Turn on preprocessor
set(CMAKE_Fortran_PREPROCESS ON)

# Add source files
file(GLOB srcs CONFIGURE_DEPENDS
  ${PROJECT_SOURCE_DIR}/*.f90
  ${PROJECT_SOURCE_DIR}/*.F
  ${PROJECT_SOURCE_DIR}/libreta_hybrid/*.f90
)
if(NOT WITH_GUI)
  list(APPEND srcs ${PROJECT_SOURCE_DIR}/noGUI/dislin_d_empty.f90)
endif()

# Add executable
add_executable(multiwfn ${srcs})
target_link_libraries(
  multiwfn
  PRIVATE
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  $<$<BOOL:${OpenMP_Fortran_FOUND}>:OpenMP::OpenMP_Fortran>
  $<$<BOOL:${WITH_GUI}>:${dislin}>
  $<$<BOOL:${WITH_GUI}>:${Xm}>
)
target_compile_options(multiwfn PRIVATE -ffree-line-length-none)

# Install
install(TARGETS multiwfn)
install(FILES settings.ini DESTINATION ${CMAKE_INSTALL_PREFIX}/etc/multiwfn)