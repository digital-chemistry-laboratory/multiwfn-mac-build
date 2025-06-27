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

# Find flint
find_package(PkgConfig)
pkg_check_modules(GMP REQUIRED gmp)
pkg_check_modules(FLINT REQUIRED flint)
# find arb for flint < 3
if(FLINT_VERSION VERSION_LESS "3")
  pkg_check_modules(ARB REQUIRED arb)
endif()
# flint has a bug with PkgConfig
if(FLINT_VERSION VERSION_LESS "3.0.2")
  find_library(FLINT_BUG_LIB NAMES flint REQUIRED)
endif()

# Find OpenMP
if(WITH_OpenMP)
  find_package(OpenMP REQUIRED)
endif()

# Find graphics libraries
if(WITH_GUI)
  find_library(
    dislin
    NAMES dislin_d
    PATHS ${DISLIN_DIR} REQUIRED)
  find_library(Xm NAMES Xm REQUIRED)
  find_library(X11 NAMES X11 REQUIRED)
  find_library(Xt NAMES Xt REQUIRED)
  find_library(GL NAMES GL REQUIRED)
  message("-- Found DISLIN: ${dislin}")
  message("-- Found X11: ${X11}")
  message("-- Found Xm: ${Xm}")
  message("-- Found Xt: ${Xt}")
  message("-- Found GL: ${GL}")
endif()

# Turn on preprocessor
set(CMAKE_Fortran_PREPROCESS ON)

# Add source files
file(
  GLOB
  srcs
  CONFIGURE_DEPENDS
  ${PROJECT_SOURCE_DIR}/*.f90
  ${PROJECT_SOURCE_DIR}/*.F
  ${PROJECT_SOURCE_DIR}/ext/2F2.f90
  ${PROJECT_SOURCE_DIR}/ext/2F2.c
  ${PROJECT_SOURCE_DIR}/ext/frj.f90
  ${PROJECT_SOURCE_DIR}/libreta_hybrid/*.f90)

if(WITH_GUI)
  list(APPEND srcs ${PROJECT_SOURCE_DIR}/ext/xlib.f90)
else()
  list(REMOVE_ITEM srcs ${PROJECT_SOURCE_DIR}/mouse_rotate.f90)
  list(APPEND srcs ${PROJECT_SOURCE_DIR}/noGUI/mouse_rotate_empty.f90
       ${PROJECT_SOURCE_DIR}/noGUI/dislin_d_empty.f90)
endif()

# Add executable
add_executable(multiwfn ${srcs})
target_link_libraries(
  multiwfn
  PRIVATE ${BLAS_LIBRARIES}
          ${LAPACK_LIBRARIES}
          ${FLINT_BUG_LIB}
          ${FLINT_LINK_LIBRARIES}
          ${ARB_LINK_LIBRARIES}
          ${GMP_LINK_LIBRARIES}
          $<$<BOOL:${OpenMP_Fortran_FOUND}>:OpenMP::OpenMP_Fortran>
          $<$<BOOL:${WITH_GUI}>:${dislin}>
          $<$<BOOL:${WITH_GUI}>:${X11}>
          $<$<BOOL:${WITH_GUI}>:${Xm}>
          $<$<BOOL:${WITH_GUI}>:${Xt}>
          $<$<BOOL:${WITH_GUI}>:${GL}>)
target_include_directories(
  multiwfn PRIVATE ${GMP_INCLUDE_DIRS} ${FLINT_INCLUDE_DIRS}
                   ${ARB_INCLUDE_DIRS})

# Add compiler flags to GFortran
target_compile_options(
  multiwfn
  PRIVATE
    $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<STREQUAL:$<Fortran_COMPILER_ID>,GNU>>:-ffree-line-length-none>
)

# Install
install(TARGETS multiwfn)
install(FILES settings.ini DESTINATION ${CMAKE_INSTALL_PREFIX}/etc/multiwfn)
