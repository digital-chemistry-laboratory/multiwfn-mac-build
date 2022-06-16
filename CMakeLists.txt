cmake_minimum_required(VERSION 3.18)
project(my_proj LANGUAGES C Fortran) # C is needed to discover Homebrew OpenBLAS

# Find BLAS and LAPACK
include(CMakeFindDependencyMacro)
find_dependency(BLAS)
find_dependency(LAPACK)
find_package(BLAS)
find_package(LAPACK)

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
target_link_libraries(multiwfn PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_compile_options(multiwfn PRIVATE -std=legacy -ffree-line-length-none)
install(TARGETS multiwfn)
