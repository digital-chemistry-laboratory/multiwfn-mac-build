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

# Append source files
set(LIBRETAPATH "./libreta_hybrid")
set(srcs)
list(APPEND srcs
    "dislin_d.f90"
    "define.f90"
    "Bspline.f90"
    "util.f90"
    "function.f90"
    "plot.f90"
    "GUI.f90"
    "DFTxclib.F"
    "Lebedev-Laikov.F"
    "sym.F"
    "edflib.f90"
    "atmraddens.f90"
    "minpack.f90"
    "fparser.f90"
    "sub.f90"
    "integral.f90"
    "fileIO.f90"
    "spectrum.f90"
    "DOS.f90"
    "Multiwfn.f90"
    "0123dim.f90"
    "LSB.f90"
    "population.f90"
    "orbcomp.f90"
    "bondorder.f90"
    "topology.f90"
    "excittrans.f90"
    "otherfunc.f90"
    "otherfunc2.f90"
    "otherfunc3.f90"
    "O1.f90"
    "surfana.f90"
    "procgriddata.f90"
    "AdNDP.f90"
    "fuzzy.f90"
    "CDA.f90"
    "basin.f90"
    "orbloc.f90"
    "visweak.f90"
    "EDA.f90"
    "CDFT.f90"
    "ETS_NOCV.f90"
    "NAONBO.f90"
    "grid.f90"
    "PBC.f90"
    "hyper_polar.f90"
    "deloc_aromat.f90"
    "noGUI/dislin_d_empty.f90"
    "${LIBRETAPATH}/libreta.f90"
    "${LIBRETAPATH}/hrr_012345.f90"
    "${LIBRETAPATH}/blockhrr_012345.f90"
    "${LIBRETAPATH}/ean.f90"
    "${LIBRETAPATH}/eanvrr_012345.f90"
    "${LIBRETAPATH}/boysfunc.f90"
    "${LIBRETAPATH}/naiveeri.f90"
    "${LIBRETAPATH}/ryspoly.f90"
)

# Add executable
add_executable(multiwfn ${srcs})
target_link_libraries(multiwfn PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_compile_options(multiwfn PRIVATE -std=legacy -ffree-line-length-none)