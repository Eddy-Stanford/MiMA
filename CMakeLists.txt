cmake_minimum_required(VERSION 3.22)
#policy CMP0076 - target_sources source files are relative to file where target_sources is run
cmake_policy (SET CMP0076 NEW)
project(MiMA VERSION 2.0.0 LANGUAGES Fortran C)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "" FORCE)
endif()

# Generate the main mima.x executable with dirs, libs, and opts
add_executable ( mima )
target_compile_definitions( mima PRIVATE use_libMPI use_netCDF gFortran ) # gFortran appears to be unused

option(INSTALL_EXEC "Install executables to exec/ directory?" OFF)

# MiMA's own code has no OpenMP directives. A downloaded FMS is built with
# OpenMP only if this is on; that needs OpenMP for C too, which AppleClang lacks.
if(APPLE)
    set(_mima_openmp_default OFF)
else()
    set(_mima_openmp_default ON)
endif()
option(MIMA_OPENMP "Build with OpenMP" ${_mima_openmp_default})

#Add cmake directory to the environment module variable
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package (MPI REQUIRED)
find_package(NetCDF REQUIRED COMPONENTS C Fortran)

# Link MPI, OpenMP, and netCDF libraries
target_link_libraries( mima PRIVATE MPI::MPI_Fortran )
if(MIMA_OPENMP)
    find_package (OpenMP REQUIRED COMPONENTS Fortran)
    target_link_libraries( mima PRIVATE OpenMP::OpenMP_Fortran )
endif()
target_link_libraries(mima PRIVATE NetCDF::NetCDF_Fortran NetCDF::NetCDF_C)

if(NOT CMAKE_Fortran_COMPILER)
    set(CMAKE_Fortran_COMPILER "gfortran")
endif()


if(APPLE)
    message(STATUS "Configuring for MacOS")
    if(NOT CMAKE_C_COMPILER)
        set(CMAKE_C_COMPILER "clang")
    endif()
    set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE STRING "Build architectures for MacOS")
endif()

if(NOT CMAKE_C_COMPILER)
    set(CMAKE_C_COMPILER "gcc")
endif()

# The install prefix must be set before FMS is configured, which reads it.
if(INSTALL_EXEC)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/exec" CACHE PATH "..." FORCE)
endif()

# FMS, found or downloaded (see cmake/MiMAFMS.cmake). This comes before
# MiMA's compiler flags are set, so that they are not applied to FMS.
include(MiMAFMS)
target_link_libraries( mima PRIVATE ${MIMA_FMS_TARGET} )

# Compiler flags for MiMA's own code. The FMS library is built separately
# with its own flags (also with 8-byte reals).
if(CMAKE_Fortran_COMPILER_ID  MATCHES "Intel|IntelLLVM")
    message(STATUS "Configuring for Intel Fortran")
    target_compile_options( mima PRIVATE -i4 -r8 -fpp -fPIC )
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -xHost ${CMAKE_Fortran_FLAGS_RELEASE} ")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    message(STATUS "Configuring for GNU Fortran")
    target_compile_options( mima PRIVATE -cpp -fallow-argument-mismatch -fbacktrace -fdefault-real-8 -fdefault-double-8 -fcray-pointer -fallow-invalid-boz -fno-range-check -ffree-line-length-none )
    if(NOT MIMA_OPENMP)
        # -fopenmp implies -frecursive (local arrays on the stack); keep that without OpenMP
        target_compile_options( mima PRIVATE -frecursive )
    endif()
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -march=native ${CMAKE_Fortran_FLAGS_RELEASE}")
endif()

# Add various subdirectories with long lists of source files
add_subdirectory( src/coupler )
add_subdirectory( src/atmos_coupled )
add_subdirectory( src/atmos_param )

set_source_files_properties (
# The following files do nothing but assign very large arrays.
# For some reason when compiling with ifort and optimisation
# the compilation will take a very long time (10s of minutes).
# Since the code doesn't actually *do* anything there's no
# need to waste time having the compiler apply probably
# meaningless optimisation.  So we disable optimisation for these
# files.
src/atmos_param/radiation/rrtm/rrtmg_sw/gcm_model/src/rrtmg_sw_k_g.f90
src/atmos_param/radiation/rrtm/rrtmg_lw/gcm_model/src/rrtmg_lw_k_g.f90
PROPERTIES
COMPILE_FLAGS -O0
)
add_subdirectory( src/atmos_shared )
add_subdirectory( src/atmos_spectral )
add_subdirectory( src/mima_shared )

if(INSTALL_EXEC)
    message(STATUS "Building exec dir instead of installing conventionally")
    install(TARGETS mima RUNTIME DESTINATION .)
    install(FILES input/diag_table DESTINATION . )
    install(FILES input/field_table DESTINATION .)
    install(DIRECTORY DESTINATION RESTART)
    install(FILES input/input.nml DESTINATION . )
    install(DIRECTORY input/INPUT/ DESTINATION INPUT)
else()
    install(TARGETS mima RUNTIME DESTINATION bin)
endif()
