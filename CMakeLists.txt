cmake_minimum_required(VERSION 3.22)

# The compilers are chosen when a build directory is first configured, from
# the FC and CC environment variables (e.g. FC=ifx CC=icx cmake -S . -B build)
# or from -DCMAKE_Fortran_COMPILER=... -DCMAKE_C_COMPILER=... To change them,
# configure a new build directory.
project(MiMA VERSION 2.0.0 LANGUAGES Fortran C)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "" FORCE)
endif()

# The model executable, mima; the subdirectories add its sources
add_executable ( mima )

option(INSTALL_EXEC "Install executables to exec/ directory?" OFF)

# OpenMP threads some of MiMA's loops. Off by default, since an unset
# OMP_NUM_THREADS would oversubscribe the cores used by MPI. A downloaded FMS
# is built with OpenMP only if this is on and the C compiler supports it
# (AppleClang does not).
option(MIMA_OPENMP "Build with OpenMP" OFF)

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
    target_compile_options( mima PRIVATE -cpp -fallow-argument-mismatch -fbacktrace -fdefault-real-8 -fdefault-double-8 -fno-range-check -ffree-line-length-none )
    if(NOT MIMA_OPENMP)
        # -fopenmp implies -frecursive (local arrays on the stack); keep that without OpenMP
        target_compile_options( mima PRIVATE -frecursive )
    endif()
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -march=native ${CMAKE_Fortran_FLAGS_RELEASE}")
endif()

# Add various subdirectories with long lists of source files
add_subdirectory( src/coupler )
add_subdirectory( src/surface )
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
    set(MIMA_RUN_CONFIG "input" CACHE STRING
        "Directory (relative to the source tree) whose input.nml, diag_table and field_table INSTALL_EXEC copies to exec/, e.g. examples/held_suarez")
    get_filename_component(_mima_run_config "${MIMA_RUN_CONFIG}" ABSOLUTE BASE_DIR "${CMAKE_SOURCE_DIR}")
    foreach(_f input.nml diag_table field_table)
        if(NOT EXISTS "${_mima_run_config}/${_f}")
            message(FATAL_ERROR "MIMA_RUN_CONFIG: ${_mima_run_config}/${_f} does not exist")
        endif()
    endforeach()
    message(STATUS "cmake --install creates a run directory in ${CMAKE_INSTALL_PREFIX} from ${MIMA_RUN_CONFIG}/")
    install(TARGETS mima RUNTIME DESTINATION .)
    install(FILES "${_mima_run_config}/input.nml" "${_mima_run_config}/diag_table"
                  "${_mima_run_config}/field_table" DESTINATION .)
    install(DIRECTORY DESTINATION RESTART)
    install(DIRECTORY input/INPUT/ DESTINATION INPUT)
else()
    install(TARGETS mima RUNTIME DESTINATION bin)
endif()
