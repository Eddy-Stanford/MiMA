# Provide the FMS library (https://github.com/NOAA-GFDL/FMS) with 8-byte reals,
# and set MIMA_FMS_TARGET to its CMake target.
#
# 1. An installed FMS 2026.01.01 or newer with 8-byte reals is used if CMake
#    finds it: set CMAKE_PREFIX_PATH or FMS_ROOT to its prefix. Both the 64BIT
#    build (FMS::fms_r8) and the default build (FMS::fms; Spack's
#    precision=mixed) have 8-byte reals. FMS 2026.01 is not supported: it
#    crashes when writing restarts on more than one PE.
# 2. Otherwise FMS 2026.02 is downloaded and built with MiMA. For an offline
#    build, set FETCHCONTENT_SOURCE_DIR_FMS to an unpacked FMS 2026.02 source tree
#    (older FMS releases cannot be built this way).
#
# FMS keeps its own compiler flags, but any -ffp-contract flag given for
# MiMA's build (CMAKE_<LANG>_FLAGS or CMAKE_<LANG>_FLAGS_<CONFIG>) is passed
# on to it, so that fused multiply-add contraction is the same in both.

set(MIMA_FMS_MIN_VERSION 2026.01.01)
set(MIMA_FMS_VERSION 2026.02)
set(MIMA_FMS_SHA256 65db44c961089c5e004dd8774cc4cfee75373c4684590d1146c8ff971f8480b7)

find_package(FMS ${MIMA_FMS_MIN_VERSION} CONFIG QUIET)

if(FMS_FOUND)
  foreach(_target FMS::fms_r8 FMS::fms)
    if(TARGET ${_target})
      message(STATUS "Using installed FMS ${FMS_VERSION} (${_target}): ${FMS_DIR}")
      set(MIMA_FMS_TARGET ${_target})
      return()
    endif()
  endforeach()
  message(WARNING "The FMS at ${FMS_DIR} has no 8-byte real library (it was "
                  "built with 32BIT only); downloading FMS ${MIMA_FMS_VERSION} instead.")
elseif(FMS_CONSIDERED_VERSIONS)
  message(WARNING "Ignoring FMS ${FMS_CONSIDERED_VERSIONS} (${FMS_CONSIDERED_CONFIGS}): "
                  "MiMA needs FMS ${MIMA_FMS_MIN_VERSION} or newer; downloading FMS ${MIMA_FMS_VERSION} instead.")
endif()

include(FetchContent)
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)  # extracted files get the time of extraction
endif()

function(mima_add_fms)
  # FMS build options (plain variables here, so they apply to FMS only)
  set(64BIT       ON)
  set(32BIT       OFF)
  set(OPENMP      ${MIMA_OPENMP})
  set(OPENACC     OFF)
  set(WITH_YAML   OFF)
  set(CONSTANTS   GFDL)
  set(UNIT_TESTS  OFF)
  set(SHARED_LIBS OFF)
  set(BUILD_TESTING OFF)

  # FMS replaces CMAKE_<LANG>_FLAGS_<CONFIG> with its own flags but appends
  # to CMAKE_<LANG>_FLAGS, so pass any -ffp-contract flags on through the latter.
  string(TOUPPER "${CMAKE_BUILD_TYPE}" _config)
  foreach(_lang Fortran C)
    string(REGEX MATCHALL "-ffp-contract=[a-z]+" _fp_contract
           "${CMAKE_${_lang}_FLAGS} ${CMAKE_${_lang}_FLAGS_${_config}}")
    if(_fp_contract)
      list(GET _fp_contract -1 _fp_contract)
      set(CMAKE_${_lang}_FLAGS "${CMAKE_${_lang}_FLAGS} ${_fp_contract}")
    endif()
  endforeach()

  # FMS is added with EXCLUDE_FROM_ALL, so that it is built only as a
  # dependency of MiMA and 'cmake --install' installs MiMA alone.
  # FetchContent_Declare accepts EXCLUDE_FROM_ALL from CMake 3.28; older
  # versions populate FMS and add it by hand.
  set(_fetch_args)
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.28)
    list(APPEND _fetch_args EXCLUDE_FROM_ALL)
  endif()
  FetchContent_Declare(FMS
    URL      https://github.com/NOAA-GFDL/FMS/archive/refs/tags/${MIMA_FMS_VERSION}.tar.gz
    URL_HASH SHA256=${MIMA_FMS_SHA256}
    ${_fetch_args})
  if(FETCHCONTENT_SOURCE_DIR_FMS)
    message(STATUS "Building FMS from ${FETCHCONTENT_SOURCE_DIR_FMS}")
  else()
    message(STATUS "FMS not found; downloading and building FMS ${MIMA_FMS_VERSION}")
  endif()
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.28)
    FetchContent_MakeAvailable(FMS)
  else()
    FetchContent_GetProperties(FMS)
    if(NOT fms_POPULATED)
      FetchContent_Populate(FMS)
      add_subdirectory(${fms_SOURCE_DIR} ${fms_BINARY_DIR} EXCLUDE_FROM_ALL)
    endif()
  endif()
endfunction()

mima_add_fms()
set(MIMA_FMS_TARGET FMS::fms_r8)
