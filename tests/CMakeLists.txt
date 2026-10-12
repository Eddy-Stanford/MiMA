# MiMA tests, added by -DMIMA_TESTING=ON. Run them with: ctest --test-dir <build>
# What each test proves, and the disk and CPU rules: tests/README.md.

cmake_host_system_information(RESULT ncores QUERY NUMBER_OF_PHYSICAL_CORES)
set(MIMA_TEST_NPES 4 CACHE STRING "MPI ranks per model test (capped at the physical cores)")
set(MIMA_TEST_REFERENCE "" CACHE PATH "Directory of .digest files of a reference build (adds the compare_* tests)")
set(npes ${MIMA_TEST_NPES})
if(npes GREATER ncores)
  set(npes ${ncores})
endif()

set(run ${CMAKE_CURRENT_SOURCE_DIR}/run_case.sh)
set(exe $<TARGET_FILE:mima>)
set(out ${CMAKE_CURRENT_BINARY_DIR})   # <build>/tests: run/ while a test runs, digests/ after
set(digests ${out}/digests)

# Every test runs the model, and run_case.sh allows one run at a time (a lock
# file), so all tests share one RESOURCE_LOCK: `ctest -j` never overlaps them.
function(mima_test name label nproc)
  add_test(NAME ${name} COMMAND ${run} ${ARGN})
  set_tests_properties(${name} PROPERTIES LABELS ${label} PROCESSORS ${nproc}
    RESOURCE_LOCK mima_run TIMEOUT 1800 ENVIRONMENT "MIMA_MPIRUN=${MPIEXEC_EXECUTABLE}")
endfunction()

# smoke: the shipped configurations run and write no NaN
foreach(case rrtm gray hs)
  mima_test(smoke_${case} smoke ${npes} run ${exe} ${case} ${out} ${npes})
  set_tests_properties(smoke_${case} PROPERTIES FIXTURES_SETUP ${case}_digest)
endforeach()

# self: restart = continuous (20 steps + 20 steps, not a multiple of the 9-step
# radiation interval), and the code defaults = input/input.nml
mima_test(restart_rrtm self ${npes} restart ${exe} rrtm ${out} ${npes} 10000,10000)
mima_test(defaults_rrtm self ${npes} run ${exe} rrtm_defaults ${out} ${npes} ${digests}/rrtm.digest)
set_tests_properties(defaults_rrtm PROPERTIES FIXTURES_REQUIRED rrtm_digest)

# guard: bad input stops with the expected message (one line of guards.txt each)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS guards.txt)
file(STRINGS guards.txt guards REGEX "^[a-z]")
foreach(line IN LISTS guards)
  if(NOT line MATCHES "^([^ ]+) +([^ ]+) +([^ ]+) +(.+)$")
    message(FATAL_ERROR "tests/guards.txt: cannot read the line '${line}'")
  endif()
  set(regex "${CMAKE_MATCH_4}")
  mima_test(guard_${CMAKE_MATCH_1} guard 1 guard ${exe} ${CMAKE_MATCH_2} ${out} ${CMAKE_MATCH_1} ${CMAKE_MATCH_3})
  set_tests_properties(guard_${CMAKE_MATCH_1} PROPERTIES PASS_REGULAR_EXPRESSION "${regex}")
endforeach()

# ref: bit for bit against the digests of another build (see README.md)
if(MIMA_TEST_REFERENCE)
  foreach(case rrtm gray hs)
    add_test(NAME compare_${case} COMMAND ${run} compare ${MIMA_TEST_REFERENCE}/${case}.digest ${digests}/${case}.digest)
    set_tests_properties(compare_${case} PROPERTIES LABELS ref FIXTURES_REQUIRED ${case}_digest)
  endforeach()
endif()
