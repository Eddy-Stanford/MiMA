# MiMA's own code that is not specific to one component: parts of the old
# vendored FMS that are no longer in FMS, and the restart file layout
set ( MIMA_SHARED_SOURCES
fft/fft99.f90
fft/fft.F90
restart/restart_file.f90
)

target_sources ( mima PRIVATE ${MIMA_SHARED_SOURCES} )
