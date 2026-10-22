# MiMA's own code that is used by several components: the restart file layout and the
# check of input.nml for removed namelist groups and variables
set ( MIMA_SHARED_SOURCES
restart/restart_file.f90
namelist_check.f90
)

target_sources ( mima PRIVATE ${MIMA_SHARED_SOURCES} )
