# MiMA's own code that is used by several components: the restart file layout
set ( MIMA_SHARED_SOURCES
restart/restart_file.f90
)

target_sources ( mima PRIVATE ${MIMA_SHARED_SOURCES} )
