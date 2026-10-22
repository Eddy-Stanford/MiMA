set ( COUPLER_SOURCES
coupler_main.f90
)

# NB this is where the executable is defined.
target_sources ( mima PRIVATE ${COUPLER_SOURCES} )
