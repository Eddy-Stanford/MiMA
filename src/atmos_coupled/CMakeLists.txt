set ( ATMOS_COUPLED_SOURCES
atmos_model.f90
)

target_sources( mima PRIVATE ${ATMOS_COUPLED_SOURCES} )
