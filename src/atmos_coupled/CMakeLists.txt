set ( ATMOS_COUPLED_SOURCES
atmos_model.f90
#atmos_model.html
)

target_sources( mima PRIVATE ${ATMOS_COUPLED_SOURCES} )
