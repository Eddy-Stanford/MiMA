set ( ATMOS_SHARED_SOURCES
tracer_driver/atmos_tracer_driver.f90
tracer_driver/atmos_radon.f90
tracer_driver/atmos_tracer_utilities.f90
interpolator/interpolator.F90
vert_advection/vert_advection.f90
)

target_sources ( mima PRIVATE ${ATMOS_SHARED_SOURCES} )
