set ( SURFACE_SOURCES
surface_types.f90
surface_flux.f90
surface_albedo.f90
slab_tendency.f90
ocean_qflux.f90
ocean_model.f90
land_model.f90
flux_exchange.f90
)

target_sources ( mima PRIVATE ${SURFACE_SOURCES} )
