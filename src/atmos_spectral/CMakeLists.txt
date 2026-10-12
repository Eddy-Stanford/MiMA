# The spectral dynamical core and its coupling to the physics:
#   atmosphere/   the atmosphere component the coupler calls, and the physics interface
#   dynamics/     the spectral time step, damping, corrections and tracer advection
#   diagnostics/  the per-step dynamics diagnostics (dynamics_every)
#   init/         the vertical coordinate, topography and the initial state
#   transforms/   the spectral and Fourier transforms, the FFT and the m/latitude decomposition
set ( ATMOS_SPECTRAL_SOURCES
atmosphere/spectral_physics.f90
atmosphere/atmosphere.f90
dynamics/leapfrog.f90
dynamics/global_integral.f90
dynamics/water_borrowing.f90
dynamics/tracer_type.f90
dynamics/fv_advection.f90
dynamics/press_and_geopot.f90
dynamics/implicit.f90
dynamics/spectral_dynamics.f90
dynamics/spectral_damping.f90
dynamics/matrix_invert.f90
diagnostics/every_step_diagnostics.f90
init/vert_coordinate.f90
init/spectral_initialize_fields.f90
init/spectral_init_cond.f90
init/topog_regularization.f90
transforms/spec_mpp.f90
transforms/gauss_and_legendre.f90
transforms/transforms.f90
transforms/spherical.f90
transforms/spherical_fourier.f90
transforms/grid_fourier.f90
transforms/fft99.f90
transforms/fft.f90
)

target_sources( mima PRIVATE ${ATMOS_SPECTRAL_SOURCES} )
