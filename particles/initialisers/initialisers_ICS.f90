!> module containing initialization subroutines for ICS (impurity coupling scheme)
!> particles, with spatial distribution proportional to background electron density
module initialisers_ICS
  use mod_particle_types
  use mod_particle_sim
  use mod_rng
  use initialisers_base
  use mod_model_settings, only: var_rho
  use mod_coronal
  implicit none

  contains

! Probability density function proportional to background electron density
! var(1) = var_rho (normalized plasma density, n_e = central_density * rho * 1e20)
pure function density_pdf(var) result(p)
  real*8, intent(in) :: var(:)
  real*8             :: p

  p = var(1)
end function density_pdf

!> Initialize ICS particles with spatial distribution proportional to
!> background electron density (n_e) via rejection sampling.
!> Velocity is set from local Maxwellian with vpar, charge state from coronal equilibrium.
subroutine ics_density_initialization(sim, group_num, concentration)
  use mod_pcg32_rng

  type(particle_sim), intent(inout) :: sim
  integer,            intent(in)    :: group_num
  real*8,             intent(in)    :: concentration

  integer :: num_part

  ! Build coronal equilibrium table from ADAS data for charge state determination
  sim%groups(group_num)%cor = coronal(sim%groups(group_num)%ad)

  ! 1. Spatial initialization: rejection sampling with PDF proportional to n_e
  call initialise_particles(sim%groups(group_num)%particles, &
      sim%fields%node_list, sim%fields%element_list, pcg32_rng(), &
      variables=[var_rho], transform=density_pdf)

  num_part = size(sim%groups(group_num)%particles, 1)

  ! 2. Velocity initialization: Maxwellian at local Ti=Te=T/2, with parallel flow and coronal charge state
  call set_velocity_from_T(sim%groups(group_num)%particles, &
      sim%groups(group_num)%mass, &
      sim%fields%node_list, sim%fields%element_list, &
      sim%groups(group_num)%cor, v_par=.true.)

  ! 3. Set particle weight = concentration
  !    (spatial distribution already reflects n_e proportionality via rejection sampling)
  sim%groups(group_num)%particles(:)%weight = concentration

end subroutine ics_density_initialization

end module initialisers_ICS
