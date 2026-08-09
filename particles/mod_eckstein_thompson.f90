!> Mod Eckstein sputter yield and sputtered energy coefficients
!> This module contains types for the needed coefficients for the Eckstein angle
!> dependent fit formulas, for the sputter yield and the sputtered energy coefficient.
module mod_eckstein_thompson
use mod_atomic_elements
use constants
use phys_module, only: central_mass
 
implicit none
private

public :: eckstein_coeff_set, eckstein_sputter_yield, eckstein_sputtered_energy_coeff
public :: chemical_sputtering, ThompsonEne
 
!> Eckstein angle dependency fit coefficients, determining the function
type :: eckstein_coeff
  real*8 :: E0 !< energy [eV]
  real*8 :: f
  real*8 :: b
  real*8 :: c
  real*8 :: Yn !< Sputtering yield/energy coeff at normal incidence at E0
  real*8 :: Esp !< [eV]
  real*8 :: theta_star !< corrected angle of incidence [degrees]
  real*8 :: theta_max !< angle of maximum sputtering yield [degrees]
  contains
    procedure :: eval => evaluate_eckstein_formula
end type eckstein_coeff

!> Base type for a set of eckstein coefficients for sputter yield or sputtered energy coefficient
type, abstract :: eckstein_coeff_set
  integer :: Z_ion !< atomic number of impacting ion
  integer :: Z_target !< atomic number of target ion

  !> Eckstein coeffs for normal incidence
  real*8 :: eps
  real*8 :: q
  real*8 :: mu
  real*8 :: lambda
  real*8 :: E_threshold

  logical :: use_Yn_func = .true. !< By default use the coefficients above to estimate Yn
  
  type(eckstein_coeff), allocatable :: yn(:)
contains
  procedure :: initialise !< initialise coefficients using formulas
  procedure(calc_E), deferred, private :: calc_E !< calculate normal incidence yield/coeff
  procedure, public :: interp !< Call interp_E and interp_theta internally as necessary
end type eckstein_coeff_set

!> Sputter yield calculation
type, extends(eckstein_coeff_set) :: eckstein_sputter_yield
contains
  procedure :: calc_E => calc_E_sputter_yield
end type eckstein_sputter_yield

!> Sputtered energy coefficient calculation
type, extends(eckstein_coeff_set) :: eckstein_sputtered_energy_coeff
  real*8 :: nu
  real*8 :: n
contains
  procedure :: calc_E => calc_E_sputtered_energy_coeff
end type eckstein_sputtered_energy_coeff

interface
  pure function calc_E(this, E) result(yield)
    import eckstein_coeff_set
    class(eckstein_coeff_set), intent(in)  :: this
    real*8, intent(in) :: E
    real*8 :: yield
  end function calc_E
end interface

contains 

!> initialise coefficients using formulas instead of reading from files
pure subroutine initialise(this, Z_ion, Z_target)
  class(eckstein_coeff_set), intent(inout) :: this
  integer, intent(in) :: Z_ion, Z_target
  integer :: i, n_energies
  real*8 :: E_min, E_max, logE_min, logE_max, logE_step
  real*8 :: M_ion, M_target, Es, Eth, ns
  real*8 :: aL, ETF, gamma1, Q
  integer :: Zion_corr
  real*8 :: Z_ion_real,Z_target_real
  if(Z_target.le.1) return   !skip for electrons, hydrogen isotopes (no Eckstein data needed)
  Zion_corr = Z_ion
  if(Z_ion.le.0) Zion_corr = 1
  this%Z_ion = Zion_corr
  this%Z_target = Z_target
  Z_ion_real=real(Zion_corr,kind=8)
  Z_target_real=real(Z_target,kind=8)
  ! Get masses
  M_ion = atomic_weights(Z_ion)
  if (Z_ion == 1) M_ion = central_mass  ! for hydrogen isotopes
  M_target = atomic_weights(Z_target)
  
  ! Get material properties
  Es = get_Es(Z_target)
  Eth = get_Eth(Zion_corr, Z_target)
  ns = get_ns(Z_target)
  
  ! Calculate parameters for normal incidence formula
  aL = 0.4685d0 / sqrt(Z_ion_real**(2.0d0/3.0d0) + Z_target_real**(2.0d0/3.0d0))
  ETF = Z_ion_real * Z_target_real * 14.4d0 / aL * ((M_ion + M_target) / M_target)
  
  ! Calculate Q (sputtering efficiency coefficient)
  gamma1 = 4.0d0 * M_ion * M_target / (M_ion + M_target)**2
  Q = 1.633d0 * Z_ion_real**(2.0d0/3.0d0) * Z_target_real**(2.0d0/3.0d0) * &
      (Z_ion_real**(2.0d0/3.0d0) + Z_target_real**(2.0d0/3.0d0))**(1.0d0/3.0d0) * &
      M_ion**(5.0d0/6.0d0) * M_target**(1.0d0/6.0d0) / (M_ion + M_target) * &
      (0.15d0 + 0.05d0 * (M_target / M_ion)) / &
      (1.0d0 + 0.05d0 * (M_target / M_ion)**1.6d0) * Es**(-2.0d0/3.0d0)
  
  ! Set parameters for normal incidence formula
  this%eps = ETF
  this%q = Q
  this%mu = 1.0d0  ! Default value, can be adjusted based on material
  this%lambda = 0.35d0  ! Default value, can be adjusted based on material
  this%E_threshold = Eth
  

  select type (c => this)
  type is (eckstein_sputtered_energy_coeff)
    c%nu = 1.0d0  ! 默认值
    if (Z_ion > 2) then
      c%n = 2.0d0
    else
      c%n = 1.0d0
    end if
  end select
  
  ! initialise energy points for interpolation
  n_energies = 100
  E_min = Eth * 1.1d0
  E_max = 10000.0d0  ! 10 keV
  
  allocate(this%yn(n_energies))
  
  logE_min = log10(E_min)
  logE_max = log10(E_max)
  logE_step = (logE_max - logE_min) / (n_energies - 1)
  
  ! Calculate coefficients at each energy point
  do i = 1, n_energies
    this%yn(i)%E0 = 10.0d0**(logE_min + (i-1) * logE_step)
    
    ! Calculate normal incidence yield for this energy point
    if (this%yn(i)%E0 > Eth) then
      this%yn(i)%Yn = this%calc_E(this%yn(i)%E0)
    else
      this%yn(i)%Yn = 0.0d0
    end if
    
    ! Calculate Eckstein parameters using formulas
    call calculate_eckstein_params(this%yn(i)%E0, Zion_corr, Z_target, M_ion, M_target, Es, ns, &
                                   this%yn(i)%f, this%yn(i)%b, this%yn(i)%c, &
                                   this%yn(i)%theta_star, this%yn(i)%theta_max)
    
    ! Surface binding energy for sputtered particles
    this%yn(i)%Esp = Es
  end do
end subroutine initialise

!> Calculate Eckstein parameters using formulas
pure subroutine calculate_eckstein_params(E, Z_ion, Z_target, M_ion, M_target, Es, ns, &
                                         f, b, c, theta_star, theta_max)
  real*8, intent(in) :: E
  integer, intent(in) :: Z_ion, Z_target
  real*8, intent(in) :: M_ion, M_target, Es, ns
  real*8, intent(out) :: f, b, c, theta_star, theta_max
  
  real*8 :: aL, ETF, epsilon1, gamma1, aopt
  
  ! Calculate surface roughness parameter f
  f = sqrt(Es) * (0.94d0 - 1.33d-3 * M_target / M_ion)
  
  ! Calculate optimum angle for maximum yield
  aL = 0.4685d0 / sqrt(Z_ion**(2.0d0/3.0d0) + Z_target**(2.0d0/3.0d0))
  ETF = Z_ion * Z_target * 14.4d0 / aL * ((M_ion + M_target) / M_target)
  epsilon1 = E / ETF
  gamma1 = 4.0d0 * M_ion * M_target / (M_ion + M_target)**2
  
  aopt = pi/2.0d0 - aL * ns**(1.0d0/3.0d0) * &
         (2.0d0 * epsilon1 * sqrt(Es/(gamma1*E)))**(-0.5d0)
  
  ! Convert aopt from radians to degrees for theta_star
  theta_star = aopt * 180.0d0 / pi
  
  ! Default values for b and c (can be adjusted based on material)
  b = 2.5d0
  c = 0.8d0
  
  ! Theta_max is typically close to theta_star
  theta_max = theta_star * 0.8d0
end subroutine calculate_eckstein_params

!> Interpolate sputter yield or sputtered energy coefficient from components
pure function interp(this, E, theta) result(yield)
  class(eckstein_coeff_set), intent(in) :: this
  real*8, intent(in)                    :: E !< energy in eV
  real*8, intent(in)                    :: theta !< angle in degrees
  real*8 :: yield
  integer :: i1, i2 !< positions to interpolate
  real*8 :: a !< interpolation factor
  
  if (.not. allocated(this%yn)) then
    yield = 0.d0
    return
  end if

  call interp_factors(E, this%yn, i1, i2, a)

  if (E .le. this%E_threshold) then
    yield = 0.0d0
  else
    if (this%use_Yn_func) then
      yield = this%calc_E(E)
    else ! interpolate log-linear
      yield = this%yn(i1)%Yn**a * this%yn(i2)%Yn**(1.d0-a)
    end if
  end if
  
  ! Apply angle dependence
  yield = yield * (this%yn(i1)%eval(theta)**a * this%yn(i2)%eval(theta)**(1.d0-a))
end function interp

!> Helper subroutine for interpolation factors
pure subroutine interp_factors(E,yn,i1,i2,a)
  real*8, intent(in)   :: E
  type(eckstein_coeff), dimension(:), intent(in) :: yn
  integer, intent(out) :: i1, i2
  real*8, intent(out)  :: a
  integer :: i_delta_E
  real*8 :: E1, E2
  real*8 :: delta_E(size(yn,1))

  delta_E = E - yn(:)%E0
  i_delta_E = minloc(delta_E, dim = 1, mask=delta_E .ge. 0.d0)

  if (i_delta_E .eq. 0) then
    i1 = 1
    i2 = 2
  else if (i_delta_E .eq. size(yn,1)) then
    i1 = i_delta_E-1
    i2 = i_delta_E
  else
    i1 = i_delta_E
    i2 = i_delta_E+1
  end if

  E1 = yn(i1)%E0
  E2 = yn(i2)%E0
  a = (log(E2)-log(E))/(log(E2)-log(E1))
  if (a .lt. 0) then
    a = 0
  end if
end subroutine interp_factors

!> Sputter yield fit formula
pure function calc_E_sputter_yield(this, E) result(yield)
  class(eckstein_sputter_yield), intent(in) :: this
  real*8, intent(in)                        :: E !< energy in eV
  real*8 :: yield
  
  real*8 :: epsilon1, Sn_KRC, w
  
  if (.not. allocated(this%yn)) then
    yield = 0.d0
    return
  end if

  if (E <= this%E_threshold) then
    yield = 0.0d0
    return
  end if

  ! Calculate normalized energy
  epsilon1 = E / this%eps
  
  ! Calculate nuclear stopping power using KrC formula
  w = epsilon1 + 0.1728d0 * sqrt(epsilon1) + 0.008d0 * epsilon1**0.1504d0
  Sn_KRC = 0.5d0 * log(1.0d0 + 1.2288d0 * epsilon1) / w
  
  ! Calculate sputter yield using the simplified formula
  ! yield = Q * Sn_KRC * (1 - (Eth/E)^(2/3)) * (1 - Eth/E)^2
  yield = this%q * Sn_KRC * (1.0d0 - (this%E_threshold/E)**(2.0d0/3.0d0)) * &
          (1.0d0 - this%E_threshold/E)**2
  
  ! Alternative: Use Eckstein original formula
  ! yield = this%q * Sn_KRC * (E/this%E_threshold - 1.0d0)**this%mu / &
  !         (this%lambda + (E/this%E_threshold - 1.0d0)**this%mu)
end function calc_E_sputter_yield

!> Sputtered energy coefficient fit formula
pure function calc_E_sputtered_energy_coeff(this, E) result(yield)
  class(eckstein_sputtered_energy_coeff), intent(in) :: this
  real*8, intent(in)                                 :: E !< energy in eV
  real*8 :: yield
  
  real*8 :: epsilon1, Sn_KRC, w
  
  if (.not. allocated(this%yn)) then
    yield = 0.d0
    return
  end if

  if (E <= this%E_threshold) then
    yield = 0.0d0
    return
  end if

  ! Calculate normalized energy
  epsilon1 = E / this%eps
  
  ! Calculate nuclear stopping power with nu and n parameters
  w = epsilon1 + 0.1728d0 * sqrt(epsilon1) + 0.008d0 * epsilon1**0.1504d0
  Sn_KRC = (0.5d0 * log(1.0d0 + 1.2288d0 * epsilon1) / w**this%nu)**this%n
  
  ! Calculate sputtered energy coefficient using Eckstein formula
  yield = this%q * Sn_KRC * (E/this%E_threshold - 1.0d0)**this%mu / &
          (this%lambda + (E/this%E_threshold - 1.0d0)**this%mu)
end function calc_E_sputtered_energy_coeff

!> Function for evaluating the angle dependency
pure function evaluate_eckstein_formula(this, theta) result(coeff)
  class(eckstein_coeff), intent(in) :: this
  real*8 ,intent(in)                :: theta !< angle in degrees
  real*8                            :: coeff
  
  ! Modified Eckstein formula with surface structure effects
  if (abs(theta) .gt. 1d-6) then
    ! Use surface structure effect formula
    coeff = cos(theta * pi/180.0d0)**(-this%f) * &
            exp(this%f * (1.0d0 - 1.0d0/cos(theta * pi/180.0d0)) * cos(this%theta_star * pi/180.0d0))
  else
    coeff = 1.d0
  end if
end function evaluate_eckstein_formula

!> Chemical sputtering for H-C only
pure function chemical_sputtering(T_target, E_in, flux_in) result(yield)
  real*8, intent(in) :: T_target, E_in, flux_in
  real*8 :: yield
  real*8 :: E_TF_H, Q_H, D_H, E_dam, E_des, E_rel, E_therm
  real*8 :: Q1, D1, E_TF1, ratio, Sn_E0, Ydam, Ydes, flux_judge
  real*8 :: c, c_sp3, Ytherm, Ysurf

  ! constants for H-C chemical sputtering
  E_TF_H = 415.0d0
  Q_H = 0.035d0
  D_H = 250.0d0
  E_dam = 15.0d0
  E_des = 2.0d0
  E_rel = 1.8d0
  E_therm = 1.7d0

  Q1 = Q_H
  D1 = D_H
  E_TF1 = E_TF_H

  ratio = E_in / E_TF1
  Sn_E0 = 0.5d0 * log(1.0d0 + 1.2288d0 * ratio) / (ratio + 0.1728d0 * sqrt(ratio) + 0.008d0 * ratio**0.1504d0)
  Ydam = Q1 * Sn_E0 * (1.0d0 - (E_dam / E_in)**(2.0d0 / 3.0d0)) * (1.0d0 - E_dam / E_in)**2.0d0
  Ydes = Q1 * Sn_E0 * (1.0d0 - (E_des / E_in)**(2.0d0 / 3.0d0)) * (1.0d0 - E_des / E_in)**2.0d0
  flux_judge = 1.0d30 * exp(-1.4d0 / T_target)

  c = 1.0d0 / (1.0d0 + 3.0d-23 * flux_in * merge(1.0d0, 0.0d0, flux_in > flux_judge) + &
               3.0d7 * exp(-1.4d0 / T_target) * merge(1.0d0, 0.0d0, flux_in <= flux_judge))
  
  c_sp3 = c * (2.0d-32 * flux_in + exp(-E_therm / T_target)) / &
          (2.0d-32 * flux_in + (1.0d0 + 2.0d29 / flux_in * exp(-E_rel / T_target)) * exp(-E_therm / T_target))
  
  Ytherm = c_sp3 * 0.033d0 * exp(-E_therm / T_target) / (2.0d-32 * flux_in + exp(-E_therm / T_target))
  Ysurf = c_sp3 * Ydes / (1.0d0 + exp((E_in - 65.0d0) / 40.0d0))
  
  if(Ydam < 0.0d0) Ydam = 0.0d0
  if(Ydes < 0.0d0) Ydes = 0.0d0
  if(Ytherm < 0.0d0) Ytherm = 0.0d0
  if(Ysurf < 0.0d0) Ysurf = 0.0d0
  
  yield = Ytherm * (1.0d0 + D1 * Ydam) + Ysurf
end function chemical_sputtering

!> Get surface binding energy
pure function get_Es(Z_imp) result(Es)
  integer, intent(in) :: Z_imp
  real*8 :: Es  ! surface binding energy (eV)
  
  select case(Z_imp)
  case(74) ! W
    Es = 8.9931d0
  case(6)  ! C
    Es = 7.42d0
  case(42) ! Mo
    Es = 6.83d0
  case(26) ! Fe
    Es = 4.29d0
  case(28) ! Ni
    Es = 4.44d0
  case(29) ! Cu
    Es = 3.49d0
  case(13) ! Al
    Es = 3.36d0
  case(14) ! Si
    Es = 4.70d0
  case default
    ! Use semi-empirical formula for other materials
    Es = 0.5d0 * atomic_weights(Z_imp) / 1000.0d0  ! eV, approximate relation
  end select
end function get_Es

!> Get threshold energy
pure function get_Eth(Z_ion, Z_imp) result(Eth)
  integer, intent(in) :: Z_ion, Z_imp
  real*8 :: Eth   ! threshold energy (eV)
  
  Eth = 0.0d0
  if (Z_ion == 1) then ! H or D or T
    select case(Z_imp)
    case(74) ! W
      Eth = 220.5d0
    case(6)  ! C
      Eth = 27.64d0
    case(42) ! Mo
      Eth = 180.0d0
    case(26) ! Fe
      Eth = 120.0d0
    end select
  elseif (Z_ion == 2) then ! He
    select case(Z_imp)
    case(74) ! W
      Eth = 62.06d0
    case(6)  ! C
      Eth = 52.98d0
    end select
  elseif (Z_ion == Z_imp) then ! Self-sputtering
    Eth = 4.0d0 * get_Es(Z_imp)  ! Rough estimate
  else
    ! Use semi-empirical formula
    Eth = 4.0d0 * get_Es(Z_imp) * &
          (atomic_weights(Z_ion) + atomic_weights(Z_imp))**2 / &
          (4.0d0 * atomic_weights(Z_ion) * atomic_weights(Z_imp))
  end if
end function get_Eth

!> Get surface coordination number
pure function get_ns(Z_imp) result(ns)
  integer, intent(in) :: Z_imp
  real*8 :: ns   ! surface coordination number
  
  select case(Z_imp)
  case(74) ! W
    ns = 0.06325d0
  case(6)  ! C
    ns = 0.11286d0
  case(42) ! Mo
    ns = 0.075d0
  case(26) ! Fe
    ns = 0.085d0
  case default
    ns = 0.1d0  ! Default value
  end select
end function get_ns

!> Thompson energy distribution for sputtered particles
subroutine ThompsonEne(E_sput, U0, rng, E)
  use mod_rng, only: type_rng
  implicit none
  
  real*8, intent(in)                    :: E_sput  !< Sputtering energy [eV]
  real*8, intent(in)                    :: U0      !< Surface binding energy [eV]
  class(type_rng), intent(inout)        :: rng     !< Random number generator
  real*8, intent(out)                   :: E       !< Energy of sputtered particle [eV]
  
  real*8 :: rng_samples(2)  !< Two random numbers for rejection sampling
  real*8 :: E_max, P_max, P, E_test
  
  ! No sputtering if incident energy is below threshold or near zero
  if (E_sput .le. U0 .or. E_sput .le. 0.d0) then
    E = 0.d0
    return
  end if

  ! Maximum energy transfer (maximum energy of sputtered particle)
  E_max = E_sput
  
  ! Generate random energy from Thompson distribution
  ! Probability density function: f(E) = E / (E + U0)^3
  
  ! Maximum of f(E) occurs at E = U0/2
  ! f_max = f(U0/2) = (U0/2) / (U0/2 + U0)^3 = (U0/2) / (3U0/2)^3 = 1/(27U0^2)
  if (U0 > 0.0d0) then
    P_max = 1.0d0 / (27.0d0 * U0**2)
  else
    P_max = 1.0d0  ! Avoid division by zero if U0=0
  end if
  
  ! Use rejection method to sample from Thompson distribution
  rejection_loop: do
    ! Get two random numbers from the RNG
    call rng%next(rng_samples)
    
    ! First random number determines trial energy (uniform in [0, E_max])
    E_test = rng_samples(1) * E_max
    
    ! Calculate probability density at trial energy
    if (E_test > 0.0d0 .and. U0 > 0.0d0) then
      P = E_test / (E_test + U0)**3
    else
      P = 0.0d0
    end if
    
    ! Second random number determines acceptance/rejection
    if (rng_samples(2) * P_max < P) then
      E = E_test
      exit rejection_loop
    end if
  end do rejection_loop
end subroutine ThompsonEne

end module mod_eckstein_thompson