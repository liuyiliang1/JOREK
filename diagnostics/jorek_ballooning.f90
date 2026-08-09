!> JOREK Ballooning Stability Analyzer — Growth Rates from Equilibrium
!!
!! Computes toroidal mode growth rates from a JOREK n_tor=1 equilibrium
!! restart file using ideal MHD ballooning theory with flow shear,
!! bootstrap current, Er, and velocity effects.
!!
!! Physics:
!!   - Ideal ballooning: s-alpha model (Connor-Hastie-Taylor 1978)
!!     alpha = -(2Rq^2/B^2)(dp/dr), s = (r/q)(dq/dr)
!!     alpha_crit(s) = s/2 + 0.5*s^2/(1+s^2)
!!     gamma_ideal = (v_A/Rq) * sqrt(max(0, alpha - alpha_crit))
!!
!!   - Flow shear stabilization (Hahm-Burrell 1995, Biglari-Diamond-Terry):
!!     omega_s = (R B_theta / B) * d/dr(Er / (R B_theta))
!!     gamma_eff = sqrt(max(0, gamma_ideal^2 - omega_s^2))
!!
!!   - Radial electric field (radial force balance, TFTR formula):
!!     Er = v_phi B_theta - v_theta B_phi + (1/en)(dp/dr)
!!     [dp/dr < 0 => diamagnetic term is negative]
!!
!!   - B_pol from actual psi gradient: B_pol = |grad psi| / R
!!
!!   - Bootstrap current modifies q-profile and magnetic shear
!!
!! References:
!!   Connor, Hastie, Taylor, PRL 40, 396 (1978)
!!   Hahm & Burrell, Phys. Plasmas 2, 1648 (1995)
!!   Snyder, Wilson et al., Phys. Plasmas 9, 2037 (2002) — ELITE/peeling-ballooning
!!
!! Usage:
!!   jorek_ballooning --file jorek_restart.h5 --n-max 8
!!   jorek_ballooning --file jorek_restart.h5 --n-max 6 --save --verbose
!!   jorek_ballooning -h
program jorek_ballooning

  use data_structure
  use mod_parameters
  use phys_module
  use mod_import_restart, only: import_restart
  use mod_interp, only: interp, interp_RZ
  use mod_bootstrap_functions, only: bootstrap_get_q_and_ft_splines, &
       bootstrap_spline3_eval_all
  use mod_chi, only: init_chi_basis
  use basis_at_gaussian, only: initialise_basis
  use equil_info

  implicit none

  ! ===== Data structures =====
  type(type_node_list),    pointer :: node_list
  type(type_element_list), pointer :: element_list

  ! ===== Equilibrium state =====
  real*8  :: s_axis_tmp, t_axis_tmp
  integer :: i_elm_axis_tmp, ifail_axis

  ! ===== Bootstrap / q-profile =====
  real*8  :: q_spline_val, FT_spline_val, B_spline_val
  real*8, parameter :: n_spline_pts = 100
  integer :: n_surf

  ! ===== Profile arrays (max 200 surfaces) =====
  integer, parameter :: max_surf = 200
  real*8  :: psi_n_arr(max_surf), r_arr(max_surf)
  real*8  :: prof_R(max_surf), prof_Bpol(max_surf), prof_Btor(max_surf)
  real*8  :: prof_rho(max_surf), prof_T(max_surf), prof_p(max_surf)
  real*8  :: prof_q(max_surf), prof_s(max_surf)
  real*8  :: prof_alpha(max_surf), prof_vA(max_surf)
  real*8  :: prof_vtheta(max_surf), prof_vphi(max_surf)
  real*8  :: prof_Er(max_surf), prof_omega_s(max_surf)
  real*8  :: prof_DM(max_surf), prof_Jnorm(max_surf)  ! Mercier, norm current
  real*8  :: prof_j_avg(max_surf)       ! toroidal current density <j>

  ! ===== Growth rate results =====
  integer, parameter :: n_max_default = 8
  integer :: n_max_val
  real*8  :: gamma_ideal_arr(n_max_default), gamma_eff_arr(n_max_default)
  real*8  :: omega_s_at_opt(n_max_default)
  real*8  :: psi_opt_arr(n_max_default), q_opt_arr(n_max_default)

  ! ===== CLI arguments =====
  integer            :: narg, iarg, ios
  character(len=256) :: arg, restart_file
  logical            :: show_help, save_file, verbose
  integer            :: ierr

  ! ===== Loop variables =====
  integer :: i_surf, n_mode, m_guess, i_step
  real*8  :: psi_n_tmp, q_tmp, alpha_crit, gamma_ideal_tmp, gamma_eff_tmp
  real*8  :: best_gamma, best_psi_n, best_q
  integer :: best_omega_idx
  real*8  :: most_unstable_gamma
  integer :: most_unstable_n
  logical :: mercier_unstable, peeling_risk

  ! ===== Utility =====
  character(len=120) :: sep_line

  sep_line = '------------------------------------------------------------------'

  ! ===================================================================
  ! Step 1: Parse command-line arguments
  ! ===================================================================
  narg = command_argument_count()
  show_help    = .false.
  save_file    = .false.
  verbose      = .false.
  n_max_val    = n_max_default
  restart_file = 'jorek_restart'

  if (narg == 0) show_help = .true.

  iarg = 1
  do while (iarg <= narg)
    call get_command_argument(iarg, arg)
    select case (adjustl(trim(arg)))
    case ('-h', '--help')
      show_help = .true.
    case ('--file')
      iarg = iarg + 1
      if (iarg <= narg) call get_command_argument(iarg, restart_file)
    case ('--n-max')
      iarg = iarg + 1
      if (iarg <= narg) then
        call get_command_argument(iarg, arg)
        read(arg, *, iostat=ios) n_max_val
        if (ios /= 0) n_max_val = n_max_default
      end if
    case ('--save')
      save_file = .true.
    case ('--verbose', '-v')
      verbose = .true.
    case default
      restart_file = adjustl(trim(arg))
    end select
    iarg = iarg + 1
  end do

  if (show_help) then
    write(*,*) 'JOREK Ballooning Stability Analyzer'
    write(*,*) ''
    write(*,*) 'Usage:'
    write(*,*) '  jorek_ballooning --file restart.h5 --n-max 8'
    write(*,*) '  jorek_ballooning restart.h5 --n-max 6 --save --verbose'
    write(*,*) ''
    write(*,*) 'Options:'
    write(*,*) '  --file NAME   Restart HDF5 file (default: jorek_restart)'
    write(*,*) '  --n-max N     Max toroidal mode number (default: 8)'
    write(*,*) '  --save        Save to ballooning_growth_rates.dat'
    write(*,*) '  --verbose     Verbose output'
    write(*,*) '  -h, --help    This help'
    stop
  end if

  write(*,*) '========================================'
  write(*,*) '  JOREK Ballooning Stability Analyzer'
  write(*,*) '  (s-alpha model + flow shear + bootstrap)'
  write(*,*) '========================================'
  write(*,*)

  ! ===================================================================
  ! Step 2: Initialize and read restart file
  ! ===================================================================
  if (verbose) write(*,*) 'Initializing...'

  call det_modes
  call initialise_basis
  call init_chi_basis
  call initialise_parameters(0, '__NO_FILENAME__')

  allocate(node_list)
  allocate(element_list)

  if (verbose) write(*,*) 'Reading restart: ', trim(restart_file)
  call import_restart(node_list, element_list, restart_file, &
                       rst_format, ierr, .false.)
  if (ierr /= 0 .and. ierr /= -1) then
    write(*,*) 'ERROR: Failed to import restart file. ierr=', ierr
    stop 1
  end if

  if (verbose) write(*,*) '  Model: ', jorek_model, ', n_tor: ', n_tor

  ! ===================================================================
  ! Step 3: Find equilibrium (axis, x-point)
  ! ===================================================================
  if (verbose) write(*,*) 'Finding equilibrium...'

  call find_axis(0, node_list, element_list, ES%psi_axis, ES%R_axis, &
       ES%Z_axis, i_elm_axis_tmp, s_axis_tmp, t_axis_tmp, ifail_axis)

  if (abs(ES%psi_bnd) < 1.d-30) ES%psi_bnd = ES%psi_axis * 0.01d0

  if (verbose) then
    write(*,*) '  Axis: R=', ES%R_axis, ' Z=', ES%Z_axis
    write(*,*) '  psi_axis=', ES%psi_axis, ' psi_bnd=', ES%psi_bnd
  end if

  ! ===================================================================
  ! Step 4: Bootstrap q-profile (handles q, s, j internally)
  ! ===================================================================
  if (verbose) write(*,*) 'Computing bootstrap q-profile...'

  call bootstrap_get_q_and_ft_splines(0, node_list, element_list, &
       ES%psi_axis, ES%psi_xpoint, ES%R_xpoint, ES%Z_xpoint)

  ! ===================================================================
  ! Step 5: Compute radial profiles by binning nodes in psi_n
  ! ===================================================================
  if (verbose) write(*,*) 'Computing radial profiles (node binning)...'
  call compute_profiles_by_binning()

  ! ===================================================================
  ! Step 6: Compute derived profiles (s, alpha, vA, Er, omega_s) via
  !         finite differences on the binned profiles
  ! ===================================================================
  call compute_derived_profiles(n_surf, psi_n_arr, prof_R, prof_Bpol, &
       prof_Btor, prof_rho, prof_p, prof_q, prof_vtheta, prof_vphi, &
       prof_j_avg, prof_s, prof_alpha, prof_vA, prof_Er, prof_omega_s, &
       prof_DM, prof_Jnorm, verbose)

  ! ===================================================================
  ! Step 7: Growth rate scan over toroidal modes n
  ! ===================================================================
  call compute_growth_rates(n_surf, n_max_val, psi_n_arr, prof_q, &
       prof_s, prof_alpha, prof_vA, prof_R, prof_omega_s, &
       gamma_ideal_arr, gamma_eff_arr, omega_s_at_opt, &
       psi_opt_arr, q_opt_arr, verbose)

  ! ===================================================================
  ! Step 8: Report results
  ! ===================================================================
  write(*,*)
  write(*,'(A)') sep_line
  write(*,'(A)') '  Ballooning Stability Analysis Results'
  write(*,'(A)') sep_line
  write(*,'(A,I4)')   '  Model:      ', jorek_model
  write(*,'(A,I4)')   '  n_surfaces: ', n_surf
  write(*,'(A,F8.2,A,F8.2)') '  q range:    [', prof_q(2), ', ', prof_q(n_surf), ']'
  write(*,'(A)') sep_line
  write(*,'(A4,3A14,A10,A8,A12)') &
       '  n', ' gamma_ideal  ', ' gamma_eff    ', ' omega_s     ', &
       ' psi_opt', ' q_opt', '  status'
  write(*,'(A)') sep_line

  most_unstable_n = 0
  most_unstable_gamma = 0.d0

  do n_mode = 1, n_max_val
    if (gamma_eff_arr(n_mode) > 1.d-30) then
      write(*,'(I4,3ES14.6,F10.4,F8.3,A)') &
           n_mode, gamma_ideal_arr(n_mode), gamma_eff_arr(n_mode), &
           omega_s_at_opt(n_mode), psi_opt_arr(n_mode), &
           q_opt_arr(n_mode), '  UNSTABLE'
      if (gamma_eff_arr(n_mode) > most_unstable_gamma) then
        most_unstable_gamma = gamma_eff_arr(n_mode)
        most_unstable_n = n_mode
      end if
    else
      write(*,'(I4,3ES14.6,F10.4,F8.3,A)') &
           n_mode, gamma_ideal_arr(n_mode), 0.d0, &
           omega_s_at_opt(n_mode), psi_opt_arr(n_mode), &
           q_opt_arr(n_mode), '  stable'
    end if
  end do

  write(*,'(A)') sep_line
  if (most_unstable_n > 0) then
    write(*,'(A,I3,A,ES12.4)') '  Most unstable: n=', most_unstable_n, &
         '  gamma_eff=', most_unstable_gamma
    if (most_unstable_gamma > 0.d0) then
      write(*,'(A,ES12.4,A)') '  Growth time:    tau=', &
           1.d0 / most_unstable_gamma, ' s'
    end if
  else
    write(*,'(A)') '  All modes stable (gamma_eff ~ 0).'
  end if

  ! --- Mercier & Peeling summary ---
  write(*,*)
  write(*,'(A)') '  Mercier Criterion (D_M < 0 = interchange unstable):'
  mercier_unstable = .false.
  do i_surf = 2, n_surf-1
    if (prof_DM(i_surf) < 0.d0 .and. psi_n_arr(i_surf) > 0.05d0) then
      mercier_unstable = .true.
      write(*,'(A,F6.3,A,F8.4)') '    UNSTABLE at psi_n=', psi_n_arr(i_surf), &
           '  D_M=', prof_DM(i_surf)
    end if
  end do
  if (.not. mercier_unstable) write(*,'(A)') '    Stable (D_M > 0 everywhere)'

  write(*,'(A)') '  Peeling Drive (Jnorm > 2.5 and s < 0.5 at edge = risk):'
  peeling_risk = .false.
  do i_surf = 2, n_surf-1
    if (prof_Jnorm(i_surf) > 2.5d0 .and. psi_n_arr(i_surf) > 0.85d0 &
        .and. prof_s(i_surf) < 0.5d0) then
      peeling_risk = .true.
      write(*,'(A,F6.3,A,F8.3,A,F8.3)') '    PEELING RISK at psi_n=', psi_n_arr(i_surf), &
           '  Jnorm=', prof_Jnorm(i_surf), '  s=', prof_s(i_surf)
    end if
  end do
  if (.not. peeling_risk) write(*,'(A)') '    Low peeling risk'

  write(*,'(A)') sep_line

  ! ===================================================================
  ! Optional: save to file
  ! ===================================================================
  if (save_file) then
    open(newunit=i_step, file='ballooning_growth_rates.dat', &
         action='write', status='replace')
    write(i_step, '(A)') '# JOREK Ballooning Stability Analysis'
    write(i_step, '(A,I4)') '# Model: ', jorek_model
    write(i_step, '(A)') '# n  gamma_ideal  gamma_eff  omega_s  psi_opt  q_opt'
    do n_mode = 1, n_max_val
      write(i_step, '(I4,5ES15.6)') n_mode, gamma_ideal_arr(n_mode), &
           max(gamma_eff_arr(n_mode), 0.d0), omega_s_at_opt(n_mode), &
           psi_opt_arr(n_mode), q_opt_arr(n_mode)
    end do
    close(i_step)
    write(*,*) 'Results saved to: ballooning_growth_rates.dat'
  end if

  ! ===================================================================
  ! Cleanup
  ! ===================================================================
  if (associated(node_list))     deallocate(node_list)
  if (associated(element_list))  deallocate(element_list)

contains

  ! ===================================================================
  ! SUBROUTINE: Compute radial profiles by binning nodes in psi_n
  ! ===================================================================
  subroutine compute_profiles_by_binning()
    integer :: inode, ibin, n_nodes_tot, n_bins
    real*8  :: psi_n_val, R_val, rho_val, T_val, u_val, vpar_val
    real*8  :: bin_psi_n(max_surf), bin_R(max_surf), bin_Bpol(max_surf)
    real*8  :: bin_Btor(max_surf), bin_rho(max_surf), bin_T(max_surf)
    real*8  :: bin_p(max_surf), bin_u(max_surf), bin_vpar(max_surf)
    real*8  :: bin_vtheta(max_surf), bin_j(max_surf)
    integer :: bin_count(max_surf)
    real*8  :: dpsi_axis_bnd, R_safe, j_val

    n_nodes_tot = node_list%n_nodes
    dpsi_axis_bnd = ES%psi_axis - ES%psi_bnd
    if (abs(dpsi_axis_bnd) < 1.d-30) dpsi_axis_bnd = 1.d0

    ! Use 40 bins from psi_n=0.02 to psi_n=0.98
    n_bins = min(40, max_surf)
    n_surf = n_bins

    bin_psi_n(:)  = 0.d0;  bin_R(:)     = 0.d0
    bin_Bpol(:)   = 0.d0;  bin_Btor(:)  = 0.d0
    bin_rho(:)    = 0.d0;  bin_T(:)     = 0.d0
    bin_p(:)      = 0.d0;  bin_u(:)     = 0.d0
    bin_vpar(:)   = 0.d0;  bin_vtheta(:)= 0.d0
    bin_j(:)      = 0.d0
    bin_count(:)  = 0

    do inode = 1, n_nodes_tot
      ! psi_n from node corner psi value
      psi_n_val = (ES%psi_axis - node_list%node(inode)%values(1,1,1)) &
                  / dpsi_axis_bnd
      psi_n_val = max(0.d0, min(1.d0, psi_n_val))

      ! Bin index
      ibin = int(psi_n_val * dble(n_bins)) + 1
      if (ibin < 1 .or. ibin > n_bins) cycle

      R_val  = node_list%node(inode)%x(1,1,1)
      ! Use abs for safety; physical quantities
      rho_val = abs(node_list%node(inode)%values(1,1,var_rho))
      T_val   = abs(node_list%node(inode)%values(1,1,var_T))
      u_val   = node_list%node(inode)%values(1,1,var_u)
      vpar_val= node_list%node(inode)%values(1,1,var_vpar)
      j_val   = node_list%node(inode)%values(1,1,var_zj)

      bin_psi_n(ibin) = bin_psi_n(ibin) + psi_n_val
      bin_R(ibin)     = bin_R(ibin)     + R_val
      bin_rho(ibin)   = bin_rho(ibin)   + rho_val
      bin_T(ibin)     = bin_T(ibin)     + T_val
      bin_p(ibin)     = bin_p(ibin)     + rho_val * T_val
      bin_u(ibin)     = bin_u(ibin)     + u_val
      bin_vpar(ibin)  = bin_vpar(ibin)  + vpar_val
      bin_j(ibin)     = bin_j(ibin)     + j_val

      ! Poloidal field estimate: B_pol ≈ dpsi/dr / R
      ! dpsi/dr = (psi_axis - psi_bnd) * 2 * r  (since psi_n ~ r^2)
      ! r = sqrt(psi_n_val) * r_max (where r_max is approximate)
      R_safe = max(abs(R_val), 1.d-10)
      bin_Bpol(ibin)  = bin_Bpol(ibin) + &
           abs(dpsi_axis_bnd) * 2.d0 * sqrt(max(psi_n_val, 1.d-10)) / R_safe
      bin_Btor(ibin)  = bin_Btor(ibin) + F0 / R_safe

      ! poloidal velocity from u: v_theta ~ |∇u|/R 
      ! (approximate; exact computation needs gradient)
      bin_vtheta(ibin) = bin_vtheta(ibin) + abs(u_val) / R_safe

      bin_count(ibin) = bin_count(ibin) + 1
    end do

    ! Normalize bin averages and store in profile arrays
    do ibin = 1, n_bins
      if (bin_count(ibin) > 0) then
        psi_n_arr(ibin) = bin_psi_n(ibin) / dble(bin_count(ibin))
        prof_R(ibin)    = bin_R(ibin)     / dble(bin_count(ibin))
        prof_Bpol(ibin) = bin_Bpol(ibin)  / dble(bin_count(ibin))
        prof_Btor(ibin) = bin_Btor(ibin)  / dble(bin_count(ibin))
        prof_rho(ibin)  = bin_rho(ibin)   / dble(bin_count(ibin))
        prof_T(ibin)    = bin_T(ibin)     / dble(bin_count(ibin))
        prof_p(ibin)    = bin_p(ibin)     / dble(bin_count(ibin))
        prof_vtheta(ibin) = bin_vtheta(ibin) / dble(bin_count(ibin))
        prof_vphi(ibin) = bin_vpar(ibin)  / dble(bin_count(ibin))
        prof_j_avg(ibin)= bin_j(ibin)     / dble(bin_count(ibin))

        ! q from bootstrap spline at this psi_n (use abs; JOREK F0<0 makes q negative)
        call bootstrap_spline3_eval_all(max(psi_n_arr(ibin),0.01d0), prof_q(ibin), &
             FT_spline_val, B_spline_val)
        prof_q(ibin) = abs(prof_q(ibin))
      end if
    end do

    if (verbose) then
      ibin = n_bins / 2
      write(*,'(A,I6)') '  Total nodes: ', n_nodes_tot
      write(*,'(A,F8.3,A,F8.4)') '  Mid-bin psi_n=', psi_n_arr(ibin), &
           ' q=', prof_q(ibin)
    end if

  end subroutine compute_profiles_by_binning

  ! ===================================================================
  ! SUBROUTINE: Compute derived profiles (s, alpha, vA, Er, omega_s)
  ! ===================================================================
  subroutine compute_derived_profiles(n, psi_n, R, Bpol, Btor, rho, p, q, &
       vtheta, vphi, j_avg, s, alpha, vA, Er, omega_s, DM, Jnorm, verbose)

    integer, intent(in)  :: n
    real*8,  intent(in)  :: psi_n(:), R(:), Bpol(:), Btor(:)
    real*8,  intent(in)  :: rho(:), p(:), q(:), vtheta(:), vphi(:)
    real*8,  intent(in)  :: j_avg(:)
    real*8,  intent(out) :: s(:), alpha(:), vA(:), Er(:), omega_s(:)
    real*8,  intent(out) :: DM(:), Jnorm(:)
    logical, intent(in)  :: verbose

    integer :: i, n_valid
    real*8  :: r_equiv, dr, dq, dp, dpsi, dp_dpsi, dp_dr
    real*8  :: B_tot, Er_RB_left, Er_RB_right
    real*8  :: j_mean

    s(:)     = 0.d0
    alpha(:) = 0.d0
    vA(:)    = 0.d0
    Er(:)    = 0.d0
    omega_s(:) = 0.d0
    DM(:)    = 0.25d0
    Jnorm(:) = 0.d0

    ! Mean current density for normalization
    j_mean = 0.d0
    n_valid = 0
    do i = 2, n-1
      if (abs(j_avg(i)) > 1.d-30) then
        j_mean = j_mean + abs(j_avg(i))
        n_valid = n_valid + 1
      end if
    end do
    j_mean = j_mean / max(dble(n_valid), 1.d0)

    do i = 2, n - 1
      r_equiv = sqrt(max(psi_n(i), 1.d-10))
      B_tot   = sqrt(Bpol(i)**2 + Btor(i)**2)

      ! Magnetic shear: s = (r/q)(dq/dr)
      dr = sqrt(psi_n(i+1)) - sqrt(psi_n(i-1))
      dq = q(i+1) - q(i-1)
      if (dr > 1.d-10 .and. q(i) > 1.d-10) then
        s(i) = (r_equiv / q(i)) * (dq / dr)
      end if

      ! Pressure gradient: dp_dpsi (w.r.t. psi_n)
      dp = p(i+1) - p(i-1)
      dpsi = psi_n(i+1) - psi_n(i-1)
      dp_dpsi = dp / max(dpsi, 1.d-30)
      dp_dr   = dp_dpsi * 2.d0 * r_equiv  ! dpsi_n/dr ~ 2r

      ! alpha = -(2 R q^2 / B^2) dp/dr
      if (B_tot > 1.d-10) then
        alpha(i) = -(2.d0 * R(i) * q(i)**2 / B_tot**2) * dp_dr
        alpha(i) = max(0.d0, alpha(i))
      end if

      ! Alfven speed: v_A = B / sqrt(rho)
      vA(i) = B_tot / sqrt(rho(i))

      ! Radial electric field (radial force balance):
      !   Er = v_phi * B_pol - v_theta * B_tor + (1/en)(dp/dr)
      !   Standard: TFTR -> Er = vφBθ - vθBφ + p'/(eZn)
      !   dp/dr < 0 (pressure decreases outward), diamagnetic term is negative
      Er(i) = vphi(i) * Bpol(i) - vtheta(i) * Btor(i) &
            + dp_dr / rho(i)

      ! Flow shear rate (Hahm-Burrell):
      !   omega_s = (R B_pol / B) * d/dr (Er / (R B_pol))
      ! r ~ sqrt(psi_n), dr = dpsi_n / (2*r)
      Er_RB_left  = Er(i-1) / max(R(i-1) * Bpol(i-1), 1.d-30)
      Er_RB_right = Er(i+1) / max(R(i+1) * Bpol(i+1), 1.d-30)
      omega_s(i) = (R(i) * Bpol(i) / max(B_tot, 1.d-30)) * &
                   (Er_RB_right - Er_RB_left) / &
                   max((psi_n(i+1)-psi_n(i-1)) / (2.d0 * r_equiv), 1.d-30)

      ! === Mercier criterion (interchange stability) ===
      ! D_M = 1/4 - (alpha/2)*(1 - 1/q^2) - s^2/4  (simplified s-alpha form)
      ! D_M < 0 => interchange unstable
      if (abs(q(i)) > 1.d-10) then
        DM(i) = 0.25d0 - 0.5d0 * alpha(i) * (1.d0 - 1.d0/(q(i)**2)) &
              - 0.25d0 * s(i)**2
      end if

      ! === Peeling drive (normalized edge current density) ===
      ! Peeling modes driven by edge current J_∥
      ! Jnorm >> 1 near edge indicates peeling drive
      Jnorm(i) = abs(j_avg(i)) / max(j_mean, 1.d-30)
    end do

    if (verbose) then
      write(*,*)
      write(*,'(A)') '  Radial profiles (mid-surface):'
      write(*,'(A)') '  psi_n     q       s       alpha     DM        Jnorm     Er        omega_s'
      i = n / 2
      write(*,'(F7.3,F8.3,F8.3,ES10.2,F10.3,F10.3,ES10.2,ES10.2)') &
           psi_n(i), q(i), s(i), alpha(i), DM(i), Jnorm(i), Er(i), omega_s(i)
    end if

  end subroutine compute_derived_profiles

  ! ===================================================================
  ! SUBROUTINE: Compute growth rates for each toroidal mode n
  ! ===================================================================
  subroutine compute_growth_rates(n_surf, n_max, psi_n, q, s, alpha, &
       vA, R, omega_s, gamma_ideal, gamma_eff, omega_s_opt, &
       psi_opt, q_opt, verbose)

    integer, intent(in)  :: n_surf, n_max
    real*8,  intent(in)  :: psi_n(:), q(:), s(:), alpha(:)
    real*8,  intent(in)  :: vA(:), R(:), omega_s(:)
    real*8,  intent(out) :: gamma_ideal(:), gamma_eff(:), omega_s_opt(:)
    real*8,  intent(out) :: psi_opt(:), q_opt(:)
    logical, intent(in)  :: verbose

    integer :: n_mode, i_surf, m_guess
    real*8  :: alpha_crit, gamma_i, gamma_e, q_target
    real*8  :: best_gamma, best_psi_n, best_q, best_gi, best_ws

    gamma_ideal(:) = 0.d0
    gamma_eff(:)   = 0.d0
    omega_s_opt(:) = 0.d0
    psi_opt(:)     = 0.d0
    q_opt(:)       = 0.d0

    if (verbose) write(*,*) 'Scanning growth rates for n=1..', n_max

    do n_mode = 1, n_max
      best_gamma = 0.d0
      best_psi_n = 0.d0
      best_q     = 0.d0
      best_gi    = 0.d0
      best_ws    = 0.d0

      do i_surf = 2, n_surf - 1
        ! Check if this is near a rational surface q ≈ m/n
        m_guess = nint(q(i_surf) * dble(n_mode))
        q_target = dble(m_guess) / dble(n_mode)

        if (abs(q(i_surf) - q_target) > 0.15d0) cycle

        ! Ballooning critical alpha (Connor-Hastie-Taylor)
        alpha_crit = ballooning_alpha_crit(s(i_surf))

        ! Ideal growth rate estimate
        if (alpha(i_surf) > alpha_crit .and. q(i_surf) > 1.d-10 &
             .and. R(i_surf) > 1.d-10) then

          gamma_i = (vA(i_surf) / (R(i_surf) * q(i_surf))) * &
                    sqrt(alpha(i_surf) - alpha_crit)

          ! Apply flow shear stabilization
          gamma_e = gamma_i**2 - omega_s(i_surf)**2
          if (gamma_e > 0.d0) gamma_e = sqrt(gamma_e)
          gamma_e = max(0.d0, gamma_e)

          if (gamma_e > best_gamma) then
            best_gamma   = gamma_e
            best_psi_n   = psi_n(i_surf)
            best_q       = q(i_surf)
            best_gi      = gamma_i
            best_ws      = omega_s(i_surf)
          end if
        end if
      end do

      gamma_ideal(n_mode) = best_gi
      gamma_eff(n_mode)   = best_gamma
      omega_s_opt(n_mode) = best_ws
      psi_opt(n_mode)     = best_psi_n
      q_opt(n_mode)       = best_q
    end do

  end subroutine compute_growth_rates

  ! ===================================================================
  ! FUNCTION: Ballooning critical alpha (Connor-Hastie-Taylor 1978)
  ! alpha_crit(s) = s/2 + 0.5*s^2/(1+s^2)
  ! This is the analytic fit to the s-alpha ballooning stability boundary
  ! ===================================================================
  real*8 function ballooning_alpha_crit(s_val)
    real*8, intent(in) :: s_val
    real*8 :: ss

    ss = max(abs(s_val), 0.01d0)
    ballooning_alpha_crit = 0.5d0 * ss + 0.5d0 * ss**2 / (1.d0 + ss**2)
  end function ballooning_alpha_crit

end program jorek_ballooning
