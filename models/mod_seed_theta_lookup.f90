!> Fast lookup of straight-field-line angle theta* for seed island injection.
!!
!! Pre-computes the mapping theta_geom -> theta_star by field line tracing
!! in the n=0 equilibrium. Results stored on an equidistant (psi_n, theta_geom)
!! grid for O(1) bilinear interpolation during element matrix assembly.
!!
!! Thread-safe: the lookup table is read-only after initialization.
module mod_seed_theta_lookup

  use constants, only: PI
  implicit none

  private
  public :: seed_theta_lookup_initialized
  public :: seed_init_theta_lookup, seed_free_theta_lookup
  public :: seed_lookup_theta_star

  logical, save :: seed_theta_lookup_initialized = .false.

  ! Grid dimensions
  integer, parameter :: N_PSI   = 101   !< radial grid points
  integer, parameter :: N_THETA = 128   !< poloidal grid points

  ! Lookup table: theta_star(j_theta, k_psi)
  real*8, save, allocatable :: theta_star_table(:,:)

  ! --- Explicit interface for external find_RZ subroutine ---
  interface
    subroutine find_RZ(node_list, element_list, R_find, Z_find, R_out, Z_out, &
                       ielm_out, s_out, t_out, ifail)
      use data_structure, only: type_node_list, type_element_list
      type(type_node_list),    intent(in)  :: node_list
      type(type_element_list), intent(in)  :: element_list
      real*8,                  intent(in)  :: R_find, Z_find
      real*8,                  intent(out) :: R_out, Z_out
      integer,                 intent(out) :: ielm_out
      real*8,                  intent(out) :: s_out, t_out
      integer,                 intent(out) :: ifail
    end subroutine find_RZ
  end interface

contains

  !=============================================================================
  !> Initialize the theta* lookup table by tracing field lines.
  !! Called once during initial_conditions, before the time loop.
  !=============================================================================
  subroutine seed_init_theta_lookup(node_list, element_list, equil_state)
    use data_structure, only: type_node_list, type_element_list
    use equil_info, only: t_equil_state
    implicit none
    type(type_node_list),    intent(in) :: node_list
    type(type_element_list), intent(in) :: element_list
    type(t_equil_state),     intent(in) :: equil_state

    real*8  :: PsiNRange(2)
    integer :: ierr

    if (seed_theta_lookup_initialized) return

    PsiNRange(1) = 0.02d0
    PsiNRange(2) = 0.98d0

    if (allocated(theta_star_table)) deallocate(theta_star_table)
    allocate(theta_star_table(N_THETA, N_PSI))

    call seed_trace_fieldlines(node_list, element_list, equil_state, &
         PsiNRange, N_PSI, N_THETA, ierr)

    if (ierr /= 0) then
      write(*,*) 'WARNING: seed_init_theta_lookup failed (ierr=', ierr, ')'
      write(*,*) '  Falling back to geometric poloidal angle.'
      deallocate(theta_star_table)
      seed_theta_lookup_initialized = .false.
      return
    end if

    seed_theta_lookup_initialized = .true.
    write(*,'(A)') 'seed theta* lookup table initialized.'
    write(*,'(A,I4,A,I4)') '  Grid: ', N_PSI, ' x ', N_THETA

  end subroutine seed_init_theta_lookup


  !=============================================================================
  !> Free the lookup table memory.
  !=============================================================================
  subroutine seed_free_theta_lookup()
    implicit none
    if (allocated(theta_star_table)) deallocate(theta_star_table)
    seed_theta_lookup_initialized = .false.
  end subroutine seed_free_theta_lookup


  !=============================================================================
  !> Fast bilinear interpolation lookup:  theta_star = f(psi_n, theta_geom)
  !!
  !! Called at every Gauss point during element matrix assembly.
  !! Thread-safe: read-only access to module-level theta_star_table.
  !=============================================================================
  function seed_lookup_theta_star(psi_n, theta_geom) result(theta_star)
    implicit none
    real*8, intent(in) :: psi_n, theta_geom
    real*8 :: theta_star

    real*8  :: psi_c, th_c, psi_frac, th_frac
    integer :: ip, ith

    if (.not. seed_theta_lookup_initialized) then
      theta_star = theta_geom
      return
    end if
    if (.not. allocated(theta_star_table)) then
      theta_star = theta_geom
      return
    end if

    ! Clamp psi_n to [0, 1]
    psi_c = min(max(psi_n, 0.d0), 1.d0)

    ! Map theta_geom to [0, 2*pi)
    th_c = mod(theta_geom, 2.d0 * PI)
    if (th_c < 0.d0) th_c = th_c + 2.d0 * PI

    ! Compute fractional indices in [1, N_PSI] and [1, N_THETA]
    psi_frac = psi_c * dble(N_PSI - 1)
    ip       = max(1, min(N_PSI - 1, int(psi_frac) + 1))
    psi_frac = psi_frac - dble(ip - 1)

    th_frac  = th_c / (2.d0 * PI) * dble(N_THETA - 1)
    ith      = max(1, min(N_THETA - 1, int(th_frac) + 1))
    th_frac  = th_frac - dble(ith - 1)

    ! Bilinear interpolation
    theta_star = (1.d0-psi_frac)*(1.d0-th_frac) * theta_star_table(ith,   ip)   &
               + (1.d0-psi_frac)*      th_frac  * theta_star_table(ith+1, ip)   &
               +       psi_frac *(1.d0-th_frac) * theta_star_table(ith,   ip+1) &
               +       psi_frac *      th_frac  * theta_star_table(ith+1, ip+1)

  end function seed_lookup_theta_star


  !=============================================================================
  ! PRIVATE routines below
  !=============================================================================


  !=============================================================================
  !> Trace field lines in the n=0 equilibrium to build the theta_geom -> theta_star
  !! mapping on each flux surface. Uses predictor-corrector stepping in phi.
  !!
  !! Algorithm adapted from mod_straight_field_line::trace_fieldlines.
  !=============================================================================
  subroutine seed_trace_fieldlines(node_list, element_list, equil_state, &
       PsiNRange, nPsiN, nTht, ierr)
    use data_structure, only: type_node_list, type_element_list
    use equil_info, only: t_equil_state, get_psi_n
    use phys_module, only: F0
    use mod_interp, only: interp_RZ
    use mod_basisfunctions, only: basisfunctions
    use mod_parameters, only: n_vertex_max, n_degrees

    implicit none

    type(type_node_list),    intent(in)    :: node_list
    type(type_element_list), intent(in)    :: element_list
    type(t_equil_state),     intent(in)    :: equil_state
    real*8,                  intent(in)    :: PsiNRange(2)
    integer,                 intent(in)    :: nPsiN, nTht
    integer,                 intent(inout) :: ierr

    ! --- Local constants ---
    real*8,  parameter :: DELTAPHI  = 0.01d0   ! toroidal step per large step [rad]
    integer, parameter :: NMAXSTEPS = 10000     ! max steps per surface

    ! --- Local variables ---
    integer :: k, j, step_count, npts
    real*8  :: psi_n_k
    real*8  :: Rleft, Rright, Rmid, Rout, Zout
    integer :: ielm, ifail
    real*8  :: s, t, P, P_s, P_t
    real*8  :: rn, zn, rh, zh, rp, zp
    real*8  :: x, x_s, x_t, y, y_s, y_t, xjac
    real*8  :: dpsi_dz, dpsi_dr
    real*8  :: theta_geom, theta_prev, theta_cumulative

    ! Temporary arrays for one surface's tracing data
    real*8, allocatable :: tt_step(:), t2_step(:)

    ! Basis function arrays
    real*8 :: G(4, n_degrees), G_s(4, n_degrees), G_t(4, n_degrees)
    real*8 :: G_st(4, n_degrees), G_ss(4, n_degrees), G_tt(4, n_degrees)
    integer :: kv, iv, kf

    ierr = 0

    allocate(tt_step(NMAXSTEPS), t2_step(NMAXSTEPS))

    ! --- Loop over flux surfaces (outer to inner for OpenMP compatibility) ---
    do k = nPsiN, 1, -1

      psi_n_k = PsiNRange(1) + (PsiNRange(2)-PsiNRange(1)) * real(k-1)/real(nPsiN-1)

      ! --- Find starting point at outer midplane via R-bisection ---
      Rleft  = equil_state%R_axis
      Rright = equil_state%R_midpl(2)
      do
        Rmid = (Rleft + Rright) / 2.d0
        if (Rright - Rleft < 1.d-7) exit

        call find_RZ(node_list, element_list, Rmid, equil_state%Z_axis, &
             Rout, Zout, ielm, s, t, ifail)
        if (ifail /= 0) exit

        call interp0_seed(node_list, element_list, ielm, 1, 1, s, t, &
             G, G_s, G_t, G_st, G_ss, G_tt, P, P_s, P_t)

        if (get_psi_n(P, Zout) < psi_n_k) then
          Rleft = Rmid
        else
          Rright = Rmid
        end if
      end do

      if (ifail /= 0) then
        write(*,*) 'WARNING: seed_trace_fieldlines bisection failed at surface', k
        cycle
      end if

      ! Starting position
      rn = Rmid
      zn = equil_state%Z_axis
      tt_step(1) = 0.d0
      step_count = 1
      theta_cumulative = 0.d0
      theta_prev = 0.d0

      ! --- Field line tracing: predictor-corrector in phi ---
      do j = 1, NMAXSTEPS

        ! -- Predictor half-step --
        call find_RZ(node_list, element_list, rn, zn, Rout, Zout, ielm, s, t, ifail)
        if (ifail /= 0) exit

        call interp0_seed(node_list, element_list, ielm, 1, 1, s, t, &
             G, G_s, G_t, G_st, G_ss, G_tt, P, P_s, P_t)

        call interp_RZ(node_list, element_list, ielm, s, t, x, x_s, x_t, y, y_s, y_t)

        xjac     = x_s*y_t - x_t*y_s
        dpsi_dz  = (-x_t*P_s + x_s*P_t) / xjac
        dpsi_dr  = ( y_t*P_s - y_s*P_t) / xjac

        rh = rn - dpsi_dz / F0 * rn * DELTAPHI / 2.d0
        zh = zn + dpsi_dr / F0 * rn * DELTAPHI / 2.d0

        ! -- Corrector full-step --
        call find_RZ(node_list, element_list, rh, zh, Rout, Zout, ielm, s, t, ifail)
        if (ifail /= 0) exit

        call interp0_seed(node_list, element_list, ielm, 1, 1, s, t, &
             G, G_s, G_t, G_st, G_ss, G_tt, P, P_s, P_t)

        call interp_RZ(node_list, element_list, ielm, s, t, x, x_s, x_t, y, y_s, y_t)

        xjac     = x_s*y_t - x_t*y_s
        dpsi_dz  = (-x_t*P_s + x_s*P_t) / xjac
        dpsi_dr  = ( y_t*P_s - y_s*P_t) / xjac

        rp = rn - dpsi_dz / F0 * rh * DELTAPHI
        zp = zn + dpsi_dr / F0 * rh * DELTAPHI

        step_count = step_count + 1

        ! Compute geometric poloidal angle at new position
        theta_geom = atan3_seed(zp - equil_state%Z_axis, rp - equil_state%R_axis)

        ! Handle theta wrapping across 2*pi
        if (theta_prev > 3.d0*PI/2.d0 .and. theta_geom < PI/2.d0) then
          theta_cumulative = theta_cumulative + 2.d0*PI
        end if
        if (theta_geom > 3.d0*PI/2.d0 .and. theta_prev < PI/2.d0) then
          theta_cumulative = theta_cumulative - 2.d0*PI
        end if

        tt_step(step_count) = theta_geom + theta_cumulative
        theta_prev = theta_geom

        ! Check if one full poloidal turn is complete
        if (abs(tt_step(step_count) - tt_step(1)) > 2.d0*PI) exit

        rn = rp
        zn = zp
      end do

      npts = step_count

      if (ifail /= 0 .or. npts < 2 .or. npts >= NMAXSTEPS) then
        write(*,*) 'WARNING: seed_trace_fieldlines incomplete at surface', &
             k, ' npts=', npts
        ! Fallback: fill with identity mapping
        do j = 1, nTht
          theta_star_table(j, k) = 2.d0*PI * real(j-1)/real(nTht-1)
        end do
        cycle
      end if

      ! --- Build theta_star: linearly proportional to step index ---
      do j = 1, npts
        t2_step(j) = 2.d0*PI * real(j-1) / real(npts - 1)
      end do

      ! --- Interpolate to equidistant theta_geom grid -> fill table ---
      do j = 1, nTht
        theta_geom = 2.d0*PI * real(j-1) / real(nTht-1)
        call interp1_seed(tt_step, t2_step, npts, theta_geom, theta_star_table(j,k))
      end do

    end do

    deallocate(tt_step, t2_step)

  end subroutine seed_trace_fieldlines


  !=============================================================================
  !> Interpolate scalar variable at (s,t) inside a Bezier element.
  !! Simplified version of mod_straight_field_line::interp0.
  !! Uses pre-computed basis functions (passed in to avoid reallocation).
  !=============================================================================
  subroutine interp0_seed(node_list, element_list, i_elm, i_var, i_harm, s, t, &
       G, G_s, G_t, G_st, G_ss, G_tt, P, P_s, P_t)
    use data_structure, only: type_node_list, type_element_list
    use mod_parameters, only: n_vertex_max, n_degrees
    use mod_basisfunctions, only: basisfunctions

    implicit none
    type(type_node_list),    intent(in)  :: node_list
    type(type_element_list), intent(in)  :: element_list
    integer,                 intent(in)  :: i_elm, i_var, i_harm
    real*8,                  intent(in)  :: s, t
    real*8, intent(out) :: G(4, n_degrees), G_s(4, n_degrees), G_t(4, n_degrees)
    real*8, intent(out) :: G_st(4, n_degrees), G_ss(4, n_degrees), G_tt(4, n_degrees)
    real*8, intent(out) :: P, P_s, P_t
    integer :: kv, iv, kf

    call basisfunctions(s, t, G(1:4,1:n_degrees), G_s(1:4,1:n_degrees), &
         G_t(1:4,1:n_degrees), G_st(1:4,1:n_degrees), G_ss(1:4,1:n_degrees), &
         G_tt(1:4,1:n_degrees))

    P   = 0.d0
    P_s = 0.d0
    P_t = 0.d0

    do kv = 1, n_vertex_max
      iv = element_list%element(i_elm)%vertex(kv)
      do kf = 1, n_degrees
        P   = P   + node_list%node(iv)%values(i_harm, kf, i_var) &
                  * element_list%element(i_elm)%size(kv, kf) * G(kv, kf)
        P_s = P_s + node_list%node(iv)%values(i_harm, kf, i_var) &
                  * element_list%element(i_elm)%size(kv, kf) * G_s(kv, kf)
        P_t = P_t + node_list%node(iv)%values(i_harm, kf, i_var) &
                  * element_list%element(i_elm)%size(kv, kf) * G_t(kv, kf)
      end do
    end do

  end subroutine interp0_seed


  !=============================================================================
  !> Linear interpolation via binary search: y(x_target).
  !! Assumes x(1:n) is strictly monotonic.
  !=============================================================================
  subroutine interp1_seed(x, y, n, x_target, y_interp)
    implicit none
    real*8,  intent(in)  :: x(n), y(n)
    integer, intent(in)  :: n
    real*8,  intent(in)  :: x_target
    real*8,  intent(out) :: y_interp

    integer :: lo, hi, mid

    lo = 1
    hi = n

    if (x_target <= x(1)) then
      y_interp = y(1)
      return
    end if
    if (x_target >= x(n)) then
      y_interp = y(n)
      return
    end if

    do while (hi - lo > 1)
      mid = (lo + hi) / 2
      if (x(mid) <= x_target) then
        lo = mid
      else
        hi = mid
      end if
    end do

    if (abs(x(hi) - x(lo)) < 1.d-14) then
      y_interp = y(lo)
    else
      y_interp = y(lo) + (x_target - x(lo)) / (x(hi) - x(lo)) * (y(hi) - y(lo))
    end if

  end subroutine interp1_seed


  !=============================================================================
  !> Like atan2 but always returns values in [0, 2*pi).
  !=============================================================================
  pure real*8 function atan3_seed(dy, dx)
    implicit none
    real*8, intent(in) :: dy, dx
    atan3_seed = atan2(dy, dx)
    if (atan3_seed < 0.d0) atan3_seed = atan3_seed + 2.d0*PI
  end function atan3_seed

end module mod_seed_theta_lookup
