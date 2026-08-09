!> Program to convert a JOREK2 restart file into 3D binary VTK format
program jorek2vtk_3d

use mod_parameters, only: n_order, n_var, variable_names
use mod_chi
use constants
use data_structure
use phys_module
use equil_info, only: get_psi_n
use mod_import_restart
use mod_interp
use basis_at_gaussian, only: initialise_basis
implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list

integer               :: nnoel, nnos, nel, nsub, inode, ielm, n_scalars, n_vectors
real*4,allocatable    :: xyz (:,:), scalars(:,:), vectors(:,:,:)
real*8,allocatable    :: HZ(:,:), HZ_p(:,:), HZ_coord(:,:), HZ_coord_p(:,:)
integer,allocatable   :: ien (:,:)
integer, parameter    :: ivtk = 22
integer               :: i, j, k, m, etype, irst, int, i_var, i_tor, index, index_node, n_points, k_tor
integer               :: n_toroidal
character             :: buffer*80, lf*1, str1*10, str2*10
character*36, allocatable :: scalar_names(:), vector_names(:)
real*4                :: float
real*8                :: s, t, phi, angle, cur_pert
real*8                :: P,P_s,P_t,P_st,P_ss,P_tt
real*8                :: R,R_s,R_t,R_phi,R_st,R_ss,R_tt,R_sp,R_tp,R_pp
real*8                :: Z,Z_s,Z_t,Z_p,Z_st,Z_ss,Z_tt,Z_sp,Z_tp,Z_pp
real*8                :: Psi,Ps_s,Ps_t,Ps_st,Ps_ss,Ps_tt, ZJ,ZJ_s,ZJ_t,ZJ_st,ZJ_ss,ZJ_tt, W,W_s,W_t,W_st,W_ss,W_tt
real*8                :: ps_p, ps_x_itor, ps_y_itor
real*8                :: U,U_s,U_t,U_st,U_ss,U_tt, RHO,RH_s,RH_t,RH_st,RH_ss,RH_tt, TT,TT_s,TT_t,TT_st,TT_ss,TT_tt
real*8                :: u_x, u_y, u_p
real*8                :: u0_x, u0_y, xjac, xjac_s, xjac_t, xjac_x, xjac_y, v_perp, Psi_J, R_p, error, zj_x, zj_y, ps_x, ps_y
real*8                :: Bx, By, Bz
real*8                :: V, Vx, Vy, Vz
real*8                :: grad_chi(3), Bv2
real*8, dimension(0:n_order-1,0:n_order-1,0:n_order-1) :: chi

integer               :: i_vec_x, i_vec_y, i_vec_z
real*8                :: angle_sign
real*8                :: rho_norm, t_norm    ! for SI unit conversion

! --- Booleans
logical               :: periodic, density_only, SI_units
integer               :: ierr, my_id
logical               :: without_n0_mode, RphiZ_coords

! --- Additional diagnostic flags
logical               :: include_pressure    ! output pressure (rho * T)
logical               :: include_psi_norm    ! output normalized poloidal flux
logical               :: include_B_scalars   ! output B_R, B_Z, B_phi, B_abs as scalars
logical               :: include_E_field     ! output E-field as a vector

! --- Scalar / vector index bookkeeping
integer               :: idx_psi, idx_U, idx_j, idx_omega, idx_rho, idx_T, idx_Vpar
integer               :: idx_pressure, idx_psin
integer               :: idx_BR, idx_BZ, idx_Bphi, idx_Babs
integer               :: idx_Ti, idx_Te
integer               :: i_vec_B, i_vec_V, i_vec_E

! --- Physics helpers
real*8                :: Br_val, Bz_val, Bp_val, Babs_val
real*8                :: E_R, E_Z, E_phi
real*8                :: psi_tmp, Z_tmp

namelist /vtk_params/ nsub, without_n0_mode, periodic, RphiZ_coords, density_only, &
                      n_toroidal, SI_units, &
                      include_pressure, include_psi_norm, include_B_scalars, include_E_field

write(*,*) 'jorek2vtk_3d'

! --- Initialise input parameters and read the input namelist.
my_id     = 0
call initialise_parameters(my_id, "__NO_FILENAME__")
call initialise_basis()

! --- Preset parameters
nsub               = 5
without_n0_mode    = .false.
periodic           = .true.
density_only       = .false.
SI_units           = .false.
n_toroidal         = 200
RphiZ_coords       = .true.
include_pressure   = .false.
include_psi_norm   = .false.
include_B_scalars  = .false.
include_E_field    = .false.

! --- Read parameters from namelist file 'vtk.nml' if it exists
open(42, file='vtk.nml', action='read', status='old', iostat=ierr)
if ( ierr == 0 ) then
  write(*,*) 'Reading parameters from vtk.nml namelist.'
  read(42,vtk_params)
  close(42)
end if

write(*,*)
write(*,*) 'Parameters:'
write(*,*) '-----------'
write(*,*) 'nsub              =', nsub
write(*,*) 'without_n0_mode   =', without_n0_mode
write(*,*) 'periodic          =', periodic
write(*,*) 'SI_units          =', SI_units
write(*,*) 'include_pressure  =', include_pressure
write(*,*) 'include_psi_norm  =', include_psi_norm
write(*,*) 'include_B_scalars =', include_B_scalars
write(*,*) 'include_E_field   =', include_E_field
write(*,*)

! --- Build scalar / vector index table and compute total counts
if (density_only) then
  idx_rho   = 1
  n_scalars = 1
  n_vectors = 0
else
  ! Base 7 scalars
  idx_psi   = 1
  idx_U     = 2
  idx_j     = 3
  idx_omega = 4
  idx_rho   = 5
  idx_T     = 6
  idx_Vpar  = 7
  n_scalars = 7

  if (include_pressure) then
    n_scalars = n_scalars + 1
    idx_pressure = n_scalars
  end if

  if (include_psi_norm) then
    n_scalars = n_scalars + 1
    idx_psin = n_scalars
  end if

  if (include_B_scalars) then
    n_scalars = n_scalars + 1; idx_BR   = n_scalars
    n_scalars = n_scalars + 1; idx_BZ   = n_scalars
    n_scalars = n_scalars + 1; idx_Bphi = n_scalars
    n_scalars = n_scalars + 1; idx_Babs = n_scalars
  end if

  ! Optional Ti / Te (only when with_TiTe)
  if (with_TiTe) then
    n_scalars = n_scalars + 1; idx_Ti = n_scalars
    n_scalars = n_scalars + 1; idx_Te = n_scalars
  end if

  ! Base 2 vectors: B_field, v_field
  i_vec_B = 1
  i_vec_V = 2
  n_vectors = 2

  if (include_E_field) then
    n_vectors = n_vectors + 1
    i_vec_E = n_vectors
  end if
end if

do i_tor=1, n_tor
  mode(i_tor) = + int(i_tor / 2) * n_period
enddo

do k_tor=1, n_coord_tor
  mode_coord(k_tor) = + int(k_tor / 2) * n_coord_period
enddo

call init_chi_basis

call import_restart(node_list, element_list, 'jorek_restart', rst_format, ierr, .true.)
nnos = n_toroidal * nsub*nsub*node_list%n_nodes

! --- Allocate and name scalars / vectors
allocate(xyz(3,nnos), scalars(nnos,1:n_scalars), scalar_names(n_scalars))

if (density_only) then
  if (SI_units) then
    scalar_names(1) = 'density_1e20m-3'
  else
    scalar_names(1) = 'density'
  end if
else
  allocate(vector_names(n_vectors), vectors(nnos,3,1:n_vectors))

  if (SI_units) then
    scalar_names(idx_psi)   = 'flux'
    scalar_names(idx_U)     = 'U_m/s'
    scalar_names(idx_j)     = 'j_MA/m2'
    scalar_names(idx_omega) = 'omega'
    scalar_names(idx_rho)   = 'density_1e20m-3'
    scalar_names(idx_T)     = 'T_keV'
    scalar_names(idx_Vpar)  = 'v_par_km/s'
  else
    scalar_names(idx_psi)   = 'flux'
    scalar_names(idx_U)     = 'U'
    scalar_names(idx_j)     = 'j'
    scalar_names(idx_omega) = 'omega'
    scalar_names(idx_rho)   = 'density'
    scalar_names(idx_T)     = 'T'
    scalar_names(idx_Vpar)  = 'v_par'
  end if

  vector_names(i_vec_B) = 'B_field'
  vector_names(i_vec_V) = 'v_field'

  if (include_pressure)   scalar_names(idx_pressure) = 'pressure'
  if (include_psi_norm)   scalar_names(idx_psin)     = 'psi_norm'
  if (include_B_scalars) then
    scalar_names(idx_BR)   = 'B_R'
    scalar_names(idx_BZ)   = 'B_Z'
    scalar_names(idx_Bphi) = 'B_phi'
    scalar_names(idx_Babs) = 'B_abs'
  end if
  if (with_TiTe) then
    scalar_names(idx_Ti) = 'Ti'
    scalar_names(idx_Te) = 'Te'
  end if
  if (include_E_field) vector_names(i_vec_E) = 'E_field'
end if

if (periodic) then
  nel   = (n_toroidal)   * (nsub-1)*(nsub-1)*element_list%n_elements
else
  nel   = (n_toroidal-1) * (nsub-1)*(nsub-1)*element_list%n_elements
end if

nnoel = 8
allocate(ien(nnoel,nel))

inode   = 0
ielm    = 0
scalars = 0.d0
vectors = 0.d0
xyz     = 0
ien     = 0
n_points = nsub*nsub*element_list%n_elements

allocate(HZ(n_tor,n_toroidal))
allocate(HZ_p(n_tor,n_toroidal))
allocate(HZ_coord(n_coord_tor,n_toroidal))
allocate(HZ_coord_p(n_coord_tor,n_toroidal))

do m=1,n_toroidal
  if (periodic) then
    phi = 2.d0 * PI * float(m-1)/float(n_toroidal)
  else
    phi = 2.d0 * PI * float(m-1)/float(n_toroidal-1) / float(n_period)
  endif
  HZ(1,m)   = 1.d0
  HZ_p(1,m) = 0.d0
  do i=1,(n_tor-1)/2
    HZ(2*i,m)     = cos(mode(2*i)  *phi)
    HZ(2*i+1,m)   = sin(mode(2*i+1)*phi)
    HZ_p(2*i,m)   = -float(mode(2*i))   * sin(mode(2*i)  *phi)
    HZ_p(2*i+1,m) =  float(mode(2*i+1)) * cos(mode(2*i+1)*phi)
  enddo

  HZ_coord(1,m) = 1.0
  HZ_coord_p(1,m) = 0.d0
  do i=1,(n_coord_tor-1)/2
    HZ_coord(2*i,m)      =                           cos(mode_coord(2*i)  *phi)
    HZ_coord_p(2*i,m)    = - float(mode_coord(2*i))      * sin(mode_coord(2*i)  *phi)
    HZ_coord(2*i+1,m)    =                         - sin(mode_coord(2*i+1)*phi)
    HZ_coord_p(2*i+1,m)  = - float(mode_coord(2*i+1))    * cos(mode_coord(2*i+1)*phi)
  enddo
enddo

do m=1, n_toroidal
  if ( mod(m,n_toroidal/40+1) == 0 ) write(*,'(" Plane ",i4.4," of ",i4.4)') m, n_toroidal

  if (periodic) then
    angle = 2.d0 * PI * float(m-1)/float(n_toroidal)
  else
    angle = 2.d0 * PI * float(m-1)/float(n_toroidal-1) / float(n_period)
  endif

  do i=1,element_list%n_elements

    do j=1,nsub
      s = float(j-1)/float(nsub-1)
      do k=1,nsub
        t = float(k-1)/float(nsub-1)

        call interp_RZP(node_list,element_list,i,s,t,angle,R,R_s,R_t,R_phi,R_st,R_ss,R_tt,R_sp,R_tp,R_pp, &
                       Z,Z_s,Z_t,Z_p,Z_st,Z_ss,Z_tt,Z_sp,Z_tp,Z_pp)

        xjac  = R_s * Z_t - R_t * Z_s
        xjac_s = R_ss*Z_t + R_s*Z_st - R_st*Z_s - R_t*Z_ss
        xjac_t = R_st*Z_t + R_s*Z_tt - R_tt*Z_s - R_t*Z_st
        xjac_x = (xjac_s*Z_t - xjac_t*Z_s) / xjac
        xjac_y = (R_s*xjac_t - R_t*xjac_s) / xjac

        chi = get_chi(R,Z,angle,node_list,element_list,i,s,t)
        grad_chi = (/ chi(1,0,0), chi(0,1,0), chi(0,0,1)/R /)
        Bv2 = dot_product(grad_chi,grad_chi)

        inode = inode+1

        if ( xjac == 0.d0 ) xjac = 1.d-8

        i_vec_x = 1
        if (RphiZ_coords) then
          xyz(1:3,inode) = (/ R * cos(angle), -R*sin(angle), Z /)
          i_vec_Z = 3
          i_vec_y = 2
          angle_sign = -1
        else
          xyz(1:3,inode) = (/ R * cos(angle), Z, R*sin(angle) /)
          i_vec_Z = 2
          i_vec_y = 3
          angle_sign = 1
        endif

        ps_x = 0.d0; ps_y = 0.d0; ps_p = 0.d0
        u_x = 0.d0; u_y = 0.d0; u_p = 0.d0
        do i_tor = 1,n_tor

          if ( ( i_tor == 1 ) .and. ( without_n0_mode ) ) cycle

          if(density_only) then
            call interp(node_list,element_list,i,var_rho,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            scalars(inode,idx_rho) = scalars(inode,idx_rho) + P * HZ(i_tor,m)
          else
            call interp(node_list,element_list,i,var_psi,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            scalars(inode,idx_psi) = scalars(inode,idx_psi) + P * HZ(i_tor,m)

            ps_x  = ps_x + (  Z_t * P_s - Z_s * P_t  ) / xjac * HZ(i_tor,m)
            ps_y  = ps_y + ( - R_t * P_s + R_s * P_t ) / xjac * HZ(i_tor,m)
                ps_p = ps_p + P*HZ_p(i_tor,m) - (   Z_t * P_s - Z_s * P_t ) / xjac * HZ(i_tor,m)*R_phi &
                                              - ( - R_t * P_s + R_s * P_t ) / xjac * HZ(i_tor,m)*Z_p

            call interp(node_list,element_list,i,var_u,i_tor,s,t,U,U_s,U_t,U_st,U_ss,U_tt)
            scalars(inode,idx_U) = scalars(inode,idx_U) + U * HZ(i_tor,m)

            u_x  = u_x   + (   Z_t * U_s - Z_s * U_t )     / xjac * HZ(i_tor,m)
            u_y  = u_y   + ( - R_t * U_s + R_s * U_t )     / xjac * HZ(i_tor,m)
            u_p = u_p + U*HZ_p(i_tor,m) - (   Z_t * U_s - Z_s * U_t ) / xjac * HZ(i_tor,m)*R_phi &
                                        - ( - R_t * U_s + R_s * U_t ) / xjac * HZ(i_tor,m)*Z_p

            call interp(node_list,element_list,i,var_zj,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            scalars(inode,idx_j) = scalars(inode,idx_j) + P * HZ(i_tor,m)

            call interp(node_list,element_list,i,var_w,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            scalars(inode,idx_omega) = scalars(inode,idx_omega) + P * HZ(i_tor,m)

            call interp(node_list,element_list,i,var_rho,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            scalars(inode,idx_rho) = scalars(inode,idx_rho) + P * HZ(i_tor,m)

            call interp(node_list,element_list,i,var_T,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            scalars(inode,idx_T) = scalars(inode,idx_T) + P * HZ(i_tor,m)

            if (with_Vpar) then
              call interp(node_list,element_list,i,var_Vpar,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            else
              P = 0.d0
            end if
            scalars(inode,idx_Vpar) = scalars(inode,idx_Vpar) + P * HZ(i_tor,m)

            ! --- Optional Ti / Te (accumulate over all toroidal harmonics)
            if (with_TiTe) then
              call interp(node_list,element_list,i,var_Ti,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
              scalars(inode,idx_Ti) = scalars(inode,idx_Ti) + P * HZ(i_tor,m)
              call interp(node_list,element_list,i,var_Te,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
              scalars(inode,idx_Te) = scalars(inode,idx_Te) + P * HZ(i_tor,m)
            end if
          endif
        enddo

        if(.not. density_only) then

          ! --- B-field local components (for scalars / E-field)
          Br_val  = ps_y / R
          Bz_val  = -ps_x / R
          Bp_val  = F0 / R
          Babs_val = sqrt(ps_x**2 + ps_y**2 + F0**2) / R

          ! --- Pressure (rho * T)
          if (include_pressure) then
            scalars(inode,idx_pressure) = scalars(inode,idx_rho) * scalars(inode,idx_T)
          end if

          ! --- Normalized flux
          if (include_psi_norm) then
            psi_tmp = scalars(inode,idx_psi)
            Z_tmp   = Z
            scalars(inode,idx_psin) = real(get_psi_n(psi_tmp, Z_tmp))
          end if

          ! --- B-field scalars
          if (include_B_scalars) then
            scalars(inode,idx_BR)   = Br_val
            scalars(inode,idx_BZ)   = Bz_val
            scalars(inode,idx_Bphi) = Bp_val
            scalars(inode,idx_Babs) = Babs_val
          end if

          ! --- B_field and v_field vectors
          if (     (jorek_model .eq. 180).or. (jorek_model .eq. 183) )then
              Bx = chi(1,0,0)      + (ps_y*chi(0,0,1) - ps_p*chi(0,1,0))/(F0*R)
              By = chi(0,1,0)      - (ps_x*chi(0,0,1) - ps_p*chi(1,0,0))/(F0*R)
              Bz = chi(0,0,1)/R    + (ps_x*chi(0,1,0) - ps_y*chi(1,0,0))/F0
              vectors(inode,i_vec_x, i_vec_B) =              Bx * cos(angle) - angle_sign * Bz * sin(angle)
              vectors(inode,i_vec_z, i_vec_B) = By
              vectors(inode,i_vec_y, i_vec_B) = angle_sign * Bx * sin(angle) +              Bz * cos(angle)

              Vx =                   ( u_y*chi(0,0,1) - u_p*chi(0,1,0))/(R*Bv2)
              Vy =                   (-u_x*chi(0,0,1) + u_p*chi(1,0,0))/(R*Bv2)
              Vz =                   ( u_x*chi(0,1,0) - u_y*chi(1,0,0))/Bv2
              vectors(inode,i_vec_x, i_vec_V) =               Vx * cos(angle) - angle_sign * Vz * sin(angle)
              vectors(inode,i_vec_z, i_vec_V) =  Vy
              vectors(inode,i_vec_y, i_vec_V) =  angle_sign * Vx * sin(angle) +              Vz * cos(angle)
          else
              vectors(inode,i_vec_x, i_vec_B) = - angle_sign * F0/R * sin(angle) +              ps_y / R * cos(angle)
              vectors(inode,i_vec_z, i_vec_B) = - ps_x / R
              vectors(inode,i_vec_y, i_vec_B) =                F0/R * cos(angle) + angle_sign * ps_y / R * sin(angle)

              V = scalars(inode,idx_Vpar)
              vectors(inode,i_vec_x, i_vec_V) = (V/R*ps_y - R * u_y) * cos(angle) - angle_sign * V*F0/R               * sin(angle)
              vectors(inode,i_vec_z, i_vec_V) = -V/R*ps_x + R * u_x
              vectors(inode,i_vec_y, i_vec_V) =  V*F0/R              * cos(angle) + angle_sign * (V/R*ps_y - R * u_y) * sin(angle)
          endif

          ! --- E-field vector (E = -grad(u) related)
          if (include_E_field) then
            E_R = -u_x
            E_Z = -u_y
            E_phi = -u_p / R
            vectors(inode,i_vec_x, i_vec_E) =               E_R * cos(angle) - angle_sign * E_phi * sin(angle)
            vectors(inode,i_vec_z, i_vec_E) = E_Z
            vectors(inode,i_vec_y, i_vec_E) = angle_sign * E_R * sin(angle) +              E_phi * cos(angle)
          end if

        end if
      enddo
    enddo

    if (m .lt. n_toroidal) then
      do j=1,nsub-1
        do k=1,nsub-1
          ielm        = ielm+1
          ien(1,ielm) = inode - nsub*nsub + nsub*(j-1) + k-1
          ien(2,ielm) = inode - nsub*nsub + nsub*(j  ) + k-1
          ien(3,ielm) = inode - nsub*nsub + nsub*(j  ) + k
          ien(4,ielm) = inode - nsub*nsub + nsub*(j-1) + k
          ien(5,ielm) = ien(1,ielm) + n_points
          ien(6,ielm) = ien(2,ielm) + n_points
          ien(7,ielm) = ien(3,ielm) + n_points
          ien(8,ielm) = ien(4,ielm) + n_points
        enddo
      enddo
    endif

    if ( (periodic) .and. (m .eq. n_toroidal)) then
      do j=1,nsub-1
        do k=1,nsub-1
          ielm        = ielm+1
          ien(1,ielm) = inode - nsub*nsub + nsub*(j-1) + k-1
          ien(2,ielm) = inode - nsub*nsub + nsub*(j  ) + k-1
          ien(3,ielm) = inode - nsub*nsub + nsub*(j  ) + k
          ien(4,ielm) = inode - nsub*nsub + nsub*(j-1) + k
          ien(5,ielm) = ien(1,ielm) - n_points * (n_toroidal-1)
          ien(6,ielm) = ien(2,ielm) - n_points * (n_toroidal-1)
          ien(7,ielm) = ien(3,ielm) - n_points * (n_toroidal-1)
          ien(8,ielm) = ien(4,ielm) - n_points * (n_toroidal-1)
        enddo
      enddo
    endif

  enddo
enddo

!--------------------------------------------------- SI unit conversion
if (SI_units) then

  rho_norm = central_density * 1.d20 * central_mass * ATOMIC_MASS_UNIT
  t_norm   = sqrt(MU_zero * rho_norm)

  do i = 1, nnos

    if (density_only) then
      scalars(i, idx_rho) = scalars(i, idx_rho) * central_density
    else
      ! Base variables
      scalars(i, idx_U)    = scalars(i, idx_U)    / t_norm                           ! m/s
      scalars(i, idx_j)    = scalars(i, idx_j)    / MU_zero * 1.e-6                  ! MA/m^2
      scalars(i, idx_rho)  = scalars(i, idx_rho)  * central_density                  ! 10^20 m^-3
      scalars(i, idx_T)    = scalars(i, idx_T)    / MU_zero / (central_density * 1.d20) / EL_CHG / 2.d0 / 1.e3  ! keV
      scalars(i, idx_Vpar) = scalars(i, idx_Vpar) / t_norm / 1.e3                    ! km/s

      ! vectors
      vectors(i, :, i_vec_V) = vectors(i, :, i_vec_V) / t_norm                       ! m/s

      ! Pressure: JOREK -> kPa
      if (include_pressure) then
        scalars(i, idx_pressure) = scalars(i, idx_pressure) / MU_zero / 1.e3
      end if

      ! Ti/Te: JOREK -> keV (separate electron/ion, factor 1 instead of 1/2)
      if (with_TiTe) then
        scalars(i, idx_Ti) = scalars(i, idx_Ti) / MU_zero / (central_density * 1.d20) / EL_CHG / 1.e3
        scalars(i, idx_Te) = scalars(i, idx_Te) / MU_zero / (central_density * 1.d20) / EL_CHG / 1.e3
      end if

      ! E-field vector: JOREK -> kV/m
      if (include_E_field) then
        vectors(i, :, i_vec_E) = vectors(i, :, i_vec_E) / t_norm / 1.e3
      end if
    endif

  enddo

  write(*,*) 'Applied SI unit conversion.'
endif

!--------------------------------------------------- write the binary VTK file
etype = 12

lf = char(10)

#ifdef IBM_MACHINE
open(unit=ivtk,file='jorek_tmp.vtk',form='unformatted',access='stream',status='replace')
#else
open(unit=ivtk,file='jorek_tmp.vtk',form='unformatted',access='stream',convert='BIG_ENDIAN',status='replace')
#endif

buffer = '# vtk DataFile Version 3.0'//lf                                             ; write(ivtk) trim(buffer)
buffer = 'vtk output'//lf                                                             ; write(ivtk) trim(buffer)
buffer = 'BINARY'//lf                                                                 ; write(ivtk) trim(buffer)
buffer = 'DATASET UNSTRUCTURED_GRID'//lf//lf                                          ; write(ivtk) trim(buffer)

! POINTS SECTION
write(str1(1:10),'(i10)') nnos
buffer = 'POINTS '//str1//'  float'//lf                                               ; write(ivtk) trim(buffer)
write(ivtk) ((xyz(i,j),i=1,3),j=1,nnos)

! CELLS SECTION
write(str1(1:10),'(i10)') nel
write(str2(1:10),'(i10)') nel*(1+nnoel)
buffer = lf//lf//'CELLS '//str1//' '//str2//lf                                        ; write(ivtk) trim(buffer)
write(ivtk) (nnoel,(ien(i,j),i=1,nnoel),j=1,nel)

! CELL_TYPES SECTION
write(str1(1:10),'(i10)') nel
buffer = lf//lf//'CELL_TYPES'//str1//lf                                               ; write(ivtk) trim(buffer)
write(ivtk) (etype,i=1,nel)

! POINT_DATA SECTION
write(str1(1:10),'(i10)') nnos
buffer = lf//lf//'POINT_DATA '//str1//lf                                              ; write(ivtk) trim(buffer)

do i_var =1, n_scalars
  buffer = 'SCALARS '//scalar_names(i_var)//' float'//lf                              ; write(ivtk) trim(buffer)
  buffer = 'LOOKUP_TABLE default'//lf                                                 ; write(ivtk) trim(buffer)
  write(ivtk) (scalars(i,i_var),i=1,nnos)
enddo

if(.not. density_only) then
  do i_var =1, n_vectors
    buffer = lf//lf//'VECTORS '//vector_names(i_var)//' float'//lf                    ; write(ivtk) trim(buffer)
    write(ivtk) ((vectors(j,i,i_var),i=1,3),j=1,nnos)
  enddo
endif

close(ivtk)

write(*,*) 'done.'

end program jorek2vtk_3d
