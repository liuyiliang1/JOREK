!**********************************************************************
!* program to convert a JOREK2 restart file into binary VTK format    *
!**********************************************************************

program jorek_diagnostics
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use mod_parameters, only: n_var
use data_structure
use phys_module
use basis_at_gaussian
use diffusivities, only: get_dperp, get_zkperp
use nodes_elements
use mod_boundary
use mod_import_restart
use equil_info, only : get_psi_n, ES
use mod_interp

implicit none

integer               :: nnoel, nnos, nel, nsub, inode, ielm, n_scalars, n_vectors, my_id
real*4,allocatable    :: xyz (:,:), scalars(:,:), vectors(:,:,:)
integer,allocatable   :: ien (:,:)
integer               :: i, j, k, m, etype, ivtk, irst, int, i_var, i_tor, index, index_node, n_points
integer               :: n_bnd, iv, iv1, iv2, inode1, inode2, ierr, i_elm_xpoint(2), ifail, i_elm_axis, inode_avg
character             :: buffer*80, lf*1, str1*8, str2*8
character*8, allocatable :: scalar_names(:), vector_names(:)
real*4                :: float
real*8                :: sg, tg, phi, angle
real*8,allocatable    :: avg6(:), avg7(:), avg8(:), avg11(:), Ravg(:), Zavg(:),B_field(:,:),yield(:)
real*8,allocatable    :: prf1(:), prf2(:), prf3(:), prf4(:), prf5(:), prf6(:), prf7(:), prf8(:),prf10(:), prf11(:), prf14(:), Rprf(:), Zprf(:),prf15(:)
  integer,allocatable   :: prf_bnd(:)  ! boundary type of each point
logical,allocatable   :: tobedone(:)
real*8                :: BigR, xjac, psi_norm
real*8                :: R,R_s,R_t,R_st,R_ss,R_tt, Z,Z_s,Z_t,Z_st,Z_ss,Z_tt
real*8                :: PS,PS_s,PS_t,PS_st,PS_ss,PS_tt, VP,VP_s,VP_t,VP_st,VP_ss,VP_tt
real*8                :: RH,RH_s,RH_t,RH_st,RH_ss,RH_tt, TT,TT_s,TT_t,TT_st,TT_ss,TT_tt
real*8                :: rho, rho_s, rho_t, T, T_s, T_t, T_p, vpar, D_prof, ZK_prof, ZKpar_T,rhon,rhon_s,rhon_t
real*8                :: Ti,Ti_s,Ti_t,Ti_p,Te,Te_s,Te_t,Te_p,TTi,TTi_s,TTi_t,TTi_ss,TTi_tt,TTi_st,TTe,TTe_s,TTe_t,TTe_ss,TTe_tt,TTe_st
real*8                :: psi, psi_s, psi_t, psi_x, psi_y, T_x, T_y, BB2, normal, rho_norm, t_norm,rhon_norm
real*8                :: distance, R_start, Z_start 
logical               :: periodic
logical               :: without_n0_mode
integer               :: i_bnd_node, i_node
  integer, allocatable   :: order_idx(:)  ! ordered indices from nearest-neighbor walk
  integer :: idx, idx1, idx2, bnd, nxt_bnd, seg_start, seg_end, jdx, jdx1, jdx2, jdx_j, jdx_k, prev_idx, nxt_idx, kdx1, kdx2
  integer, allocatable :: target_ord_idx(:), plate_idx(:)
  real*8, allocatable :: psi2d(:,:), s_arc(:)  ! 2D surface output: psi per (point, toroidal plane), arc length
  integer :: n_target_pts, plate_id, p_start, p_end, pdx, n_plate_pts
  real*8  :: dist_min_strike_s, strike_dist_sorted
  integer :: m2d, i2d, idx2d, inode2d
  real*8  :: psi_n2d, rho2d, phi_deg2d, s_cum2d
  logical :: strike_found_sorted
  real*8  :: strike_dist_at, strike_R, strike_Z, strike_psi, cumdist, sdist, dist_to_strike, dgap
  logical :: is_last, is_gap
  real*8  :: frac, psi1, psi2
  integer :: bnd1, bnd2, assigned_k
  logical :: strike_used(10)  ! max 10 strike points
  logical :: strike_found
  integer :: ibndelm
  real*8  :: dist_min, dist_min_strike, dist_ps
  real*8, parameter :: jump_threshold = 0.3d0  ! m, > this => new plate

namelist /vtk_params/ nsub, periodic, without_n0_mode

write(*,*) '*********************************'
write(*,*) '* jorek_target2vtk              *'
write(*,*) '*********************************'

! --- Initialise input parameters and read the input namelist.
my_id = 0
call initialise_parameters(my_id, "__NO_FILENAME__")

! --- Preset parameters
nsub      = 5             ! Number of subdivisions of the cubic finite elements into linear pieces
periodic  = .false.
without_n0_mode = .false. ! If .true., the axisymmetric part will not be included
gamma_sheath = gamma_i_stangeby
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
write(*,*) 'nsub            =', nsub
write(*,*) 'periodic        =', periodic
write(*,*) 'without_n0_mode =', without_n0_mode
write(*,*) '-----------'
write(*,*) 'n_tor           =', n_tor
write(*,*) 'n_period        =', n_period
write(*,*) 'F0        =', F0
write(*,*)
call flush_it(6)

ivtk = 21                 ! an arbitrary unit number for the VTK output file

n_scalars = 15             ! number of scalars to write to the VTK output file
n_vectors = 3

allocate(scalar_names(n_scalars), vector_names(n_vectors))

scalar_names = (/ 'flux    ','density ','T       ','Vpar    ','nV.n    ','nTV.n   ','KparT.n ', &
                  'Kperp.T ','Dperp.n ','B.n     ','nvT_gam ','n       ','nV3.n   ','Te      ','sputter'/)

vector_names = (/ 'B_field ','Velocity','normal  '/)

do i_tor=1, n_tor
  mode(i_tor) = + int(i_tor / 2) * n_period
enddo

call import_restart(node_list,element_list, 'jorek_restart', rst_format, ierr, .true.)

! Fix ES with correct values from HDF5 time series (update_equil_state corrupts them)
ES%psi_bnd = psi_bnd_t(index_start)
do i = 1, 2
  ES%psi_xpoint(i) = psi_xpoint_t(index_start, i)
  ES%R_xpoint(i)   = R_xpoint_t(index_start, i)
  ES%Z_xpoint(i)   = Z_xpoint_t(index_start, i)
end do

call initialise_basis                              ! define the basis functions at the Gaussian points

rho_norm = central_density*1.d20 * central_mass * 1.67d-27
rhon_norm = rho_norm
t_norm   = sqrt(MU_zero*rho_norm)

! --- Find the lowest point on the outer divertor target (required later)
call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)
R_start = 1.d99
Z_start = 1.d99
do i_bnd_node = 1, bnd_node_list%n_bnd_nodes
  i_node = bnd_node_list%bnd_node(i_bnd_node)%index_jorek
  R = node_list%node(i_node)%x(1,1,1)
  Z = node_list%node(i_node)%x(1,1,2)
  if ( ( node_list%node(i_node)%boundary == 3 ) .and. ( R > ES%R_xpoint(2) ) .and. ( Z < Z_start ) ) then
    R_start = R
    Z_start = Z
  end if
end do

n_bnd = 0
do i=1,element_list%n_elements  
  do iv = 1, n_vertex_max

    iv1 = iv
    iv2 = mod(iv,4) + 1

    inode1 = element_list%element(i)%vertex(iv1)
    inode2 = element_list%element(i)%vertex(iv2)

    if   (( ((node_list%node(inode1)%boundary .eq. 1) .or. (node_list%node(inode1)%boundary .eq. 3))   &
       .and.  ((node_list%node(inode2)%boundary .eq. 1) .or. (node_list%node(inode2)%boundary .eq. 3)) ) &

       .or. ( ((node_list%node(inode1)%boundary .eq. 19) .or. (node_list%node(inode1)%boundary .eq. 11))   &
       .and.  ((node_list%node(inode2)%boundary .eq. 19) .or. (node_list%node(inode2)%boundary .eq. 11)) ) &

       .or. ( ((node_list%node(inode1)%boundary .eq. 4) .or. (node_list%node(inode1)%boundary .eq. 9))   &
       .and.  ((node_list%node(inode2)%boundary .eq. 4) .or. (node_list%node(inode2)%boundary .eq. 9)) ) &

       .or. ( ((node_list%node(inode1)%boundary .eq. 4) .or. (node_list%node(inode1)%boundary .eq. 1))   &
       .and.  ((node_list%node(inode2)%boundary .eq. 4) .or. (node_list%node(inode2)%boundary .eq. 1)) ) &

       .or. ( ((node_list%node(inode1)%boundary .eq. 5) .or. (node_list%node(inode1)%boundary .eq. 9))   &
       .and.  ((node_list%node(inode2)%boundary .eq. 5) .or. (node_list%node(inode2)%boundary .eq. 9)) )) then 

      n_bnd = n_bnd + 1

    endif
  enddo
enddo

write(*,*) ' number of boundary points : ',n_bnd
write(*,*) ' number of toroidal planes : ',n_plane

nnos = n_plane * nsub * n_bnd
allocate(xyz(3,nnos),scalars(nnos,1:n_scalars),vectors(nnos,3,1:n_vectors))

nnoel = 4

if (periodic) then
  nel   = (n_plane)   * (nsub-1)*n_bnd
else
  nel   = (n_plane-1) * (nsub-1)*n_bnd
endif

allocate(ien(nnoel,nel))

inode   = 0
ielm    = 0
scalars = 0.d0
vectors = 0.d0
xyz     = 0
ien     = 0
n_points = nsub*n_bnd       ! number of points in one poloidal plane

allocate(avg6(n_points),avg7(n_points),avg8(n_points),avg11(n_points),Ravg(n_points),Zavg(n_points),B_field(n_points,4))
avg6 = 0.d0; avg7 = 0.d0; avg8 = 0.d0; avg11 = 0.d0; Ravg = 0.d0; Zavg = 0.d0

allocate(prf1(n_points),prf2(n_points),prf3(n_points),prf4(n_points),prf5(n_points),Rprf(n_points),Zprf(n_points))
prf2 = 0.d0; prf3 = 0.d0; prf4 = 0.d0; prf5 = 0.d0; Rprf = 0.d0; Zprf = 0.d0
allocate(prf6(n_points),prf7(n_points),prf8(n_points),prf10(n_points),prf11(n_points),prf14(n_points),prf15(n_points),yield(n_points))
prf6 = 0.d0; prf7 = 0.d0; prf8 = 0.d0; prf10(:)=0.d0; prf11 = 0.d0; prf14 = 0.d0; prf15 = 0.d0
allocate(prf_bnd(n_points))
prf_bnd(:) = 0
allocate(psi2d(n_points, n_plane))
psi2d = 0.d0
allocate(s_arc(n_points))
s_arc = 0.d0

do i_tor=1, n_tor
  write(*,*) ' toroidal mode numbers : ',i_tor,mode(i_tor)
enddo

do m=1,n_plane
  if (periodic) then
    phi = 2.d0 * PI * float(m-1)/float(n_plane)
  else
    phi = 2.d0 * PI * float(m-1)/float(n_plane-1) / float(n_period)
  endif
enddo

do m=1, n_plane

   inode_avg = 0

  if (periodic) then
    angle = 2.d0 * PI * float(m-1)/float(n_plane)
  else
    angle = 2.d0 * PI * float(m-1)/float(n_plane-1) / float(n_period)
  endif

  do i=1,element_list%n_elements

    do iv = 1, n_vertex_max

      iv1 = iv
      iv2 = mod(iv,4) + 1

      inode1 = element_list%element(i)%vertex(iv1)
      inode2 = element_list%element(i)%vertex(iv2)

      if (m .eq.1) then

        if ((node_list%node(inode1)%boundary .eq. 4) .and. (node_list%node(inode2)%boundary .eq. 5)) then
          write(*,*) ' PROBLEM mixed boundary',node_list%node(inode1)%boundary,node_list%node(inode2)%boundary
        endif
        if ((node_list%node(inode1)%boundary .eq. 5) .and. (node_list%node(inode2)%boundary .eq. 4)) then
          write(*,*) ' PROBLEM2 mixed boundary',node_list%node(inode1)%boundary,node_list%node(inode2)%boundary
        endif
        if ((node_list%node(inode1)%boundary .eq. 1) .and. (node_list%node(inode2)%boundary .eq. 5)) then
          write(*,*) ' PROBLEM4 mixed boundary',node_list%node(inode1)%boundary,node_list%node(inode2)%boundary
        endif
        if ((node_list%node(inode1)%boundary .eq. 5) .and. (node_list%node(inode2)%boundary .eq. 1)) then
          write(*,*) ' PROBLEM6 mixed boundary',node_list%node(inode1)%boundary,node_list%node(inode2)%boundary
        endif
        if ((node_list%node(inode1)%boundary .eq. 2) .and. (node_list%node(inode2)%boundary .eq. 4)) then
          write(*,*) ' PROBLEM7 mixed boundary',node_list%node(inode1)%boundary,node_list%node(inode2)%boundary
        endif
        if ((node_list%node(inode1)%boundary .eq. 2) .and. (node_list%node(inode2)%boundary .eq. 5)) then
          write(*,*) ' PROBLEM8 mixed boundary',node_list%node(inode1)%boundary,node_list%node(inode2)%boundary
        endif
        if ((node_list%node(inode1)%boundary .eq. 4) .and. (node_list%node(inode2)%boundary .eq. 2)) then
          write(*,*) ' PROBLEM9 mixed boundary',node_list%node(inode1)%boundary,node_list%node(inode2)%boundary
        endif
        if ((node_list%node(inode1)%boundary .eq. 5) .and. (node_list%node(inode2)%boundary .eq. 2)) then
          write(*,*) ' PROBLEM10 mixed boundary',node_list%node(inode1)%boundary,node_list%node(inode2)%boundary
        endif

      endif

      if   (( ((node_list%node(inode1)%boundary .eq. 1) .or. (node_list%node(inode1)%boundary .eq. 3))   &
       .and.  ((node_list%node(inode2)%boundary .eq. 1) .or. (node_list%node(inode2)%boundary .eq. 3)) ) &
       .or. ( ((node_list%node(inode1)%boundary .eq. 4) .or. (node_list%node(inode1)%boundary .eq. 9))   &
       .and.  ((node_list%node(inode2)%boundary .eq. 4) .or. (node_list%node(inode2)%boundary .eq. 9)) ) &
       .or. ( ((node_list%node(inode1)%boundary .eq. 4) .or. (node_list%node(inode1)%boundary .eq. 1))   &
       .and.  ((node_list%node(inode2)%boundary .eq. 4) .or. (node_list%node(inode2)%boundary .eq. 1)) ) &
       .or. ( ((node_list%node(inode1)%boundary .eq. 19) .or. (node_list%node(inode1)%boundary .eq. 11))   &
       .and.  ((node_list%node(inode2)%boundary .eq. 19) .or. (node_list%node(inode2)%boundary .eq. 11)) ) &
       .or. ( ((node_list%node(inode1)%boundary .eq. 5) .or. (node_list%node(inode1)%boundary .eq. 9))   &
       .and.  ((node_list%node(inode2)%boundary .eq. 5) .or. (node_list%node(inode2)%boundary .eq. 9)))) then 



        do j=1,nsub

         if   (  ((node_list%node(inode1)%boundary .eq. 1)  &
              .or.(node_list%node(inode1)%boundary .eq. 4)  &
              .or.(node_list%node(inode1)%boundary .eq. 9)  &
              .or.(node_list%node(inode1)%boundary .eq. 11)  &
              .or.(node_list%node(inode1)%boundary .eq. 3)) &
           .and. ((node_list%node(inode2)%boundary .eq. 1)  &
              .or.(node_list%node(inode2)%boundary .eq. 4)  &
              .or.(node_list%node(inode2)%boundary .eq. 9)  &
              .or.(node_list%node(inode2)%boundary .eq. 11)  &
              .or.(node_list%node(inode2)%boundary .eq. 3)) ) then
 
           sg = float(j-1)/float(nsub-1)

            if (iv1 .eq. 1) then
              tg     =  0.d0
              normal = -1.d0
            elseif (iv .eq. 3) then
              tg     = 1.d0
              normal = 1.d0
            else
              if ((m.eq.1) .and. (j.eq.1)) then
                write(*,*) ' problem 4/9 : ',i,iv1,iv2,inode1,inode2
                write(*,*) ' node1 R/Z   : ',node_list%node(inode1)%x(1,1,:),node_list%node(inode1)%boundary
                write(*,*) ' node2 R/Z   : ',node_list%node(inode2)%x(1,1,:),node_list%node(inode2)%boundary
              endif
            endif

         elseif (((node_list%node(inode1)%boundary .eq. 2) &
              .or.(node_list%node(inode1)%boundary .eq. 3)  &
              .or.(node_list%node(inode1)%boundary .eq. 5)  &
              .or.(node_list%node(inode1)%boundary .eq. 11)  &
              .or.(node_list%node(inode1)%boundary .eq. 9)) &
           .and. ((node_list%node(inode2)%boundary .eq. 2)  &
              .or.(node_list%node(inode2)%boundary .eq. 3)  &
              .or.(node_list%node(inode2)%boundary .eq. 5)  &
              .or.(node_list%node(inode2)%boundary .eq. 11)  &
              .or.(node_list%node(inode2)%boundary .eq. 9)) ) then

            tg = float(j-1)/float(nsub-1)

            if (iv1 .eq. 4) then
              sg     =  0.d0
              normal = -1.d0
            elseif (iv1 .eq.2) then
              sg     =  1.d0
              normal = +1.d0
            else
               if ((m .eq.1) .and. (j.eq.1)) write(*,*) ' problem 5/9 : ',i,iv1,iv2,inode1,inode2
            endif

          else
            if ((m .eq.1) .and. (j.eq.1)) then
              write(*,*) 'UNTREATED BOUNDARY : ',iv1,iv2,inode1,inode2,node_list%node(inode1)%boundary,node_list%node(inode2)%boundary
            endif
          endif

          call interp_RZ(node_list,element_list,i,sg,tg,R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt)
          BigR = R
          xjac = R_s * Z_t - R_t * Z_s

          inode = inode+1

          inode_avg = inode_avg + 1

          xyz(1:3,inode) = (/ R * cos(angle), Z, R*sin(angle) /)

          Ravg(inode_avg) = R
          Zavg(inode_avg) = Z


          psi = 0.d0; psi_s = 0.d0; psi_t = 0.d0;
          rho = 0.d0; rho_s = 0.d0; rho_t = 0.d0;
          rhon= 0.d0; rhon_s= 0.d0; rhon_t= 0.d0;
          T   = 0.d0; T_s = 0.d0;   T_t = 0.d0; T_p = 0.d0
          Vpar = 0.d0;Ti=0.0;Te=0.0;Ti_s=0.0;Ti_t=0.0;Ti_p=0.0;Te_s=0.0;Te_t=0.0;Te_p=0.0

          do i_tor = 1,n_tor
            
            if ( ( i_tor == 1 ) .and. ( without_n0_mode ) ) cycle ! Do not include the n=0 mode

            call interp(node_list,element_list,i,var_psi,i_tor,sg,tg,PS,PS_s,PS_t,PS_st,PS_ss,PS_tt)
            psi   = psi   + PS   * HZ(i_tor,m)
            psi_s = psi_s + PS_s * HZ(i_tor,m)
            psi_t = psi_t + PS_t * HZ(i_tor,m)

            call interp(node_list,element_list,i,var_rho,i_tor,sg,tg,RH,RH_s,RH_t,RH_st,RH_ss,RH_tt)
            rho   = rho   + RH   * HZ(i_tor,m)
            rho_s = rho_s + RH_s * HZ(i_tor,m)
            rho_t = rho_t + RH_t * HZ(i_tor,m)
            if (with_TiTe) then
              call interp(node_list,element_list,i,var_Ti,i_tor,sg,tg,TTi,TTi_s,TTi_t,TTi_st,TTi_ss,TTi_tt)
              Ti   = Ti   + TTi   * HZ(i_tor,m)
              Ti_s = Ti_s + TTi_s  * HZ(i_tor,m)
              Ti_t = Ti_t + TTi_t * HZ(i_tor,m)
              Ti_p = Ti_p + TTi   * HZ_p(i_tor,m)
              call interp(node_list,element_list,i,var_Te,i_tor,sg,tg,TTe,TTe_s,TTe_t,TTe_st,TTe_ss,TTe_tt)
              Te   = Te   + TTe   * HZ(i_tor,m)
              Te_s = Te_s + TTe_s  * HZ(i_tor,m)
              Te_t = Te_t + TTe_t * HZ(i_tor,m)
              Te_p = Te_p + TTe   * HZ_p(i_tor,m)
            end if

            call interp(node_list,element_list,i,var_T,i_tor,sg,tg,TT,TT_s,TT_t,TT_st,TT_ss,TT_tt)
            T   = T   + TT   * HZ(i_tor,m)                    
            T_s = T_s + TT_s * HZ(i_tor,m)
            T_t = T_t + TT_t * HZ(i_tor,m)
            T_p = T_p + TT   * HZ_p(i_tor,m)

            call interp(node_list,element_list,i,var_vpar,i_tor,sg,tg,VP,VP_s,VP_t,VP_st,VP_ss,VP_tt)
            Vpar = Vpar + VP * HZ(i_tor,m)              !V parallel

            call interp(node_list,element_list,i,var_rhon,i_tor,sg,tg,RH,RH_s,RH_t,RH_st,RH_ss,RH_tt)
            rhon   = rhon   + RH   * HZ(i_tor,m)
            rhon_s = rhon_s + RH_s * HZ(i_tor,m)
            rhon_t = rhon_t + RH_t * HZ(i_tor,m)

          enddo

          if (.not. with_TiTe) then
            Ti   = T / 2.d0
            Te   = T / 2.d0
            Ti_s = T_s / 2.d0
            Ti_t = T_t / 2.d0
            Ti_p = T_p / 2.d0
            Te_s = T_s / 2.d0
            Te_t = T_t / 2.d0
            Te_p = T_p / 2.d0
          end if

          psi_x = (   Z_t * psi_s - Z_s * psi_t ) / xjac
          psi_y = ( - R_t * psi_s + R_s * psi_t ) / xjac
          T_x   = (   Z_t * T_s   - Z_s * T_t )   / xjac
          T_y   = ( - R_t * T_s   + R_s * T_t )   / xjac

          BB2 = (F0**2 + (psi_x*psi_x+psi_y*psi_y)) / BigR**2
          
          psi_norm = get_psi_n(psi, Z)

          D_prof   = get_dperp (psi_norm)
          ZK_prof  = get_zkperp(psi_norm)

          if (ZKpar_T_dependent) then
            ZKpar_T  = ZK_i_par * abs(max(Ti,T_min)/T_0)**2.5
          else
            ZKpar_T  = ZK_i_par
          endif
          if (T .lt. ZK_prof_neg_thresh) then
            ZK_prof = ZK_prof_neg
          endif
          if (T .lt. ZK_par_neg_thresh) then
            ZKpar_T = ZK_par_neg
          endif

          scalars(inode,1) = psi
          scalars(inode,2) = rho
          scalars(inode,3) = Ti
          scalars(inode,4) = Vpar  * sqrt(BB2)
          scalars(inode,14)= Te
          psi2d(inode_avg, m) = psi   ! full precision psi for the 2D surface output

         if   (  ((node_list%node(inode1)%boundary .eq. 1)  &
              .or.(node_list%node(inode1)%boundary .eq. 4)  &
              .or.(node_list%node(inode1)%boundary .eq. 9)  &
              .or.(node_list%node(inode1)%boundary .eq. 11)  &
              .or.(node_list%node(inode1)%boundary .eq. 3)) &
           .and. ((node_list%node(inode2)%boundary .eq. 1)  &
              .or.(node_list%node(inode2)%boundary .eq. 4)  &
              .or.(node_list%node(inode2)%boundary .eq. 9)  &
              .or.(node_list%node(inode2)%boundary .eq. 11)  &
              .or.(node_list%node(inode2)%boundary .eq. 3)) ) then

            scalars(inode,5) = - (rho * Vpar * psi_s * normal)                  / R / sqrt(R_s**2 + Z_s**2)* sqrt(BB2)


            scalars(inode,6) = - (gamma_sheath -1.d0+gamma)*(rho * Ti * Vpar * psi_s * normal) / R / sqrt(R_s**2 + Z_s**2)

            if (abs(xjac) .gt. 1.d-7) then

            scalars(inode,7) = - ZKpar_T * (((T_x*psi_y - T_y*psi_x) + F0/BigR * T_p) / BigR / BB2 &
                             * (-psi_s*normal)) /BigR /sqrt(R_s**2+Z_s**2) 

            scalars(inode,8) = - ZK_prof * ( T_t   * (R_s**2 + Z_s**2) - T_s   * (R_s*R_t+Z_s*Z_t)) / Xjac &
                             / sqrt(R_s**2+Z_s**2) * normal

            scalars(inode,9) = - D_prof  * ( rho_t * (R_s**2 + Z_s**2) - rho_s * (R_s*R_t+Z_s*Z_t)) / Xjac &
                             / sqrt(R_s**2+Z_s**2) * normal
            endif

            scalars(inode,10) = - psi_s / BigR /sqrt(R_s**2+Z_s**2) * normal
            scalars(inode,11) = gamma_sheath * rho * T * Vpar * sqrt(BB2) 
            scalars(inode,12) = normal
            scalars(inode,13) = - rho * (Vpar * sqrt(BB2))**2 * Vpar * psi_s * normal / R / sqrt(R_s**2 + Z_s**2)

            vectors(inode,:,1) = (/ + psi_y /BigR * cos(angle), - psi_x /BigR, + psi_y /BigR * sin(angle) /) 
            vectors(inode,:,2) = (/ + vpar * psi_y /BigR* cos(angle), - vpar * psi_x /BigR, + vpar * psi_y /BigR * sin(angle) /) 
            vectors(inode,:,3) = (/ - Z_s * cos(angle), + R_s, -Z_s * sin(angle)  /) / sqrt(R_s**2+Z_s**2) * normal

            avg6(inode_avg) = avg6(inode_avg) + scalars(inode,6)
            avg7(inode_avg) = avg7(inode_avg) + scalars(inode,7)
            avg8(inode_avg) = avg8(inode_avg) + scalars(inode,8)
            avg11(inode_avg) = avg11(inode_avg) + scalars(inode,11)

         elseif (((node_list%node(inode1)%boundary .eq. 2) &
              .or.(node_list%node(inode1)%boundary .eq. 3)  &
              .or.(node_list%node(inode1)%boundary .eq. 5)  &
              .or.(node_list%node(inode1)%boundary .eq. 11)  &
              .or.(node_list%node(inode1)%boundary .eq. 9)) &
           .and. ((node_list%node(inode2)%boundary .eq. 2)  &
              .or.(node_list%node(inode2)%boundary .eq. 3)  &
              .or.(node_list%node(inode2)%boundary .eq. 5)  &
              .or.(node_list%node(inode2)%boundary .eq. 11)  &
              .or.(node_list%node(inode2)%boundary .eq. 9)) ) then

            scalars(inode,5) = (rho * Vpar * psi_t * normal)                    / R / sqrt(R_t**2 + Z_t**2)

            scalars(inode,6) = gamma_sheath * (rho * T * Vpar * psi_t * normal) / R / sqrt(R_t**2 + Z_t**2)
 
            if (abs(xjac) .gt. 1.d-7) then

            scalars(inode,7) = - ZKpar_T * (((T_x*psi_y - T_y*psi_x) + F0/BigR * T_p) / BigR / BB2 &
                             * (psi_t*normal)) /BigR /sqrt(R_t**2+Z_t**2) 

            scalars(inode,8) = - ZK_prof * ( T_s   * (R_t**2 + Z_t**2) - T_t   * (R_s*R_t+Z_s*Z_t)) / Xjac &
                             / sqrt(R_t**2+Z_t**2) * normal

            scalars(inode,9) = - D_prof  * ( rho_s * (R_t**2 + Z_t**2) - rho_t * (R_s*R_t+Z_s*Z_t)) / Xjac &
                             / sqrt(R_t**2+Z_t**2) * normal

            endif

            scalars(inode,10) = psi_t / BigR /sqrt(R_t**2+Z_t**2) * normal
            scalars(inode,11) = gamma_sheath * rho * Ti * Vpar * sqrt(BB2) 
            scalars(inode,12) = normal
            scalars(inode,13) = -rho * (Vpar * sqrt(BB2))**2 * Vpar * psi_t * normal / R / sqrt(R_t**2 + Z_t**2)

            vectors(inode,:,1) = (/ + psi_y /BigR*cos(angle)+F0*sin(angle)/BigR,-psi_x/BigR,+psi_y/BigR *sin(angle)+F0*cos(angle)/BigR /) 
            vectors(inode,:,2) = (/ + vpar * psi_y /BigR* cos(angle), - vpar * psi_x /BigR, + vpar * psi_y /BigR * sin(angle)  /) 
            vectors(inode,:,3) = (/ + Z_t * cos(angle), - R_t, Z_t * sin(angle)  /) / sqrt(R_t**2+Z_t**2) * normal

            avg6(inode_avg) = avg6(inode_avg) + scalars(inode,6)
            avg7(inode_avg) = avg7(inode_avg) + scalars(inode,7)
            avg8(inode_avg) = avg8(inode_avg) + scalars(inode,8)
            avg11(inode_avg) = avg11(inode_avg) + scalars(inode,11)

          else

            !write(*,*) 'warning untreated boundary : ',inode1,inode2,node_list%node(inode1)%boundary,node_list%node(inode2)%boundary

          endif

          if (m.eq.1) then
            Rprf(inode_avg)  = R
            Zprf(inode_avg)  = Z
	        prf1(inode_avg)  = scalars(inode,1)
            prf7(inode_avg)  = scalars(inode,7) / MU_zero / t_norm * 1.5
            prf6(inode_avg)  = scalars(inode,6) / MU_zero / t_norm * 1.5
            prf2(inode_avg)  = scalars(inode,2) * central_density
            prf3(inode_avg)  = scalars(inode,3) / MU_zero / (central_density * 1d20) / 1.602d-19
            prf4(inode_avg)  = scalars(inode,4) / t_norm
            prf5(inode_avg)  = scalars(inode,5) * central_density / t_norm
            prf10(inode_avg) = acos(scalars(inode,10))*180/PI
            prf11(inode_avg) = scalars(inode,11) / MU_zero / t_norm
            prf14(inode_avg) = scalars(inode,14) / MU_zero / (central_density * 1d20) / 1.602d-19
            B_field(inode_avg,1:3) = vectors(inode,:,1)
            B_field(inode_avg,4) = sqrt(BB2)
            ! Determine target plate ID: require BOTH nodes to have the same target type
            ! 1 = target (type 1 or 3),  0 = wall/other
            bnd1 = node_list%node(inode1)%boundary
            bnd2 = node_list%node(inode2)%boundary
            if ( ((bnd1 == 1 .or. bnd1 == 3) .and. (bnd2 == 1 .or. bnd2 == 3)) ) then
              prf_bnd(inode_avg) = 1
            else
              prf_bnd(inode_avg) = 0
            end if
          endif

        enddo

        if (m .lt. n_plane) then

          do j=1,nsub-1

            ielm        = ielm + 1
            ien(1,ielm) = inode - nsub + (j-1)      ! 0 based indices for VTK
            ien(2,ielm) = inode - nsub + (j  ) 
            ien(3,ielm) = ien(2,ielm) + n_points
            ien(4,ielm) = ien(1,ielm) + n_points

          enddo

        endif

        if ( (periodic) .and. (m .eq. n_plane)) then

          do j=1,nsub-1

            ielm        = ielm+1
            ien(1,ielm) = inode - nsub + (j-1)       ! 0 based indices for VTK
            ien(2,ielm) = inode - nsub + (j  ) 

            ien(3,ielm) = ien(2,ielm) - n_points * (n_plane-1)
            ien(4,ielm) = ien(1,ielm) - n_points * (n_plane-1)

          enddo

        endif
	
      endif

    enddo
  enddo
enddo

close(22)


!scalar_names = (/ 'flux    ','density ','T       ','Vpar    ','nV.n    ','nTV.n   ','KparT.n ', &
!                  'Kperp.T ','Dperp.n ','B.n     '/)

scalars(:,2) = scalars(:,2) * central_density
scalars(:,3) = scalars(:,3) / MU_zero / (central_density * 1d20) / 1.602d-19
scalars(:,4) = scalars(:,4) / t_norm
scalars(:,5) = scalars(:,5) * central_density / t_norm
scalars(:,6) = scalars(:,6) / MU_zero / t_norm * 1.5
scalars(:,7) = scalars(:,7) / MU_zero / t_norm * 1.5
scalars(:,8) = scalars(:,8) / MU_zero / t_norm * 1.5
scalars(:,9) = scalars(:,9) * central_density / t_norm
scalars(:,11) = scalars(:,11) / MU_zero / t_norm
scalars(:,13) = scalars(:,13) / MU_zero / t_norm * 0.5
scalars(:,14) = scalars(:,14) / MU_zero /(central_density * 1d20) / 1.602d-19
  ! Sputter calculation removed (use Python script instead)
  ! call Sputter(n_points,scalars(:,3),prf10(:),scalars(:,5),yield(:))
  scalars(:,15) = 0.d0
  prf15(:) = 0.d0
  yield(:) = 0.d0

avg6(:) = avg6(:) / MU_zero / t_norm / float(n_plane) * 1.5
avg7(:) = avg7(:) / MU_zero / t_norm / float(n_plane) * 1.5
avg8(:) = avg8(:) / MU_zero / t_norm / float(n_plane) * 1.5
avg11(:)= avg11(:) / MU_zero / t_norm / float(n_plane)


open(22,file='target_profile')
write(22,'(A132)') "R               Z               angle           KparT_normal    gam_nVT_normal  density         Ti              Vpar            nv_normal       gam_nvT      Te"

open(23,file='average_target_profile')
!write(23,'(A132)') '      time         step            Length         R               Z              nTV.n           KparT.n        Kperp.T        nvT_gam'
open(24,file='B_field.dat')
write(24,'(A132)') 'B_R B_Z B_phi     B0  '
open(25,file='node_boundary')

do i=1,node_list%n_nodes
  write(25,'(2e16.8,I3)') node_list%node(i)%x(1,1,1),node_list%node(i)%x(1,1,2),node_list%node(i)%boundary
enddo


allocate(tobedone(n_points))

allocate(order_idx(n_points))  ! ordered indices from nearest-neighbor walk

tobedone = .true.
distance = 0.d0
print *,n_points

do i=1, n_points

  index = minval(minloc((Rprf-R_start)**2+(Zprf-Z_start)**2,tobedone))

  distance = distance + sqrt((Rprf(index)-R_start)**2+(Zprf(index)-Z_start)**2)

  R_start = Rprf(index)
  Z_start = Zprf(index)
  tobedone(index) = .false.

  order_idx(i) = index  ! save ordering

  write(22,'(12e16.8)') Rprf(index),Zprf(index),prf10(index),prf7(index),prf6(index), &
                        prf2(index), prf3(index),  prf4(index), prf5(index), prf11(index) ,prf14(index),prf15(index)

  write(23,'(12e16.8)') xtime(index_start)*t_norm,xtime(index_start), distance, &
                        Ravg(index),Zavg(index),avg6(index),avg7(index),avg8(index),avg11(index),normal

 write(24,'(4e16.8)') B_field(index,1),B_field(index,2),B_field(index,3),B_field(index,4)
enddo

close(23)
close(22)
close(24)

! --- Arc length along the boundary in the ordered walk (for the 2D surface output)
s_cum2d = 0.d0
s_arc(1) = 0.d0
do i2d = 2, n_points
  s_cum2d  = s_cum2d + sqrt((Rprf(order_idx(i2d))-Rprf(order_idx(i2d-1)))**2 &
                          + (Zprf(order_idx(i2d))-Zprf(order_idx(i2d-1)))**2)
  s_arc(i2d) = s_cum2d
enddo

! ============================================================================
! 2D surface output (midplane2d-style text file): target_surface2d.dat
! Layout: one '# Column m of n_plane' block per toroidal plane m (like the
! jorek2_postproc midplane2d files), n_points rows per column ordered along
! the boundary by order_idx, each row holds all values (es23.15).
! Columns: s=arc length along boundary [m], R, Z, psi (raw poloidal flux),
! Psi_N = normalized poloidal flux, rho = sqrt(Psi_N), density, Ti, Te,
! Vpar, angle (B.n incidence [deg]), nv, nvT_gam, KparT (target fluxes,
! only meaningful where bnd==1), bnd (1 = target plate, 0 = wall/other).
! ============================================================================
open(28, file='target_surface2d.dat')
write(28,'(A,15A23)') '# ', 's', 'R', 'Z', 'psi', 'Psi_N', 'rho', 'density', &
                      'Ti', 'Te', 'Vpar', 'angle', 'nv', 'nvT_gam', 'KparT', 'bnd'

do m2d = 1, n_plane
  if (periodic) then
    phi_deg2d = 360.d0 * float(m2d-1) / float(n_plane)
  else
    phi_deg2d = 360.d0 * float(m2d-1) / float(n_plane-1) / float(n_period)
  endif
  write(28,'(A,I4.4,A,I4.4,A,F8.3)') '  # Column ', m2d, ' of ', n_plane, '  phi(deg)=', phi_deg2d

  do i2d = 1, n_points
    idx2d   = order_idx(i2d)
    inode2d = (m2d-1)*n_points + idx2d
    psi_n2d = get_psi_n(psi2d(idx2d, m2d), Zprf(idx2d))
    rho2d   = sqrt(max(psi_n2d, 0.d0))
    write(28,'(9999es23.15)') s_arc(i2d), Rprf(idx2d), Zprf(idx2d), &
      psi2d(idx2d, m2d), psi_n2d, rho2d, &
      scalars(inode2d,2), scalars(inode2d,3), scalars(inode2d,14), scalars(inode2d,4), &
      acos(max(-1.d0,min(1.d0,scalars(inode2d,10))))*180.d0/PI, &
      scalars(inode2d,5), scalars(inode2d,11), scalars(inode2d,7), dble(prf_bnd(idx2d))
  end do
end do
close(28)
write(*,'(A)') ' target_surface2d.dat: 2D surface output written.'

! ============================================================================
! Output target plate profiles with signed distance from strike point
! Approach (like MATLAB script): filter target-only points, then segment
! plates by geometric distance jumps between consecutive points.
! strike_dist < 0 : private flux region (inside separatrix)
! strike_dist > 0 : scrape-off layer (outside separatrix)
! ============================================================================

! --- Step 0: Find strike points on boundary elements using corrected ES%psi_bnd
ES%num_strike = 0
do ibndelm = 1, bnd_elm_list%n_bnd_elements
  iv1 = bnd_elm_list%bnd_element(ibndelm)%vertex(1)
  iv2 = bnd_elm_list%bnd_element(ibndelm)%vertex(2)
  psi1 = node_list%node(iv1)%values(1,1,1)
  psi2 = node_list%node(iv2)%values(1,1,1)
  if ((min(psi1,psi2) < ES%psi_bnd) .and. (max(psi1,psi2) > ES%psi_bnd)) then
    ES%num_strike = ES%num_strike + 1
    frac = (ES%psi_bnd - psi1) / (psi2 - psi1)
    frac = max(0.d0, min(1.d0, frac))
    ES%R_strike(ES%num_strike) = (1.d0-frac)*node_list%node(iv1)%x(1,1,1) + frac*node_list%node(iv2)%x(1,1,1)
    ES%Z_strike(ES%num_strike) = (1.d0-frac)*node_list%node(iv1)%x(1,1,2) + frac*node_list%node(iv2)%x(1,1,2)
  end if
end do
write(*,'(A,I2,A)') ' Found ', ES%num_strike, ' strike points on boundary elements.'
open(27, file='strike_points.dat')
do i = 1, ES%num_strike
  write(*,'(A,I1,A,F10.4,A,F10.4)') '   Strike ', i, ': R=', ES%R_strike(i), '  Z=', ES%Z_strike(i)
  write(27,'(F12.6,2x,F12.6)') ES%R_strike(i), ES%Z_strike(i)
end do
close(27)

! --- Step 1: Extract target-only points from the ordered list
allocate(target_ord_idx(n_points))
n_target_pts = 0
do i = 1, n_points
  idx = order_idx(i)
  if (prf_bnd(idx) == 1) then
    n_target_pts = n_target_pts + 1
    target_ord_idx(n_target_pts) = idx
  end if
end do

if (n_target_pts == 0) then
  write(*,*) 'WARNING: No target plate points found!'
else

open(26, file='target_strike_profile.dat')
write(26,'(A)') '# plate_id  strike_dist  R  Z  angle  density  Ti  Te  Vpar  nv  nvT_gam  KparT  psi  sputter'

! --- Step 2: Segment plates by distance jumps
strike_used(:) = .false.
plate_id = 0
p_start = 1

do i = 1, n_target_pts
  is_last = (i == n_target_pts)

  if (.not. is_last) then
    idx1 = target_ord_idx(i)
    idx2 = target_ord_idx(i+1)
    dgap = sqrt((Rprf(idx2)-Rprf(idx1))**2 + (Zprf(idx2)-Zprf(idx1))**2)
    is_gap = (dgap > jump_threshold)
  else
    is_gap = .true.
  end if

  if (is_gap) then
    p_end = i
    plate_id = plate_id + 1

    ! --- Match this plate to nearest strike point (not yet assigned to another plate)
    strike_R = 0.d0
    strike_Z = 0.d0
    strike_found = .false.
    dist_min_strike = 1.d99
    do k = 1, ES%num_strike
      if (strike_used(k)) cycle  ! already assigned to a previous plate
      do j = p_start, p_end
        jdx = target_ord_idx(j)
        dist_ps = sqrt((Rprf(jdx)-ES%R_strike(k))**2 + (Zprf(jdx)-ES%Z_strike(k))**2)
        if (dist_ps < dist_min_strike) then
          dist_min_strike = dist_ps
          strike_R = ES%R_strike(k)
          strike_Z = ES%Z_strike(k)
          strike_found = .true.
          assigned_k = k
        end if
      end do
    end do

    if (strike_found) strike_used(assigned_k) = .true.

    if (.not. strike_found) then
      strike_R = Rprf(target_ord_idx(p_start))
      strike_Z = Zprf(target_ord_idx(p_start))
    end if

    write(26,'(A,I2,A,F10.4,A,F10.4)') '# Plate ', plate_id, &
      '  strike_R=', strike_R, '  strike_Z=', strike_Z
    write(*,'(A,I2,A,F10.4,A,F10.4)') ' Plate ', plate_id, &
      ' matched to strike at R=', strike_R, ' Z=', strike_Z

    ! --- Copy plate points to temp array, sort by Z for monotonic arc-length
    n_plate_pts = p_end - p_start + 1
    allocate(plate_idx(n_plate_pts))
    do j = 1, n_plate_pts
      plate_idx(j) = target_ord_idx(p_start + j - 1)
    end do

    ! Sort by Z (ascending for upper targets, descending for lower)
    do j = 1, n_plate_pts-1
      do k = j+1, n_plate_pts
        jdx_j = plate_idx(j)
        jdx_k = plate_idx(k)
        if (strike_Z > 0.d0) then
          if (Zprf(jdx_j) > Zprf(jdx_k)) then
            plate_idx(j) = jdx_k; plate_idx(k) = jdx_j
          end if
        else
          if (Zprf(jdx_j) < Zprf(jdx_k)) then
            plate_idx(j) = jdx_k; plate_idx(k) = jdx_j
          end if
        end if
      end do
    end do

    ! --- Compute strike point arc-length in sorted order
    dist_min_strike_s = 1.d99
    strike_dist_sorted = 0.d0
    cumdist = 0.d0
    do j = 1, n_plate_pts
      jdx = plate_idx(j)
      if (j > 1) then
        pdx = plate_idx(j-1)
        cumdist = cumdist + sqrt((Rprf(jdx)-Rprf(pdx))**2 + (Zprf(jdx)-Zprf(pdx))**2)
      end if
      dist_ps = sqrt((Rprf(jdx)-strike_R)**2 + (Zprf(jdx)-strike_Z)**2)
      if (dist_ps < dist_min_strike_s) then
        dist_min_strike_s = dist_ps
        strike_dist_sorted = cumdist
      end if
    end do

    ! --- Write sorted plate with signed distances
    cumdist = 0.d0
    do j = 1, n_plate_pts
      jdx = plate_idx(j)

      if (j > 1) then
        pdx = plate_idx(j-1)
        cumdist = cumdist + sqrt((Rprf(jdx)-Rprf(pdx))**2 + (Zprf(jdx)-Zprf(pdx))**2)
      end if

      sdist = cumdist - strike_dist_sorted
      ! Sign: private region -> negative
      if (strike_Z > 0.d0) then
        if (abs(Zprf(jdx) - strike_Z) > 0.01d0) then
          if (Zprf(jdx) > strike_Z) sdist = -abs(sdist)
          if (Zprf(jdx) < strike_Z) sdist =  abs(sdist)
        else
          if ((Rprf(jdx) - strike_R) * (ES%R_xpoint(1) - strike_R) > 0.d0) sdist = -abs(sdist)
          if ((Rprf(jdx) - strike_R) * (ES%R_xpoint(1) - strike_R) < 0.d0) sdist =  abs(sdist)
        end if
      else
        if (abs(Zprf(jdx) - strike_Z) > 0.01d0) then
          if (Zprf(jdx) < strike_Z) sdist = -abs(sdist)
          if (Zprf(jdx) > strike_Z) sdist =  abs(sdist)
        else
          if ((Rprf(jdx) - strike_R) * (ES%R_xpoint(1) - strike_R) > 0.d0) sdist = -abs(sdist)
          if ((Rprf(jdx) - strike_R) * (ES%R_xpoint(1) - strike_R) < 0.d0) sdist =  abs(sdist)
        end if
      end if

      write(26,'(I3,2x,F12.6,12E14.6)') plate_id, sdist, Rprf(jdx), Zprf(jdx), prf10(jdx), &
        prf2(jdx), prf3(jdx), prf14(jdx), prf4(jdx), prf5(jdx), prf11(jdx), prf7(jdx), &
        prf1(jdx), prf15(jdx)
    end do

    deallocate(plate_idx)

    write(26,*)  ! blank line between plates
    write(26,*)

    p_start = i + 1
  end if
end do

close(26)
write(*,'(A,I2,A,I5,A)') ' target_strike_profile.dat: ', plate_id, ' plates, ', n_target_pts, ' points written.'

end if

deallocate(order_idx)
deallocate(target_ord_idx)
deallocate(psi2d)
deallocate(s_arc)



end

! Sputter subroutine removed — use Python script (util/plot_target_strike.py)
! for physically correct Eckstein sputtering calculation.


