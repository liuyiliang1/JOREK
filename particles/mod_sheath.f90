!A module about sheath calculation
!We need to calculate the sheath potential and the sheath thickness
!The sheath potential is calculated by the Boltzmann relation
!The sheath thickness is calculated by the Child-Langmuir law
!We also need to calculate the electric field in the sheath and the density of the electrons in the sheath
!To calculate these parameters, We need to know the distance between particles and the target
!So We calculate the distance between particles and the target
!We need to know the position of the the target
module mod_sheath
    implicit none
    public :: calc_E_sheath,calc_dis
contains

!We assume that the sheath is the Brooks model
!We need to know the temperature of the electrons and the magnetic field
!We need to know the density of the electrons
!We need to know the distance between the particle and the target
!We need to calculate the sheath potential and the sheath thickness
!We also need to calculate the electric field in the sheath and the density of the electrons in the sheath
!You can call Ealc_E_sheath(feedback_nodelist, feedback_element_list,edge_elm_template,particle_tmp%i_elm,particle_tmp%x,particle_tmp%st,kTb/EL_CHG,kTb/EL_CHG,B,n_b,E_sheath) to get the electric field and the sheath potential
subroutine calc_E_sheath(node_list,element_list,edge_element_template,i_elm,x,st,Te,Ti,B,ne,E_sheath,phi_optional)
  use data_structure
  use mod_edge_domain
  use mod_edge_elements
  use mod_boundary, only: wall_normal_vector
  use mod_interp
  use constants
  use phys_module,only: R_limiter,Z_limiter,n_limiter,T_min_neg,central_density
  use equil_info, only : ES
    implicit none
    type(edge_elements),intent(in)                               :: edge_element_template
    type(type_node_list), intent(in)                             :: node_list
    type(type_element_list), intent(in)                          :: element_list
    integer, intent(in) :: i_elm !the number of elements which is the particle in
    real*8, dimension(3), intent(in) :: x !the position of the particle
    real*8, dimension(2), intent(in) :: st !the position of the target
    real*8, intent(in) :: Te !the temperature of the electrons,eV
    real*8, intent(in) :: Ti !the temperature of the ions, eV
    real*8, intent(in) :: B(3) !the magnetic field
    real*8, intent(in) :: ne !the density of the electrons,m^-3
    real*8, intent(inout) :: E_sheath(3) !the electric field in the sheath
    real*8 :: E_sheath_0 !the electric field in the sheath
    real*8,intent(inout),optional :: phi_optional !the sheath potential
    real*8 :: phi_0,phi_1,phi_2 !the sample sheath potential
    real*8 :: phi !the sheath potential
    real*8 :: R   !the radius of the particle
    real*8 :: lamuda !the debye length
    real*8 :: dis,dis1
    real*8 :: Te_safe, ne_safe, Ti_safe, T_min_neg_eV
    integer:: i_0,k_0,i_limiter,i_elm_edge
    integer :: i,j,k,i_patch !loop index
    real*8 :: vector_normal(3),r_vector(2) !< vector normal to the wall
    real*8:: l_limiter(2,n_limiter),vec_limiter(2,n_limiter)
    E_sheath=0.d0
    ! convert T_min_neg from JOREK to eV (same as calc_NeTeTi: T_K = T_jorek * T_norm)
    ! T_norm = 1/(K_BOLTZ * 2 * MU_ZERO * central_density * 1d20) for single-T
    T_min_neg_eV = T_min_neg / (2.d0 * MU_ZERO * central_density * 1.d20 * EL_CHG)
    if (ne .le. 0.d0) return
    Te_safe = Te
    if (Te_safe .le. 0.d0) then
      if (T_min_neg_eV .le. 0.d0) return
      Te_safe = T_min_neg_eV
    end if
    Ti_safe = max(Ti, T_min_neg_eV)
    ne_safe = ne
    vector_normal=0.d0
    if(present(phi_optional)) phi_optional=0
    l_limiter(1,1:n_limiter)=R_limiter(1:n_limiter)
    l_limiter(2,1:n_limiter)=Z_limiter(1:n_limiter)
    vec_limiter(1:2,1:n_limiter-1)=l_limiter(1:2,2:n_limiter)-l_limiter(1:2,1:n_limiter-1)
    vec_limiter(1:2,n_limiter)=l_limiter(1:2,1)-l_limiter(1:2,n_limiter)
    if (i_elm .le. 0) then
        write(*,*) 'Something wrong in the sheath calculation with i_elm'
         return!< .not. .lt.!< if this is not a lost particle go to next particle
    end if
    call calc_dis(node_list,element_list,i_elm,x,edge_element_template,dis,i_elm_edge)
    if(dis.gt.0.005) return !< if the distance is larger than 0.5 cm, we need not calculate the sheath, set E_sheath=0

vector_normal=-wall_normal_vector(node_list, element_list, i_elm_edge, st(1), st(2))
!now we have the distance between the particle and the target
!We have already know the temperature of the electrons and the magnetic field  and the density of the electrons
!now We can calculate the sheath potential and the sheath thickness
phi_0=-3*Te_safe  !the sample sheath potential
lamuda=sqrt(EPS_ZERO*Te_safe/(ne_safe*EL_CHG)) !the debye length
R=0.5d0*sqrt(2*Ti_safe*2*ATOMIC_MASS_UNIT*EL_CHG)/(norm2(B)*EL_CHG) !the radius of the particle
phi_1=0.25*phi_0
phi_2=0.75*phi_0
phi=phi_1*exp(-dis/(2*lamuda)) + phi_2*exp(-dis/R) !the sheath potential
E_sheath_0=-(phi_1/(2*lamuda)*exp(-dis/(2*lamuda))+phi_2/R*exp(-dis/R)) !the electric field in the sheath
E_sheath=E_sheath_0*vector_normal !the electric field in the sheath
end subroutine calc_E_sheath

subroutine calc_dis(node_list, element_list, i_elm, x, edge_element_template, dis0, i_elm_edge)
    use data_structure
    use mod_edge_domain
    use mod_edge_elements
    use mod_boundary, only: wall_normal_vector
    use mod_interp
    use constants
    implicit none
    type(edge_elements), intent(in) :: edge_element_template
    type(type_node_list), intent(in) :: node_list
    type(type_element_list), intent(in) :: element_list
    integer, intent(in) :: i_elm
    real*8, dimension(3), intent(in) :: x
    real*8, intent(inout) :: dis0
    integer, optional, intent(out) :: i_elm_edge
    real*8 :: dis
    integer :: i_patch, k(1), i_patch0, k0
    real*8 :: vector_normal(3), r_vector(2), r_vector0(2)
    real*8 :: min_dis
    real*8, allocatable :: dis_array(:)

    if (i_elm <= 0) return

    min_dis = 1.0e6
    i_patch0 = 0
    k0 = 0
    do i_patch = 1,size(edge_element_template%patch,1)
        allocate(dis_array(size(edge_element_template%patch(i_patch)%xyz, 2)))
        dis_array = sqrt((edge_element_template%patch(i_patch)%xyz(1, :) - x(1))**2 + &
                         (edge_element_template%patch(i_patch)%xyz(2, :) - x(2))**2)
        k = MINLOC(dis_array)
        dis=dis_array(k(1))
        if (dis < min_dis) then
            min_dis = dis
            i_patch0 = i_patch
            k0 = k(1)
            r_vector0 = x(1:2) - edge_element_template%patch(i_patch)%xyz(1:2, k0)
        endif
        deallocate(dis_array)
    enddo
    dis0 = min_dis

    if (min_dis > 0.01) return

    vector_normal = wall_normal_vector(node_list, element_list, &
                                       edge_element_template%patch(i_patch0)%i_elm_jorek_edge(k0), &
                                       edge_element_template%patch(i_patch0)%st(1, k0), &
                                       edge_element_template%patch(i_patch0)%st(2, k0))
    dis0 = abs(dot_product(r_vector0, vector_normal(1:2)))
    if (present(i_elm_edge)) i_elm_edge = edge_element_template%patch(i_patch0)%i_elm_jorek_edge(k0)
end subroutine calc_dis

end module mod_sheath
