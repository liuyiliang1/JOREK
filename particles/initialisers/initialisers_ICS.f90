!> module containing initialization subroutines for ICS (impurity coupling scheme)
!> particles, with spatial distribution proportional to background electron density
module initialisers_ICS
  use mod_particle_types
  use mod_particle_sim
  use mod_rng
  use initialisers_base
  use mod_model_settings, only: var_rho
  use mod_coronal
  use mod_sampling, only: transform_uniform_cylindrical
  use mod_interp,  only: interp_0
  use data_structure
  use phys_module, only: central_density, n_period
  use mpi
  implicit none

  contains

!> Initialize ICS particles with spatial positions distributed ∝ n_e via rejection sampling.
!> Entirely sequential (no OpenMP) to avoid allocatable-array issues in parallel regions.
!> Weight = concentration * total_main_ions / n_particles.
subroutine ics_density_initialization(sim, group_num, concentration)
  use mod_pcg32_rng

  type(particle_sim), intent(inout), target :: sim
  integer,            intent(in)    :: group_num
  real*8,             intent(in)    :: concentration

  class(particle_base), pointer     :: particles(:)
  type(type_node_list),    pointer  :: node_list
  type(type_element_list), pointer  :: element_list
  class(type_rng), allocatable      :: rng
  integer  :: i, j, n_particles_global, ifail, i_elm, n_found, n_left, seed, ierr
  real*8   :: Rbox(2), Zbox(2), Phibox(2)
  real*8   :: R, Z, phi, s, t, DUMMY_REAL
  real*8   :: ran(7), P_interp(1)
  real*8   :: t0, t1, total_ions, weight_per_particle
  integer, allocatable :: i_to_find(:)
  logical, allocatable :: not_found(:)

  node_list    => sim%fields%node_list
  element_list => sim%fields%element_list
  particles    => sim%groups(group_num)%particles
  !> use the global total number of particles (particles are distributed across MPI ranks)
  n_particles_global = int(sim%groups(group_num)%n_particles)

  ! Build coronal equilibrium table
  sim%groups(group_num)%cor = coronal(sim%groups(group_num)%ad)

  ! N_imp = c/(1-c) * N_main, where c = N_imp/(N_imp+N_main) is the concentration
  ! total_main_ions = central_density * 1d20 * ∫ rho dV  [physical ions]
  total_ions = compute_total_ions(node_list, element_list)
  weight_per_particle = concentration / (1.d0 - concentration) * total_ions / dble(n_particles_global)
  if (sim%my_id == 0) write(*,'(A,ES12.4,A,ES12.4)') &
      ' total_main_ions = ', total_ions, '  weight_per_particle = ', weight_per_particle

  ! Domain bounding box
  call domain_bounding_box(node_list, element_list, Rbox(1), Rbox(2), Zbox(1), Zbox(2))
  Phibox = [0.d0, 6.283185307179586d0]  ! 2*pi

  ! Setup RNG
  allocate(rng, source=pcg32_rng())
  seed = 123456789
  call rng%initialize(7, seed, 1, 1, ifail)

  ! Particles to find (local count on this MPI rank)
  allocate(i_to_find(size(particles)), not_found(size(particles)))
  i_to_find = [(i, i=1, size(particles))]
  not_found = .true.

  call cpu_time(t0)
  n_left = size(particles)
  do while (n_left > 0)
    do i = 1, n_left
      j = i_to_find(i)

      call rng%next(ran)
      call transform_uniform_cylindrical(ran(1:3), Rbox, Zbox, Phibox, R, Z, phi)

      call find_RZ(node_list, element_list, R, Z, DUMMY_REAL, DUMMY_REAL, i_elm, s, t, ifail)
      if (ifail == 0) then
        call interp_0(node_list, element_list, i_elm, [var_rho], 1, s, t, phi, P_interp)

        if (ran(4) < P_interp(1)) then
          ! Verify via find_RZ_nearby that the particle is robustly placed
          call verify_nearest_position(node_list, element_list, R, Z, s, t, i_elm, ifail)
          if (ifail /= 0) cycle  ! retry with new random position
          particles(j)%x     = [R, Z, phi]
          particles(j)%i_elm = i_elm
          particles(j)%st    = [s, t]
          select type (pa => particles(j))
          type is (particle_kinetic_leapfrog)
            pa%v = ran(5:7)
          end select
          not_found(i) = .false.
        end if
      end if
    end do

    i_to_find = pack(i_to_find, not_found)
    deallocate(not_found)
    n_left = size(i_to_find)
    allocate(not_found(n_left))
    not_found = .true.
  end do

  call cpu_time(t1)
  if (sim%my_id == 0) write(*,'(A,2f12.4)') ' ICS density init cpu/wall :', t1-t0

  ! Velocity: Maxwellian at local Ti=Te=T/2, with vpar and coronal charge state
  call set_velocity_from_T(particles, sim%groups(group_num)%mass, &
      node_list, element_list, sim%groups(group_num)%cor, v_par=.true.)

  ! Weight = concentration * total_main_ions / n_particles
  particles(:)%weight = weight_per_particle

  ! Write particle positions, charge, and kinetic energy to VTK for visualisation
  call write_particles_vtk(particles, sim%groups(group_num)%mass, sim%my_id, sim%groups(group_num)%id)

end subroutine ics_density_initialization

!> Write particle (R,Z) positions, charge and kinetic energy to ASCII VTK PolyData file
subroutine write_particles_vtk(particles, mass, my_id, group_id)
  use constants, only: EL_CHG, ATOMIC_MASS_UNIT
  class(particle_base), intent(in) :: particles(:)
  real*8,               intent(in) :: mass
  integer,              intent(in) :: my_id
  character(len=3),     intent(in) :: group_id
  integer :: i, n, iunit
  real*8  :: E_kin_eV
  character(len=128) :: fname

  if (my_id /= 0) return
  n = size(particles)
  write(fname, '(A,A,A)') 'particles_init_', trim(group_id), '.vtk'
  open(newunit=iunit, file=trim(fname), status='replace')
  write(iunit,'(A)') '# vtk DataFile Version 3.0'
  write(iunit,'(A)') 'ICS initial particle distribution'
  write(iunit,'(A)') 'ASCII'
  write(iunit,'(A)') 'DATASET POLYDATA'
  write(iunit,'(A,I12,A)') 'POINTS ', n, ' float'
  do i = 1, n
    write(iunit,'(3ES15.5)') particles(i)%x(1), particles(i)%x(2), 0.0
  end do
  write(iunit,'(A,I12,A,I12)') 'VERTICES ', n, ' ', 2*n
  do i = 1, n
    write(iunit,'(A,I12)') '1 ', i-1
  end do
  write(iunit,'(A,I12)') 'POINT_DATA ', n
  write(iunit,'(A)') 'SCALARS charge float 1'
  write(iunit,'(A)') 'LOOKUP_TABLE default'
  do i = 1, n
    select type (pa => particles(i))
    type is (particle_kinetic_leapfrog)
      write(iunit,'(I4)') int(pa%q)
    class default
      write(iunit,'(I4)') 0
    end select
  end do
  write(iunit,'(A)') 'SCALARS E_kin_eV float 1'
  write(iunit,'(A)') 'LOOKUP_TABLE default'
  do i = 1, n
    select type (pa => particles(i))
    type is (particle_kinetic_leapfrog)
      E_kin_eV = 0.5d0 * mass * ATOMIC_MASS_UNIT * dot_product(pa%v, pa%v) / EL_CHG
      write(iunit,'(ES15.5)') E_kin_eV
    class default
      write(iunit,'(ES15.5)') 0.d0
    end select
  end do
  write(iunit,'(A)') 'SCALARS weight float 1'
  write(iunit,'(A)') 'LOOKUP_TABLE default'
  do i = 1, n
    write(iunit,'(ES15.5)') particles(i)%weight
  end do
  close(iunit)
  if (my_id == 0) write(*,*) 'Wrote particle VTK: ', trim(fname)
end subroutine write_particles_vtk

!> Verify that a particle placed by find_RZ can be re-located by find_RZ_nearby
!> after a tiny perturbation. Rejects positions too close to element boundaries
!> that would immediately fail during particle pushing.
subroutine verify_nearest_position(node_list, element_list, R, Z, s, t, i_elm, ifail)
  use mod_find_rz_nearby, only: find_RZ_nearby
  type(type_node_list),    intent(in)    :: node_list
  type(type_element_list), intent(in)    :: element_list
  real*8,                  intent(in)    :: R, Z
  real*8,                  intent(inout) :: s, t
  integer,                 intent(inout) :: i_elm
  integer,                 intent(out)   :: ifail

  real*8 :: R2, Z2, s2, t2, eps
  integer :: i_elm2

  ! Tiny perturbation (~1 micron) to test robustness
  eps = 1.d-6
  R2 = R + eps
  Z2 = Z + eps
  call find_RZ_nearby(node_list, element_list, R, Z, s, t, i_elm, &
                      R2, Z2, s2, t2, i_elm2, ifail)
end subroutine verify_nearest_position

!> Compute total main ion count = central_density * 1d20 * ∫ rho dV
!> using element-centre approximation for the volume integral
function compute_total_ions(node_list, element_list) result(total_ions)
  type(type_node_list),    intent(in) :: node_list
  type(type_element_list), intent(in) :: element_list
  real*8 :: total_ions

  integer :: i_elm, iv, ierr
  real*8  :: R_v(4), Z_v(4), R_c, area, rho_c, local_sum

  local_sum = 0.d0
  do i_elm = 1, element_list%n_elements
    do iv = 1, 4
      R_v(iv) = node_list%node(element_list%element(i_elm)%vertex(iv))%x(1,1,1)
      Z_v(iv) = node_list%node(element_list%element(i_elm)%vertex(iv))%x(1,1,2)
    end do
    R_c  = 0.25d0 * sum(R_v)
    area = 0.5d0 * abs((R_v(3)-R_v(1))*(Z_v(4)-Z_v(2)) - (R_v(4)-R_v(2))*(Z_v(3)-Z_v(1)))
    rho_c = 0.d0
    do iv = 1, 4
      rho_c = rho_c + node_list%node(element_list%element(i_elm)%vertex(iv))%values(1,1,var_rho)
    end do
    rho_c = 0.25d0 * rho_c
    local_sum = local_sum + 6.283185307179586d0 * R_c * area * rho_c  ! 2*pi*R*area*rho
  end do

  !> node_list/element_list are replicated on all MPI ranks, no reduction needed
  total_ions = dble(n_period) * local_sum * central_density * 1.d20
end function compute_total_ions

end module initialisers_ICS
