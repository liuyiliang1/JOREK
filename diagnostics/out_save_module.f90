!-------------------------------------------------
! file : output_saving.f90
! date : 19/01/2012
!  used for saving the physics results in HDF5 format
!------------------------------------------------
module out_save_module
  use data_structure
  use phys_module
  use basis_at_gaussian
#ifdef USE_HDF5
  use hdf5_io_module

  implicit none
  
  !******************************
  contains
  !******************************
      

  !***********************************************************
  !   SAVING DIRECTLY IN HDF5 FORMAT
  !***********************************************************

  !----------------------------------------------
  ! Output saving of the basic model quantities
  !----------------------------------------------
  subroutine HDF5_basics_save(Sindex_now,St_now)
#ifdef USE_HDF5
    use hdf5
#endif
    use mod_parameters,  only: jorek_model, n_tor, n_plane, n_period
    use phys_module, only: tstep, heatsource, particlesource, &
         eta, eta_num, visco, visco_num, visco_par, visco_par_num, &
         tauIC, gamma_sheath, numfmt
    !type(type_node_list)   , intent(in) :: node_list
    !type(type_element_list), intent(in) :: element_list
    integer                , intent(in) :: Sindex_now
    real*8                 , intent(in) :: St_now

#ifdef USE_HDF5
    character(LEN=100) :: basics_file
    integer(HID_T)     :: file_id
    integer            :: i, ierr
    integer            :: Nbr, Nbc

    if (Sindex_now.ne.-1) then
       write(basics_file,'(A,'//numfmt//',A)') &
            "basics", Sindex_now, ".h5"
    else
       basics_file = "basics_init.h5"
    end if
    if (pglobal_id == 0) then

    call HDF5_create(trim(basics_file),file_id,ierr)
    if (ierr.ne.0) then
       print*,'pglobal_id = ',pglobal_id, &
            ' ==> error for opening of ',basics_file
    end if

    Nbr = 1
    Nbc = 1
    !*** the basic model parameters ***
    ! -> save : 'jorek_model   '
    call HDF5_integer_saving(file_id,jorek_model,'BASIC_jorek_model'//char(0))
    ! -> save : 'n_tor   '
    call HDF5_integer_saving(file_id,n_tor,'BASIC_n_tor'//char(0))
    ! -> save : 'n_plane   '
    call HDF5_integer_saving(file_id,n_plane,'BASIC_n_plane'//char(0))
    ! -> save : 'n_period   '
    call HDF5_integer_saving(file_id,n_period,'BASIC_n_period'//char(0))

    !*** the source, sink, dissipation, boundary condition, etc. parameters ***
    ! -> save : 'time_diag   '
    call HDF5_real_saving(file_id,St_now,'time_diag'//char(0))
    ! -> save : 'tstep   '
    call HDF5_real_saving(file_id,tstep,'BASIC_deltat_diag'//char(0))
    ! -> save : 'iteration_global   '
    call HDF5_integer_saving(file_id,Sindex_now,'BASIC_it_diag'//char(0))
    ! -> save : 'heatsource   '
    call HDF5_real_saving(file_id,heatsource,'BASIC_heatsource'//char(0))
    ! -> save : 'particlesource   '
    call HDF5_real_saving(file_id,particlesource,'BASIC_particlesource'//char(0))
    ! -> save : 'eta   '
    call HDF5_real_saving(file_id,eta,'BASIC_eta'//char(0))
    ! -> save : 'eta_num   '
    call HDF5_real_saving(file_id,eta_num,'BASIC_eta_num'//char(0))
    ! -> save : 'visco   '
    call HDF5_real_saving(file_id,visco,'BASIC_visco'//char(0))
    ! -> save : 'visco_num   '
    call HDF5_real_saving(file_id,visco_num,'BASIC_visco_num'//char(0))
    ! -> save : 'visco_par   '
    call HDF5_real_saving(file_id,visco_par,'BASIC_visco_par'//char(0))
    ! -> save : 'visco_par_num   '
    call HDF5_real_saving(file_id,visco_par_num,'BASIC_visco_par_num'//char(0))
    ! -> save : 'tauIC   '
    call HDF5_real_saving(file_id,tauIC,'BASIC_tauIC'//char(0))
    ! -> save : 'gamma_sheath   '
    call HDF5_real_saving(file_id,gamma_sheath,'BASIC_gamma_sheath'//char(0))

    call HDF5_close(file_id)

    end if  ! (pglobal_id == 0)
#else
    print*,' ==> no savings of the equilibrium HDF5 files, check the USE_HDF5 option'
#endif

  end subroutine HDF5_basics_save

  !----------------------------------------------
  ! Output savings for different n_tor values
  !----------------------------------------------
  subroutine hdf5_ntor_profiles_save(Sindex_now)
#ifdef USE_HDF5
    use hdf5
#endif
    use mod_parameters,  only: n_tor
    use phys_module, only: tstep, index_start, energies
    !type(type_node_list)   , intent(in) :: node_list
    !type(type_element_list), intent(in) :: element_list
    integer                , intent(in) :: Sindex_now
    
#ifdef USE_HDF5
    character(LEN=100) :: ntor_profiles_file
    integer(HID_T)     :: file_id
    integer            :: i, ierr
    integer            :: Nbr, Nbc
    real*8             :: Emag_ntor_glob(n_tor), Ekin_ntor_glob(n_tor), &
         Growthrate_mag_ntor_glob(n_tor), Growthrate_kin_ntor_glob(n_tor)

    if (Sindex_now.ne.-1) then
       write(ntor_profiles_file,'(A,'//numfmt//',A)') &
            "ntor_profiles", Sindex_now, ".h5"
    else
       ntor_profiles_file = "ntor_profiles_init.h5"
    end if
    if (pglobal_id == 0) then

    call HDF5_create(trim(ntor_profiles_file),file_id,ierr)
    if (ierr.ne.0) then
       print*,'pglobal_id = ',pglobal_id, &
            ' ==> error for opening of ',ntor_profiles_file
    end if

    !*** Compute energies(a, b, c), where:      ***
    !***   a = 1:n_tor  ==> for each n harmonic ***
    !***   b = (1, magnetic) or (2, kinetic)    ***
    !***   c = 1:index_now ==> all times        ***
    Emag_ntor_glob          = 0.d0
    Emag_ntor_glob(1:n_tor) = energies(1:n_tor,1,Sindex_now)
    Ekin_ntor_glob          = 0.d0
    Ekin_ntor_glob(1:n_tor) = energies(1:n_tor,2,Sindex_now)
    !*** Compute growth rates            ***
    !***   computed from energies(a,b,c) ***
    Growthrate_mag_ntor_glob          = 0.d0
    Growthrate_kin_ntor_glob          = 0.d0
    if (Sindex_now > index_start+1) then
       do i=1,n_tor
          Growthrate_mag_ntor_glob(i) = 0.5d0*log(abs(energies(i,1,Sindex_now)/energies(i,1,Sindex_now-1))) / tstep
          Growthrate_kin_ntor_glob(i) = 0.5d0*log(abs(energies(i,2,Sindex_now)/energies(i,2,Sindex_now-1))) / tstep
       end do
    end if

    ! -> save : 'kin. and mag. energies   '
    Nbc = n_tor
    call HDF5_array1D_saving(file_id,Emag_ntor_glob(1:),Nbc, &
         'Emag_ntor_glob'//char(0))
    call HDF5_array1D_saving(file_id,Ekin_ntor_glob(1:),Nbc, &
         'Ekin_ntor_glob'//char(0))
    ! -> save : 'kin. and mag. growth rates   '
    call HDF5_array1D_saving(file_id,Growthrate_mag_ntor_glob(1:),Nbc, &
         'Growthrate_mag_ntor_glob'//char(0))
    call HDF5_array1D_saving(file_id,Growthrate_kin_ntor_glob(1:),Nbc, &
         'Growthrate_kin_ntor_glob'//char(0))

    end if  ! (pglobal_id == 0)
#else
    print*,' ==> no savings of the ntor profiles HDF5 files, check the USE_HDF5 option'
#endif

  end subroutine hdf5_ntor_profiles_save

  !----------------------------------------------
  ! Output saving of the radial profiles
  !----------------------------------------------
  subroutine HDF5_radial_profiles_save(node_list,element_list,Sindex_now,St_now)
#ifdef USE_HDF5
    use hdf5
#endif
    type(type_node_list)   , intent(in) :: node_list
    type(type_element_list), intent(in) :: element_list
    integer                , intent(in) :: Sindex_now
    real*8                 , intent(in) :: St_now
    
#ifdef USE_HDF5
    character(LEN=100) :: radial_profiles_file
    integer(HID_T)     :: file_id
    integer            :: i, ierr
    integer            :: Nbr, Nbc


    if (Sindex_now.ne.-1) then
       write(radial_profiles_file,'(A,'//numfmt//',A)') &
            "radial_profiles", Sindex_now, ".h5"
    else
       radial_profiles_file = "radial_profiles_init.h5"
    end if
    if (pglobal_id == 0) then

    call HDF5_create(trim(radial_profiles_file),file_id,ierr)
    if (ierr.ne.0) then
       print*,'pglobal_id = ',pglobal_id, &
            ' ==> error for opening of ',radial_profiles_file
    end if

    !*** geometry ***
    end if  ! (pglobal_id == 0)
#else
    print*,' ==> no savings of the radial profiles HDF5 files, check the USE_HDF5 option'
#endif

  end subroutine HDF5_radial_profiles_save

#endif
end module out_save_module
