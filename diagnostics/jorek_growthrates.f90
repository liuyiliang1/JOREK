!> JOREK Growth Rate Analyzer — Toroidal Mode Stability from Energy Principle
!!
!! Reads JOREK restart files and computes growth rates of toroidal modes
!! using the magnetic and kinetic energy time evolution (energy principle).
!!
!! Growth rate formula:  gamma[n] = 0.5 * d(log|E[n]|)/dt
!!   where E[n] is the magnetic or kinetic energy of toroidal harmonic n.
!!
!! Two modes of operation:
!!   Multi-file: Read a series of restart files, call energy() for each.
!!   Single-file: Read 'energies' and 'xtime' directly from HDF5 restart.
!!
!! Usage:
!!   jorek_growthrates --first 100 --last 200              ! steps 100-200
!!   jorek_growthrates --first 100 --nstep 50              ! 50 steps from 100
!!   jorek_growthrates --first 100 --last 200 --step 5     ! every 5th step
!!   jorek_growthrates --file restart.h5                   ! single file (energies)
!!   jorek_growthrates -h                                  ! help
!!
!! Output:
!!   - Formatted table of growth rates per toroidal mode
!!   - Most unstable mode identification
!!   - Optionally: energies_growth_rates.dat file
!!
!! References:
!!   - energy.f90 (per-toroidal-harmonic magnetic/kinetic energy)
!!   - recalc_egr.f90 (existing but hard-coded growth rate tool)
!!   - out_save_module.f90:152-161 (runtime growth rate computation)
!!   - mod_jorek_timestepping.f90:445-465 (energies storage during run)
program jorek_growthrates

  use hdf5_io_module
  use data_structure
  use mod_parameters, only: n_var, n_tor, n_period, n_degrees, &
                              n_vertex_max
  use phys_module, only: tstep, xtime, energies, mode, &
                          var_psi, var_u, var_rho, F0, &
                          index_now, n_coord_period, rst_format
  use mod_import_restart, only: import_restart
  use mod_chi, only: init_chi_basis
  use basis_at_gaussian, only: initialise_basis

  implicit none

  ! ===== Data structures =====
  type(type_node_list),    pointer :: node_list
  type(type_element_list), pointer :: element_list

  ! ===== Energy and time arrays =====
  real*8, allocatable :: W_mag(:), W_kin(:)        ! per-toroidal-harmonic
  real*8, allocatable :: times(:)                   ! time at each step
  real*8, allocatable :: mag_hist(:,:)              ! (n_tor, n_files) magnetic
  real*8, allocatable :: kin_hist(:,:)              ! (n_tor, n_files) kinetic

  ! ===== CLI arguments =====
  integer              :: narg, iarg
  character(len=256)   :: arg
  integer              :: first_step, last_step, step_inc, nstep_max
  character(len=256)   :: single_filename
  logical              :: single_mode, multi_mode, show_help, save_file

  ! ===== Loop / file I/O =====
  integer              :: istep, itor, ierr, n_files, ios, iunit
  character(len=14)    :: filein
  integer              :: rst_format_val
  character(len=6)     :: numfmt
  real*8               :: dt_val, gamma_mag, gamma_kin, max_gamma
  integer              :: most_unstable_n
  character(len=8)     :: mode_label

  ! ===== Output format =====
  character(len=120)   :: sep_line, header_line
  sep_line = '------------------------------------------------------------------'
  
  ! =====================================================================
  ! Step 1: Parse command-line arguments
  ! =====================================================================
  narg = command_argument_count()
  rst_format_val = rst_format   ! use phys_module default
  single_mode = .false.
  multi_mode  = .false.
  show_help   = .false.
  save_file   = .false.
  first_step  = 0
  last_step   = 0
  step_inc    = 1
  nstep_max   = 0
  single_filename = ''

  if (narg == 0) show_help = .true.

  iarg = 1
  do while (iarg <= narg)
    call get_command_argument(iarg, arg)
    select case (adjustl(trim(arg)))
    case ('-h', '--help')
      show_help = .true.
    case ('--first')
      iarg = iarg + 1
      if (iarg <= narg) then
        call get_command_argument(iarg, arg)
        read(arg, *, iostat=ios) first_step
        if (ios /= 0) stop 'Error: invalid --first value'
        multi_mode = .true.
      end if
    case ('--last')
      iarg = iarg + 1
      if (iarg <= narg) then
        call get_command_argument(iarg, arg)
        read(arg, *, iostat=ios) last_step
        if (ios /= 0) stop 'Error: invalid --last value'
        multi_mode = .true.
      end if
    case ('--step')
      iarg = iarg + 1
      if (iarg <= narg) then
        call get_command_argument(iarg, arg)
        read(arg, *, iostat=ios) step_inc
        if (ios /= 0) step_inc = 1
      end if
    case ('--nstep')
      iarg = iarg + 1
      if (iarg <= narg) then
        call get_command_argument(iarg, arg)
        read(arg, *, iostat=ios) nstep_max
        if (ios /= 0) stop 'Error: invalid --nstep value'
        multi_mode = .true.
      end if
    case ('--file')
      iarg = iarg + 1
      if (iarg <= narg) then
        call get_command_argument(iarg, arg)
        single_filename = adjustl(trim(arg))
        single_mode = .true.
      end if
    case ('--save')
      save_file = .true.
    case ('--format')
      iarg = iarg + 1
      if (iarg <= narg) then
        call get_command_argument(iarg, arg)
        select case (adjustl(trim(arg)))
        case ('HDF5', 'hdf5', '1')
          rst_format_val = 1
        case ('bin', 'binary', '0')
          rst_format_val = 0
        end select
      end if
    case default
      ! Assume it's a filename (positional)
      if (single_filename == '') then
        single_filename = adjustl(trim(arg))
        single_mode = .true.
      end if
    end select
    iarg = iarg + 1
  end do

  if (show_help) then
    call print_usage()
    stop
  end if

  if (.not. single_mode .and. .not. multi_mode) then
    write(*,*) 'Error: specify either --file or --first/--last/--nstep'
    call print_usage()
    stop 1
  end if

  ! =====================================================================
  ! Step 2: Initialize (det_modes → basis → chi → parameters)
  ! =====================================================================
  write(*,*) '========================================'
  write(*,*) '  JOREK Growth Rate Analyzer'
  write(*,*) '========================================'
  write(*,*)

  ! CRITICAL: det_modes must be called first (reviewer finding #1)
  call det_modes

  ! Initialize basis functions (must complete before energy() calls)
  call initialise_basis
  call init_chi_basis

  ! Initialize parameters (skip namelist file)
  call initialise_parameters(0, '__NO_FILENAME__')

  ! =====================================================================
  ! Step 3: Allocate data structures
  ! =====================================================================
  allocate(node_list)
  allocate(element_list)

  allocate(W_mag(n_tor))
  allocate(W_kin(n_tor))

  ! =====================================================================
  ! Step 4: Process data and compute growth rates
  ! =====================================================================
  if (single_mode) then
    call process_single_file()
  else
    call process_multi_files()
  end if

  ! =====================================================================
  ! Step 5: Cleanup
  ! =====================================================================
  if (allocated(W_mag))      deallocate(W_mag)
  if (allocated(W_kin))      deallocate(W_kin)
  if (allocated(times))      deallocate(times)
  if (allocated(mag_hist))   deallocate(mag_hist)
  if (allocated(kin_hist))   deallocate(kin_hist)
  if (associated(node_list))    deallocate(node_list)
  if (associated(element_list)) deallocate(element_list)

contains

  ! ------------------------------------------------------------------
  ! Print usage information
  ! ------------------------------------------------------------------
  subroutine print_usage()
    write(*,*) 'JOREK Growth Rate Analyzer'
    write(*,*) ''
    write(*,*) 'Usage (multi-file mode):'
    write(*,*) '  jorek_growthrates --first 100 --last 200'
    write(*,*) '  jorek_growthrates --first 100 --nstep 50'
    write(*,*) '  jorek_growthrates --first 100 --last 200 --step 5'
    write(*,*) '  jorek_growthrates --first 100 --last 200 --save'
    write(*,*) ''
    write(*,*) 'Usage (single-file mode):'
    write(*,*) '  jorek_growthrates --file restart.h5'
    write(*,*) '  jorek_growthrates restart.h5'
    write(*,*) ''
    write(*,*) 'Options:'
    write(*,*) '  --first N     First restart step number'
    write(*,*) '  --last N      Last restart step number'
    write(*,*) '  --step N      Step increment (default: 1)'
    write(*,*) '  --nstep N     Number of steps from --first'
    write(*,*) '  --file NAME   Single restart HDF5 file'
    write(*,*) '  --format FMT  Restart format (HDF5 or bin, default: HDF5)'
    write(*,*) '  --save        Save results to energies_growth_rates.dat'
    write(*,*) '  -h, --help    This help message'
    write(*,*) ''
    write(*,*) 'Output:'
    write(*,*) '  Growth rates per toroidal harmonic (magnetic and kinetic)'
    write(*,*) '  Identifies the most unstable mode'
  end subroutine print_usage

  ! ------------------------------------------------------------------
  ! Single-file mode: read energies/xtime from HDF5
  ! ------------------------------------------------------------------
  subroutine process_single_file()
#ifdef USE_HDF5
    use hdf5
    integer(HID_T) :: file_id
    integer        :: ierr_h5, n_times, n_tor_file
    integer        :: ed1, ed2, ed3      ! temp dims for HDF5 reads
    real*8, allocatable :: t_energies(:,:,:)
    real*8, allocatable :: t_xtime(:)

    write(*,'(A,A)') 'Reading single restart file: ', trim(single_filename)

    ! Open HDF5 file
    call H5Fopen_f(trim(single_filename)//char(0), H5F_ACC_RDONLY_F, &
                    file_id, ierr_h5)
    if (ierr_h5 /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(single_filename)
      stop 1
    end if

    ! Read scalar metadata for verification
    n_tor_file = h5_read_integer(file_id, 'n_tor', n_tor)
    write(*,'(A,I4)') '  n_tor from file: ', n_tor_file

    ! Read xtime (1D: time at each save index)
    call h5_read_1d_double_alloc(file_id, 'xtime', t_xtime, n_times)
    if (.not. allocated(t_xtime)) then
      write(*,*) 'ERROR: No xtime dataset found in restart file.'
      write(*,*) '  Try multi-file mode: --first N --last M'
      call H5Fclose_f(file_id, ierr_h5)
      stop 1
    end if
    write(*,'(A,I6)') '  Time steps found: ', n_times
    write(*,'(A,ES12.4,A,ES12.4)') '  Time range: [', &
         t_xtime(1), ', ', t_xtime(n_times), ']'

    if (n_times < 2) then
      write(*,*) 'ERROR: Need at least 2 time steps. Only found ', n_times
      call H5Fclose_f(file_id, ierr_h5)
      stop 1
    end if

    ! Read energies (3D: n_tor x 2 x n_times)
    call h5_read_3d_double_alloc(file_id, 'energies', t_energies, ed1, ed2, ed3)
    if (.not. allocated(t_energies)) then
      write(*,*) 'ERROR: No energies dataset found.'
      write(*,*) '  JOREK must be compiled with energy tracking enabled.'
      call H5Fclose_f(file_id, ierr_h5)
      stop 1
    end if
    if (ed1 /= n_tor .or. ed2 /= 2) then
      write(*,'(A,3I6)') '  WARNING: energies shape mismatch, got ', ed1, ed2, ed3
    end if

    call H5Fclose_f(file_id, ierr_h5)

    ! Copy to module-level arrays
    allocate(times(n_times))
    times = t_xtime

    allocate(mag_hist(n_tor_file, n_times))
    allocate(kin_hist(n_tor_file, n_times))
    mag_hist = t_energies(:, 1, :)
    kin_hist = t_energies(:, 2, :)

    deallocate(t_xtime)
    deallocate(t_energies)

    n_files = n_times

    ! Compute and report growth rates
    call compute_and_report()
#else
    write(*,*) 'ERROR: JOREK was compiled without HDF5 support.'
    write(*,*) '  Single-file mode requires USE_HDF5.'
    write(*,*) '  Use multi-file mode instead: --first N --last M'
    stop 1
#endif
  end subroutine process_single_file

  ! ------------------------------------------------------------------
  ! Multi-file mode: read restart files, call energy() for each
  ! ------------------------------------------------------------------
  subroutine process_multi_files()
    integer :: i, step_val, n_tot
    logical :: file_exists

    ! Resolve step range
    if (last_step == 0 .and. nstep_max > 0) then
      last_step = first_step + (nstep_max - 1) * step_inc
    end if
    if (last_step == 0) then
      write(*,*) 'ERROR: specify --last or --nstep'
      stop 1
    end if

    n_tot = (last_step - first_step) / step_inc + 1
    if (n_tot <= 1) then
      write(*,*) 'ERROR: Need at least 2 restart files.'
      stop 1
    end if

    write(*,'(A,I6,A,I6,A,I3)') 'Multi-file mode: ', first_step, &
         ' -> ', last_step, '  step=', step_inc
    write(*,'(A,I6)') '  Expected files: ', n_tot

    ! Allocate history arrays
    allocate(times(n_tot))
    allocate(mag_hist(n_tor, n_tot))
    allocate(kin_hist(n_tor, n_tot))
    times     = 0.d0
    mag_hist  = 0.d0
    kin_hist  = 0.d0

    ! Determine number format for restart filenames
    if (last_step >= 100000) then
      numfmt = 'i6.6'
    else
      numfmt = 'i6.6'
    end if

    n_files = 0
    do i = 1, n_tot
      step_val = first_step + (i - 1) * step_inc

      ! Build restart filename
      write(filein, '(''jorek'',' // numfmt // ')') step_val

      ! Check if file exists
      inquire(file=trim(filein)//'.h5', exist=file_exists)
      if (.not. file_exists) then
        inquire(file=trim(filein)//'.rst', exist=file_exists)
      end if
      if (.not. file_exists) then
        write(*,'(A,I6,A)') '  WARNING: step ', step_val, ' not found, skipping'
        cycle
      end if

      ! Import restart
      call import_restart(node_list, element_list, filein, &
                           rst_format_val, ierr, .false.)

      if (ierr /= 0 .and. ierr /= -1) then
        write(*,'(A,I6,A,I4)') '  WARNING: step ', step_val, &
             ' import error=', ierr
        cycle
      end if

      n_files = n_files + 1

      ! Compute energy for this step
      W_mag = 0.d0
      W_kin = 0.d0
      call energy(W_mag, W_kin)

      ! Store
      times(n_files) = xtime(index_now)
      mag_hist(:, n_files) = W_mag(:)
      kin_hist(:, n_files) = W_kin(:)

      write(*,'(A,I6,A,ES12.4)') '  [', step_val, '] t=', times(n_files)
    end do

    if (n_files < 2) then
      write(*,*) 'ERROR: Only ', n_files, ' files successfully read.'
      write(*,*) '  Need at least 2 for growth rate analysis.'
      stop 1
    end if

    write(*,*)
    write(*,'(A,I6,A)') 'Successfully processed ', n_files, ' files.'

    ! Truncate arrays to actual number of files read
    ! (Keep the allocations as-is; n_files tracks the valid count)

    call compute_and_report()

  end subroutine process_multi_files

  ! ------------------------------------------------------------------
  ! Compute growth rates from energy history and report
  ! ------------------------------------------------------------------
  subroutine compute_and_report()
    integer :: i, it, n_used
    real*8  :: gm, gk, em_i, ek_i

    ! Use only the valid entries
    n_used = n_files

    write(*,*)
    write(*,'(A)') sep_line
    write(*,'(A)') '  Growth Rate Analysis Results'
    write(*,'(A)') sep_line
    write(*,'(A,I4,A,I4,A,I4)') '  n_tor: ', n_tor, &
         '  n_period: ', n_period, '  n_steps: ', n_used
    write(*,'(A,ES12.4,A,ES12.4)') '  Time range: [', times(1), &
         ', ', times(n_used), ']'
    write(*,'(A)') sep_line
    write(*,'(A8,A3,A14,A3,A14,A3,A14,A3,A14)') &
         '  Mode', ' |', '  gamma_mag   ', ' |', '  gamma_kin   ', &
         ' |', '  E_mag(last) ', ' |', '  E_kin(last) '
    write(*,'(A)') sep_line

    most_unstable_n = 0
    max_gamma = -1.d30

    do it = 1, n_tor
      ! Compute growth rate from log-derivative
      ! gamma = 0.5 * (log|E(t_last)| - log|E(t_first)|) / (t_last - t_first)
      ! Use linear regression for robustness (simple two-point here for
      ! clarity; multi-point average across all intervals for better stats)

      ! --- Magnetic growth rate ---
      gm = 0.d0
      if (abs(mag_hist(it, n_used)) > 1.d-40 .and. &
          abs(mag_hist(it, 1)) > 1.d-40) then
        dt_val = times(n_used) - times(1)
        if (dt_val > 0.d0) then
          gm = 0.5d0 * log(abs(mag_hist(it, n_used) / &
                               mag_hist(it, 1))) / dt_val
        end if
      end if

      ! Try multi-point average if single-point is zero/noisy
      if (abs(gm) < 1.d-40 .and. n_used > 2) then
        gm = fit_growth_rate(times(1:n_used), mag_hist(it, 1:n_used))
      end if

      ! --- Kinetic growth rate ---
      gk = 0.d0
      if (abs(kin_hist(it, n_used)) > 1.d-40 .and. &
          abs(kin_hist(it, 1)) > 1.d-40) then
        dt_val = times(n_used) - times(1)
        if (dt_val > 0.d0) then
          gk = 0.5d0 * log(abs(kin_hist(it, n_used) / &
                               kin_hist(it, 1))) / dt_val
        end if
      end if

      if (abs(gk) < 1.d-40 .and. n_used > 2) then
        gk = fit_growth_rate(times(1:n_used), kin_hist(it, 1:n_used))
      end if

      ! --- Mode label ---
      call get_mode_label(it, mode_label)

      ! --- Final energies ---
      em_i = mag_hist(it, n_used)
      ek_i = kin_hist(it, n_used)

      ! --- Print ---
      if (it == 1) then
        write(*,'(A8,A3,ES14.6,A3,ES14.6,A3,ES14.6,A3,ES14.6)') &
             mode_label, ' |', 0.d0, ' |', 0.d0, ' |', em_i, ' |', ek_i
      else
        write(*,'(A8,A3,ES14.6,A3,ES14.6,A3,ES14.6,A3,ES14.6)') &
             mode_label, ' |', gm, ' |', gk, ' |', em_i, ' |', ek_i
      end if

      ! Track most unstable (n>0 only)
      if (it > 1 .and. gm > max_gamma) then
        max_gamma = gm
        most_unstable_n = mode(it)
      end if
    end do

    write(*,'(A)') sep_line
    if (most_unstable_n > 0) then
      write(*,'(A,I3,A,ES12.4)') '  Most unstable: n=', &
           most_unstable_n, '  gamma_mag=', max_gamma
      if (max_gamma > 0.d0) then
        write(*,'(A,ES12.4,A)') '  Growth time:   tau=', &
             1.d0 / max_gamma, ' s'
      end if
    else
      write(*,'(A)') '  No unstable modes detected (all gamma <= 0)'
    end if
    write(*,'(A)') sep_line

    ! --- Optionally save to file ---
    if (save_file) then
      open(newunit=iunit, file='energies_growth_rates.dat', &
           action='write', status='replace')
      write(iunit, '(A)') '# JOREK Growth Rate Analysis'
      write(iunit, '(A,I4,A,I4)') '# n_tor=', n_tor, ' n_period=', n_period
      write(iunit, '(A)') '#   itor  n_mode       time_start       time_end', &
           '    gamma_mag       gamma_kin      E_mag_final      E_kin_final'
      do it = 1, n_tor
        call get_mode_label(it, mode_label)
        gm = 0.d0; gk = 0.d0
        if (it > 1 .and. n_used > 1) then
          dt_val = times(n_used) - times(1)
          if (dt_val > 0.d0 .and. abs(mag_hist(it, n_used)) > 1.d-40) &
               gm = 0.5d0 * log(abs(mag_hist(it, n_used) / &
               max(abs(mag_hist(it, 1)), 1.d-40))) / dt_val
          if (dt_val > 0.d0 .and. abs(kin_hist(it, n_used)) > 1.d-40) &
               gk = 0.5d0 * log(abs(kin_hist(it, n_used) / &
               max(abs(kin_hist(it, 1)), 1.d-40))) / dt_val
        end if
        write(iunit, '(I6,A10,4ES15.6,2ES16.8)') it, adjustr(mode_label), &
             times(1), times(n_used), gm, gk, &
             mag_hist(it, n_used), kin_hist(it, n_used)
      end do
      close(iunit)
      write(*,*) 'Results saved to: energies_growth_rates.dat'
    end if

  end subroutine compute_and_report

  ! ------------------------------------------------------------------
  ! Get human-readable mode label from toroidal index
  ! ------------------------------------------------------------------
  subroutine get_mode_label(itor, label)
    integer,          intent(in)  :: itor
    character(len=8), intent(out) :: label

    if (n_period <= 1) then
      ! Tokamak: index 1 = n=0, 2=n=1 cos, 3=n=1 sin, 4=n=2 cos, ...
      if (itor == 1) then
        write(label, '(A)') 'n=0'
      else if (mod(itor, 2) == 0) then
        write(label, '(A,I1,A)') 'n=', (itor / 2), '(c)'
      else
        write(label, '(A,I1,A)') 'n=', (itor / 2), '(s)'
      end if
    else
      ! Stellarator: mode family based
      write(label, '(A,I3)') 'mf=', mode(itor)
    end if
  end subroutine get_mode_label

  ! ------------------------------------------------------------------
  ! Fit growth rate using linear regression on log(|E|) vs time
  ! ------------------------------------------------------------------
  real*8 function fit_growth_rate(t, e)
    real*8, intent(in) :: t(:), e(:)
    integer :: n, i, n_valid
    real*8  :: t_mean, slope, s_tt, s_te
    real*8, allocatable :: log_e(:), t_work(:)

    n = size(t)
    if (n < 2) then
      fit_growth_rate = 0.d0
      return
    end if

    ! Filter: only use points with positive energy
    n_valid = 0
    do i = 1, n
      if (abs(e(i)) > 1.d-40) n_valid = n_valid + 1
    end do
    if (n_valid < 2) then
      fit_growth_rate = 0.d0
      return
    end if

    allocate(log_e(n_valid), t_work(n_valid))
    n_valid = 0
    do i = 1, n
      if (abs(e(i)) > 1.d-40) then
        n_valid = n_valid + 1
        t_work(n_valid) = t(i)
        log_e(n_valid) = log(abs(e(i)))
      end if
    end do

    ! Linear regression: logE = a + 2*gamma * (t - t_mean)
    t_mean = sum(t_work) / dble(n_valid)
    s_tt = sum((t_work - t_mean)**2)
    s_te = sum((t_work - t_mean) * (log_e - sum(log_e) / dble(n_valid)))

    if (s_tt > 0.d0) then
      slope = s_te / s_tt
      fit_growth_rate = slope * 0.5d0   ! gamma = slope/2
    else
      fit_growth_rate = 0.d0
    end if

    deallocate(log_e, t_work)
  end function fit_growth_rate

  ! ------------------------------------------------------------------
  ! HDF5 helper: read scalar integer
  ! ------------------------------------------------------------------
#ifdef USE_HDF5
  integer function h5_read_integer(file_id, dset_name, default_val)
    use hdf5
    integer(HID_T),   intent(in) :: file_id
    character(len=*), intent(in) :: dset_name
    integer,          intent(in) :: default_val

    integer(HID_T) :: dset_id
    integer        :: ierr, val
    logical        :: exists

    h5_read_integer = default_val

    call H5Lexists_f(file_id, trim(dset_name)//char(0), exists, ierr)
    if (.not. exists) return

    call H5Dopen_f(file_id, trim(dset_name)//char(0), dset_id, ierr)
    if (ierr /= 0) return

    call H5Dread_f(dset_id, H5T_NATIVE_INTEGER, val, &
                    [integer(HID_T) :: 0], ierr)
    if (ierr == 0) h5_read_integer = val

    call H5Dclose_f(dset_id, ierr)
  end function h5_read_integer

  ! ------------------------------------------------------------------
  ! HDF5 helper: read 1D double array (allocatable)
  ! ------------------------------------------------------------------
  subroutine h5_read_1d_double_alloc(file_id, dset_name, arr, n)
    use hdf5
    integer(HID_T),                  intent(in)  :: file_id
    character(len=*),                intent(in)  :: dset_name
    real*8, allocatable,             intent(out) :: arr(:)
    integer,                         intent(out) :: n

    integer(HID_T) :: dset_id, dspace_id
    integer(HSIZE_T), allocatable :: dims(:)
    integer :: ierr, ndims
    logical :: exists

    n = 0
    if (allocated(arr)) deallocate(arr)

    call H5Lexists_f(file_id, trim(dset_name)//char(0), exists, ierr)
    if (.not. exists) return

    call H5Dopen_f(file_id, trim(dset_name)//char(0), dset_id, ierr)
    if (ierr /= 0) return

    call H5Dget_space_f(dset_id, dspace_id, ierr)
    call H5Sget_simple_extent_ndims_f(dspace_id, ndims, ierr)
    allocate(dims(ndims))
    call H5Sget_simple_extent_dims_f(dspace_id, dims, dims, ierr)

    n = int(dims(1))
    allocate(arr(n))
    call H5Dread_f(dset_id, H5T_NATIVE_DOUBLE, arr, dims, ierr)

    call H5Sclose_f(dspace_id, ierr)
    call H5Dclose_f(dset_id, ierr)
    deallocate(dims)
  end subroutine h5_read_1d_double_alloc

  ! ------------------------------------------------------------------
  ! HDF5 helper: read 3D double array (allocatable)
  ! ------------------------------------------------------------------
  subroutine h5_read_3d_double_alloc(file_id, dset_name, arr, d1, d2, d3)
    use hdf5
    integer(HID_T),                  intent(in)  :: file_id
    character(len=*),                intent(in)  :: dset_name
    real*8, allocatable,             intent(out) :: arr(:,:,:)
    integer,                         intent(out) :: d1, d2, d3

    integer(HID_T) :: dset_id, dspace_id
    integer(HSIZE_T), allocatable :: dims(:), maxdims(:)
    integer :: ierr, ndims
    logical :: exists

    d1 = 0; d2 = 0; d3 = 0
    if (allocated(arr)) deallocate(arr)

    call H5Lexists_f(file_id, trim(dset_name)//char(0), exists, ierr)
    if (.not. exists) return

    call H5Dopen_f(file_id, trim(dset_name)//char(0), dset_id, ierr)
    if (ierr /= 0) return

    call H5Dget_space_f(dset_id, dspace_id, ierr)
    call H5Sget_simple_extent_ndims_f(dspace_id, ndims, ierr)
    allocate(dims(ndims), maxdims(ndims))
    call H5Sget_simple_extent_dims_f(dspace_id, dims, maxdims, ierr)

    d1 = int(dims(1)); d2 = int(dims(2)); d3 = int(dims(3))
    allocate(arr(d1, d2, d3))
    call H5Dread_f(dset_id, H5T_NATIVE_DOUBLE, arr, dims, ierr)

    call H5Sclose_f(dspace_id, ierr)
    call H5Dclose_f(dset_id, ierr)
    deallocate(dims, maxdims)
  end subroutine h5_read_3d_double_alloc
#endif

end program jorek_growthrates
