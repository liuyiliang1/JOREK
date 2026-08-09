      program eqdsk
!--------------------------------------------------------------------
! little program to construct an input file for jorek out of data
! in a eqdsk file
!                         Guido Huysmans,          date : 14-12-2010
!
! Some documentation can be found here: https://www.jorek.eu/wiki/doku.php?id=eqdsk2jorek.f90
!--------------------------------------------------------------------
implicit none

real*8,allocatable :: psi(:),p(:),f(:),q(:),rlim(:),zlim(:),rbnd(:), zbnd(:)
real*8,allocatable :: apsi(:),bpsi(:),cpsi(:),dpsi(:)
real*8,allocatable :: ap(:),bp(:),cp(:),dp(:)
real*8,allocatable :: af(:),bf(:),cf(:),df(:)
real*8,allocatable :: radius(:),theta(:),rad(:)
real*8,allocatable :: ar(:),br(:),cr(:),dr(:)
real*8,allocatable :: dpr(:),df2(:),dg(:),work(:),psirz(:,:)
real*8,allocatable :: xx(:),yy(:),zc(:), r_bnd(:), z_bnd(:), psi_bnd(:)
real*8,allocatable :: df2_ext(:),rho_ext(:),T_ext(:),psi_ext(:),p_ext(:),T_ext_i(:),T_ext_e(:)
real*8,allocatable :: tx(:),ty(:),c(:,:),wrk(:)
integer,allocatable :: iwrk(:)
real*8             :: angle, ellip, tria_u, tria_l, quad_u, quad_l, r0, z0, a0, PI
real*8             :: dummy(3), xdim,zdim,rzero,rgrid1,zmid,rmaxis,zmaxis,ssimag,ssibry,bcentr
real*8             :: xip,xdum1,xdum2,xdum3,xdum4,xdum5
real*8             :: psi_sep, sig_sep, tanh1, zmu0, zn0, zmd, rho_bnd
real*8             :: sig_sol, pres_sol
real*8             :: xb ,xe, yb, ye, smth, fp, fout
real*8             :: B_scale, I_scale, R_scale, F_axis, factor, dfactor
integer            :: mx,my,kx,ky,nxest,nyest,lwrk,kwrk,ier,iopt,nx,ny, i1, j1
integer            :: nr, nz, n_psi, nbbs, limitr, i,j, nc, n_tht, n_sol, n_ext, ivtk
character          :: AA*52, tokamak_name*50
character          :: buffer*80, lf*1, str1*12, str2*24
real*8             :: coef_Ti_Te                  !coef_Ti_Te=Ti/Te
real*8,allocatable :: ne_spline(:)
integer            :: err_alloc
logical            :: ferr
real*8             :: T_bnd

!----------------------------- read eqdsk file -----------

B_scale = 1.d0/1.d0  ! scaling factor for the vacuum toroidal field 
I_scale = 1.d0/1.d0  ! scaling factor for the toroidal current
R_scale = 1.d0/1.d0  ! scaling factor for the space coordinates 

write(*,*) ' EQDSK to JOREK2 '

tokamak_name = 'EXL-50U' !'ITER' 'DIII-D' 'JET' 'CFETR' 'HL-2M' 'EAST''EXL-50U','EHL2'

write(*,*) '   Tokamak = ', tokamak_name

read(5,'(A52,2i4)') AA,nr,nz

write(*,*) AA
write(*,'(A,2i5)') ' nr, nz : ',nr,nz

!----------------------------- read eqdsk file -----------
read(5,'(5e16.9)') xdim,zdim,rzero,rgrid1,zmid
read(5,'(5e16.9)') rmaxis,zmaxis,ssimag,ssibry,bcentr
read(5,'(5e16.9)') xip,ssimag,xdum1,rmaxis,xdum2
read(5,'(5e16.9)') zmaxis,xdum3,ssibry,xdum4,xdum5

write(*,'(A,3f10.5,A)') ' xdim,rgrid1,zdim : ',xdim,rgrid1,zdim,' m'
write(*,'(A,2f10.5,A)') ' rzero, zmid : ',rzero,zmid, ' m'
write(*,'(A,f10.5,A)')  ' xip         : ',xip/1e6,' MA'
write(*,'(A,f10.5,A)')  ' rmaxis      : ',rmaxis,' m'
write(*,'(A,f10.5,A)')  ' zmaxis      : ',zmaxis,' m'
write(*,'(A,f10.5,A)')  ' Bvac        : ',bcentr,' T'
write(*,'(A,3f10.5,A)') ' psi (axis,bnd) :',ssimag,ssibry,ssibry-ssimag,' Wb'

write(*,*) '   reading profiles'
  
n_psi=nr
allocate(f(n_psi),p(n_psi),df2(n_psi),dpr(n_psi),psirz(nr,nz),q(n_psi))

read(5,'(5e16.9)') (f(i),i=1,n_psi)
read(5,'(5e16.9)') (p(i),i=1,n_psi)
read(5,'(5e16.9)') (df2(i),i=1,n_psi)
read(5,'(5e16.9)') (dpr(i),i=1,n_psi)
read(5,'(5e16.9)') ((psirz(i,j),i=1,nr),j=1,nz)
read(5,'(5e16.9)') (q(i),i=1,n_psi)

allocate(psi(n_psi))
open(21,file='profiles_eqdsk')
do i=1,n_psi
  psi(i) = real(i-1)/real(n_psi-1)
  write(21,'(6e14.6)') psi(i),f(i),p(i),df2(i),dpr(i),q(i)
enddo
close(21)


write(*,*) '   reading limiter'

read(5,*)  nbbs,limitr
allocate(rbnd(nbbs),zbnd(nbbs))
read(5,'(5e16.9)') (rbnd(i),zbnd(i),i=1,nbbs)
allocate(rlim(limitr),zlim(limitr))
read(5,'(5e16.9)') (rlim(i),zlim(i),i=1,limitr)

write(*,*) ' done reading'

!=============== scaling of equilibrium with B
bcentr = bcentr * B_scale
p      = p      * B_scale**2
dpr    = dpr    * B_scale
F      = F      * B_scale
dF2    = dF2    * B_scale 
psirz  = psirz  * B_scale 
xip    = xip    * B_scale

!=============== scaling of equilibrium with space dimension
xip    = xip    * R_scale
p      = p
dpr    = dpr    / R_scale**2
F      = F      * R_scale
dF2    = dF2
psirz  = psirz  * R_scale**2

rgrid1 = rgrid1 * R_scale
rzero  = rzero  * R_scale
xdim   = xdim   * R_scale
zmid   = zmid   * R_scale
zdim   = Zdim   * R_scale

!=============== scaling of equilibrium with plasma current
xip    = xip    * I_scale
p      = p      * I_scale**2
dpr    = dpr    * I_scale
dF2   = dF2     * I_scale
psirz  = psirz  * I_scale

F_axis = f(1)

do i=1, nr
  factor  = sqrt(1.d0 + F_axis**2/F(i)**2 * (1.d0/I_scale**2 - 1.d0))
  dfactor = -1.d0/(factor) * (1.d0/I_scale**2 - 1.d0) * F_axis**2/F(i)**4 * dF2(i)
  q(i)   = factor * q(i)
  F(i)   = factor * I_scale * F(i)
enddo

if ((R_scale .ne. 1.d0) .or. (B_scale .ne. 1.d0) .or. (I_scale .ne. 1.d0)) then
  write(*,'(A)')               '********************************************************'
  write(*,'(A)')               '  Equilibrium scaling : '
  write(*,'(A,f8.4,A,f8.4,A)') '    R_scale : ',R_scale,'    major radius   :',rzero,' [m]'
  write(*,'(A,f8.4,A,f8.4,A)') '    B_scale : ',B_scale,'    vacuum field   :',bcentr,' [T]'
  write(*,'(A,f8.4,A,f8.4,A)') '    I_scale : ',I_scale,'    plasma current :',xip/1d6,' [MA]'
  write(*,'(A)')               '********************************************************'
endif

allocate(xx(nr),yy(nz))
do i=1,nr
  xx(i) = rgrid1 + xdim*real(i-1)/real(nr-1)
enddo     
do i=1,nz
  yy(i) = zmid + zdim*(real(i-1)/real(nz-1)-0.5)
enddo     

open(21,file='psiRZ_eqdsk')
do i=1,nr
  do j=1,nz
    write(21,'(3e14.6)') xx(i),yy(j),psirz(i,j)
  enddo
  write(21,'(A)') ''
enddo
close(21)


if (tokamak_name == 'ITER') then

  !--------------------close fit to ITER wall
  ellip  = 2.0
  tria_u = 0.55
  tria_l = 0.65
  quad_u = -0.1
  quad_l = 0.15
  n_tht  = 257
  r0     = 6.2  * R_scale
  z0     = 0.1  * R_scale
  a0     = 2.25 * R_scale

  !-------------------- contour outside ITER wall
  ellip  = 2.1
  tria_u = 0.58
  tria_l = 0.65
  quad_u = -0.12
  quad_l = -0.
  n_tht  = 257
  r0     = 6.2   * R_scale
  z0     = -0.05 * R_scale
  a0     = 2.34  * R_scale

else if (tokamak_name == 'JET') then
  
  !-------------------- contour outside JET wall
  ! blue contour in https://www.jorek.eu/wiki/doku.php?id=eqdsk2jorek.f90
  ellip  = 1.85
  tria_u = 0.4
  tria_l = 0.4
  quad_u = -0.2
  quad_l = -0.2
  n_tht  = 257
  r0     = 2.9  * R_scale
  z0     = 0.1  * R_scale
  a0     = 1.08 * R_scale
else if (tokamak_name == 'DIII-D') then

  !-------------------- contour outside DIII-D wall
  ellip  = 1.85
  tria_u = 0.4
  tria_l = 0.4
  quad_u = -0.2
  quad_l = -0.2
  n_tht  = 257
  r0     = 1.7 * R_scale
  z0     = 0.  * R_scale
  a0     = 0.7 * R_scale
  
  !-------------------- Atomic physics JOREK/NIMROD/M3D-C1 benchmark case (paper by B. Lyons)
  ellip  = 1.35/0.7
  tria_u = 0.3
  tria_l = 0.3
  quad_u = 0.
  quad_l = 0.
  n_tht   = 257
  r0     = 1.7 * R_scale
  z0     = 0.  * R_scale
  a0     = 0.7 * R_scale

else if (tokamak_name == 'CFETR') then

  ellip  = 2.43
  tria_u = 0.18
  tria_l = 0.18
  quad_u = -0.55
  quad_l = -0.55
  n_tht  = 257
  r0     = 7.12 * R_scale
  z0     = 0.15  * R_scale
  a0     = 2.43 * R_scale

else if (tokamak_name == 'EAST') then

  ellip  = 1.96
  tria_u = 0.55
  tria_l = 0.55
  quad_u = -0.0
  quad_l = -0.0
  n_tht  = 257
  r0     = 1.85 * R_scale
  z0     = 0.00  * R_scale
  a0     = 0.63 * R_scale


else if (tokamak_name == 'HL-2M') then

  ! ellip  = 2.10
  ! tria_u = 0.38
  ! tria_l = 0.68
  ! quad_u = -0.0
  ! quad_l = 0.2
  ! n_tht  = 257
  ! r0     = 1.80 * R_scale
  ! z0     = -0.2  * R_scale
  ! a0     = 0.65 * R_scale
  
else if (tokamak_name == 'EXL-50U') then

  ellip  = 2.8
  tria_u = 0.1
  tria_l = 0.1
  quad_u = -0.85
  quad_l = -0.85
  n_tht  = 257
  r0     =  0.85251* R_scale !0.8 !0.9
  z0     = -0.0  * R_scale
  a0     = 0.61 * R_scale !0.54 !0.66
else if (tokamak_name == 'EHL2') then

  ! ellip  = 2.1
  ! tria_u = 0.56
  ! tria_l = 0.56
  ! quad_u = 0.0
  ! quad_l = 0.0
  ! n_tht  = 257
  ! r0     =  1.17* R_scale
  ! z0     = -0.0  * R_scale
  ! a0     = 0.68 * R_scale
  ellip  = 2.4
  tria_u = 0.4
  tria_l = 0.4
  quad_u = 0.0 !0.22
  quad_l = 0.0 !0.22
  n_tht  = 257
  r0     =  1.1* R_scale
  z0     = -0.0  * R_scale
  a0     = 0.64 * R_scale

else

  write(*,*) 'Tokamak name not or wrongly specified, stopping'
  stop

end if  
  
PI = 2.d0 * asin(1.d0)

allocate(r_bnd(n_tht),z_bnd(n_tht),psi_bnd(n_tht))

!--------------------------------- interpolate flux using Dierckx spline routine
iopt= 0 
mx = nr
my = nz
xb = xx(1)
xe = xx(nr)
yb = yy(1)
ye = yy(nz)
kx = 3
ky = 3
smth = 1.d-6 ! Controls the tradeoff between closeness of fit and smoothness of fit. When too small, can lead to noise pick-up. When too large, can lead to inaccurate fit.
             ! May need hand tuning, based on a visual inspection of the output.
             ! For more details, see the documentation of regrid.f in libdierckx or the "Hard-coded parameters" section of the Wiki page https://www.jorek.eu/wiki/doku.php?id=eqdsk2jorek.f90. 
nxest = 3*nr/2!3*nr/4 ! Upper bound for the number of knots used for the splines. We set it a bit smaller than nr to test the quality of the fit.
nyest = 3*nz/2!3*nz/4
lwrk  = 4+nxest*(my+2*kx+5)+nyest*(2*ky+5)+mx*(kx+1)+my*(ky+1)+my+nxest
kwrk  = 3+mx+my+nxest+nyest

allocate(tx(nxest),ty(nyest),c(nxest,nyest),wrk(lwrk),iwrk(kwrk))

call regrid(iopt,mx,xx,my,yy,transpose(psirz),xb,xe,yb,ye,kx,ky,smth,nxest,nyest,nx,tx,ny,ty,c,fp,wrk,lwrk,iwrk,kwrk,ier)

if (ier > 0) then
  write(*,*) '!!!!! WARNING: Problem with the Dierckx spline interpolation !!!!!'
  write(*,*) '!!!!! You may need to tune smth and/or nxest and nyest.      !!!!!'
end if 
write(*,*) ' Dierckx ier   : ',ier
write(*,*) ' Dierckx fp    : ',fp
write(*,*) ' Dierckx nx,ny : ',nx,ny
if (ier > 0) then
  write(*,*) '!!!!! Exiting                                                !!!!!'
  stop
end if 

lwrk = mx*(kx+1)+my*(ky+1)
kwrk = mx+my
deallocate(wrk,iwrk)
allocate(wrk(lwrk),iwrk(kwrk))

!do i=1,nr
!  call bispev(tx,nx,ty,ny,c,kx,ky,xx(i),1,yy(nz/2),1,fout,wrk,lwrk,iwrk,kwrk,ier)
!  write(*,'(4e16.8,i3)') xx(i),yy(nz/2),fout,psirz(i,nz/2),ier
!enddo

do i=1,n_tht/2
  angle = 2.d0 * PI * float(i-1)/float(n_tht-1)
  r_bnd(i) = r0 + a0 * cos(angle + tria_u*sin(angle) + quad_u*sin(2.d0*angle))
  z_bnd(i) = z0 + a0 * ellip * sin(angle)
  call bispev(tx,nx,ty,ny,c,kx,ky,r_bnd(i),1,z_bnd(i),1,psi_bnd(i),wrk,lwrk,iwrk,kwrk,ier)
enddo
do i=n_tht/2+1,n_tht
  angle = 2.d0 * PI * float(i-1)/float(n_tht-1)
  r_bnd(i) = r0 + a0 * cos(angle + tria_l*sin(angle) + quad_l*sin(2.d0*angle))
  z_bnd(i) = z0 + a0 * ellip * sin(angle)
  call bispev(tx,nx,ty,ny,c,kx,ky,r_bnd(i),1,z_bnd(i),1,psi_bnd(i),wrk,lwrk,iwrk,kwrk,ier)
enddo

! write r_bnd(j), z_bnd(j)                     
open(21,file='JOREK_bnd.dat')
  do j=1,n_tht
    write(21,'(2e16.8)') r_bnd(j), z_bnd(j)
  enddo
close(21)


write(*,*) ' plotting results'  
nc = 51
allocate(zc(nc))
call begplt('eqdsk.ps')
call lblbot('eqdsk data',10)

call cplot(22,1,0,xx,yy,nr,nz,1,1,psirz,n_psi,zc,-nc,'fluxcontours',12,'R [m]',5,'Z [m]',5)
call lincol(3)
call lplot6(2,1,rlim,zlim,-limitr,'limiter')
call lincol(1)
call lplot6(2,1,rbnd,zbnd,-nbbs,'boundary')
call lincol(2)
call lplot6(2,1,r_bnd,z_bnd,-n_tht,'JOREK boundary')
call lincol(0)
call lplot6(3,2,psi,p,n_psi,'pressure')
call lplot6(3,3,psi,q,n_psi,'q')

call lplot6(2,2,psi,df2,n_psi,'df2')
call lplot6(3,2,psi,p,n_psi,'pressure')
call lplot6(2,3,psi,f,n_psi,'f')
call lplot6(3,3,psi,q,n_psi,'q')


!---------------------------- write JOREK input files
n_sol = (n_psi-1)/2
n_ext = n_psi + n_sol

write(*,*) ' n_psi, n_sol, n_ext : ',n_psi, n_sol, n_ext

allocate(df2_ext(n_ext),rho_ext(n_ext),T_ext(n_ext),T_ext_i(n_ext),T_ext_e(n_ext),psi_ext(n_ext),p_ext(n_ext))

df2_ext(1:n_psi) = df2(1:n_psi)
rho_ext(1:n_psi) = 1.d0

df2_ext(n_psi-1:n_ext) = df2_ext(n_psi)
rho_ext(n_psi-1:n_ext) = rho_ext(n_psi)

psi_sep = 1.0d0     ! in normalised psi units
sig_sep = 0.005     ! in normalised psi units
rho_bnd = 1d-4      ! in jorek units
T_bnd   = 1.d-6     ! in jorek units

!------------------------read the ne_profile-------------------------!
if (allocated(ne_spline)) then
  deallocate(ne_spline)
end if
allocate(ne_spline(n_ext),stat=err_alloc)
if (err_alloc /= 0) then
  write(6,*) "Error when trying to dynamically allocate memeries for ne_spline."
else
  inquire(file="./ne_profile_interp.dat",exist=ferr)
  if (ferr) then
    open(42,file="./ne_profile_interp.dat",status="OLD",action="READ")
    read(42,*) ne_spline(1:n_ext)
    close(42)
  else
    write(6,*) "Warning!!! ne_spline file does not found!"
    deallocate(ne_spline)
  end if
end if
!--------------------End of reading the ne_profile---------------------!

psi_ext(1:n_psi) = psi(1:n_psi)
do i=n_psi+1,n_ext
  psi_ext(i) = 1.d0 + 0.4 * float(i-n_psi)/float(n_sol)
enddo

! Taper FF' to zero in the SOL (no toroidal current outside LCFS)
do i = n_psi+1, n_ext
  df2_ext(i) = df2_ext(n_psi) * exp(-100.d0 * (psi_ext(i) - 1.d0))
end do

zmu0 = 4.d-7 * PI
! Extend pressure into SOL using monotonic cubic Hermite spline
! (Steffen's method) to guarantee C1 continuity at psi=1 and
! suppress oscillations. sig_sol controls the SOL decay width.

sig_sol  = 0.03d0           ! SOL pressure decay width in psi_n
pres_sol = p(n_psi)*0.01d0  ! SOL pressure floor (1% of separatrix)

call extend_pressure_cubic(n_psi, psi, p, n_ext, psi_ext, p_ext, sig_sol, pres_sol)
do i=1,n_ext
  tanh1 = tanh((psi_ext(i) - psi_sep)/sig_sep)
  !rho_ext(i) = (rho_ext(i) - rho_bnd) * (0.5d0 - 0.5d0*tanh1) + rho_bnd
  
  if (allocated(ne_spline)) then
    rho_ext(i) = rho_ext(i) * ne_spline(i)!+rho_bnd * (0.5 + 0.5*tanh1)
   T_ext(i)   = p_ext(i)*zmu0 / rho_ext(i) +T_bnd * (0.5 + 0.5*tanh1)
!   print *,'psi=',psi_ext(i)
!   print *,'tanh1=',tanh1
  else
    rho_ext(i) = (rho_ext(i) - rho_bnd) * (0.5d0 - 0.5d0*tanh1) + rho_bnd
    T_ext(i)   = p_ext(i)*zmu0 / rho_ext(i)!+T_bnd * (0.5 + 0.5*tanh1)
  end if
enddo
!p_ext(1:n_psi) = T_ext(1:n_psi)*rho_ext(1:n_psi)/zmu0
coef_Ti_Te= 2.0
T_ext_e=T_ext/(1+coef_Ti_Te)
T_ext_i=T_ext-T_ext_e

call lplot6(2,2,psi_ext,df2_ext,n_ext,'df2')
call lplot6(3,2,psi_ext,p_ext,n_ext,'pressure')
call lplot6(2,3,psi_ext,rho_ext,n_ext,'density')
call lplot6(3,3,psi_ext,T_ext,n_ext,'T')

call lincol(1)
call lplot6(2,2,psi,df2,-n_psi,'df2')
call lplot6(3,2,psi,p*zmu0,-n_psi,'pressure')

open(21,file='jorek_ffprime')
! We change or not the sign of ff' depending on the sign of Ip because (we assume that) in EQDSK files, 
! psi_axis is always < psi_boundary, whatever the direction of Ip.
if (xip>0) then
  do i=1,n_ext
    write(21,*) psi_ext(i), df2_ext(i) ! The minus sign is because ff' in JOREK is opposite to the usual ff' for historical reasons.
  enddo  
else
  do i=1,n_ext
    write(21,*) psi_ext(i),-df2_ext(i) 
  enddo  
end if
close(21)

open(21,file='jorek_density')
do i=1,n_ext
  write(21,*) psi_ext(i),rho_ext(i)
enddo
close(21)
open(21,file='jorek_temperature')
do i=1,n_ext
  write(21,*) psi_ext(i),T_ext(i)
enddo
close(21)
open(21,file='jorek_temperature_i')
do i=1,n_ext
  write(21,*) psi_ext(i),T_ext_i(i)
enddo
close(21)
open(21,file='jorek_temperature_e')
do i=1,n_ext
  write(21,*) psi_ext(i),T_ext_e(i)
enddo
close(21)



open(21,file='jorek_namelist')


write(21,*)             '***************************************'
write(21,'(A)')        '*  namelist produced by eqdsk2jorek   *'
write(21,*)             '***************************************'
write(21,'(A)')        '   TEQ G      05/14  /2019     DThmode24    500ms   3'
write(21,'(A,f8.3,A)') '   magnetic field   : ',Bcentr,' [T]'
write(21,'(A,f8.3,A)') '   current          : ',xip/1d6,' [MA]'
write(21,'(A,e14.6,A)')'   central pressure : ',p(1), '[Pa]'
write(21,*)             '***************************************'
write(21,*)
write(21,*) ' &in1'
write(21,*) ' restart = .f.'
write(21,*) ' regrid  = .f.'
write(21,*) ' equil   = .t.'
write(21,*) ' regrid_from_rz=.f.'
write(21,*) ' tstep_n   = 1!0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10'
write(21,*) ' nstep_n   = 0!100,200,100,200,100,100,200,100,100,1000'
write(21,*) ' nout = 1!100'
write(21,*)
write(21,*) ' time_evol_scheme= "Gears"'
write(21,*)  
write(21,*) ' tgnum_u   = 0.5'
write(21,*) ' tgnum_zj  = 0.5'
write(21,*) ' tgnum_w   = 0.5'
write(21,*) ' tgnum_rho = 0.5'
write(21,*) ' tgnum_T   = 0.5'
write(21,*) ' tgnum_vpar= 0.5'
write(21,*) ' tgnum_rhon= 0.5'
write(21,*) ' tgnum_nre = 0.5'
write(21,*) ' tgnum_AR  = 0.5'
write(21,*) ' tgnum_AZ  = 0.5'
write(21,*) ' tgnum_A3  = 0.5'
write(21,*)
write(21,*) ' gmres_4     = 1.d4'
write(21,*) ' iter_precon = 21'
write(21,*) ' gmres_tol   = 1.d-8'
write(21,*) ' gmres_m     = 20'
write(21,*)
write(21,*) ' !use_mumps =.true.'
write(21,*) ' use_pastix = .true.'
write(21,*) ' use_strumpack=.false.'
write(21,*) ' use_strumpack_eq=.false.'
write(21,*) ' write_ps = .false.'
write(21,*)
write(21,*) ' !_____________________________________limiter definition'
write(21,*) ' mf = 0'
!write(21,*) ' n_limiter = ',limitr
write(21,*) ' n_boundary = ',n_tht
write(21,*)
do j=1,n_tht
  ! write(21,'(A,i3,A,e16.8,A,i3,A,e16.8,A)'), &
  !          '  R_limiter(',j,') =',rlim(j), &
  !          ', Z_limiter(',j,') =',zlim(j), &
  !          ','
  write(21,'(A,i3,A,e16.8,A,i3,A,e16.8,A,i3,A,e16.8,A)'), &
           '  R_boundary(',j,') =',r_bnd(j), &
           ', Z_boundary(',j,') =',z_bnd(j), &
           ', psi_boundary(',j,') =',psi_bnd(j), &
           ','
enddo
write(21,*)
write(21,*) ' ellip  = ',ellip
write(21,*) ' tria_u = ',tria_u
write(21,*) ' tria_l = ',tria_l
write(21,*) ' quad_u = ',quad_u
write(21,*) ' quad_l = ',quad_l
write(21,*)
write(21,*) ' xampl  = +0.'
write(21,*) ' xpoint = .t.'
write(21,*)
write(21,*) ' freeboundary = .f.'
write(21,*) ' resistive_wall = .f.'
write(21,*)
write(21,*) ' R_geo = ',r0
write(21,*) ' Z_geo = ',z0
if (tokamak_name=='JET') then
  write(21,*) ' F0    = ',-2.96*bcentr ! By convention, the vacuum toroidal field is given at 2.96m in JET eqdsk files. 					
else
  write(21,*) ' F0    = ',-rzero*bcentr
end if
write(21,*) ' amin  = 1.d0 ! scale factor for plasma size only'
write(21,*)
write(21,*) ' psi_axis_init  = ', ssimag 
write(21,*) ' RZ_grid_inside_wall = .t.'
write(21,*) ' grid_to_wall        = .t.'
write(21,*)
write(21,*) ' !_____________________________________grid parameters'
write(21,*) ' n_R      = 200'
write(21,*) ' n_Z      = 600'
write(21,*) ' n_radial = 0'
write(21,*) ' n_pol    = 0'
write(21,*) ' n_flux   = 120'
write(21,*) ' n_tht    = 160'
write(21,*) ' n_open   = 15'
write(21,*) ' n_leg    = 20'
write(21,*) ' n_private = 2'
write(21,*) ' SIG_closed = 0.2'
write(21,*) ' SIG_theta = 0.1'
write(21,*) ' SIG_open  = 0.3'
write(21,*) ' SIG_private= 0.3'
write(21,*) ' dPSI_open    = 0.023'
write(21,*) ' dPSI_private = 0.001'
write(21,*) ' amix = 0.1d0'
write(21,*) ' !_____________________extend nodes_____________'
write(21,*) '   n_wall_blocks =0'
write(21,*) ' !----First block'
write(21,*) '   n_ext_block         ( 1)    = 8'
write(21,*) ' n_block_points_left ( 1)    =  2'
write(21,*) ' R_block_points_left ( 1,1)  = +4.21506'
write(21,*) ' Z_block_points_left ( 1,1)  = -3.80080'
write(21,*) ' R_block_points_left ( 1,2)  = 4.1794'
write(21,*) ' Z_block_points_left ( 1,2)  = -3.8725'
write(21,*) ' n_block_points_right( 1)    = 3'
write(21,*) ' R_block_points_right( 1,1)  = 4.69925'
write(21,*) ' Z_block_points_right( 1,1)  = -3.63420'
write(21,*) ' R_block_points_right( 1,2)  = 4.6453'
write(21,*) ' Z_block_points_right( 1,2)  = -3.7332'
write(21,*) ' R_block_points_right( 1,3)  = 4.4951'
write(21,*) ' Z_block_points_right( 1,3)  = -3.8965'
write(21,*)
write(21,*) ' !_____________________________________physics parameters'
write(21,*) ' eta   = 1.4957E-9'
write(21,*) ' eta_ohmic= 1.4957E-9'
write(21,*) ' visco_par = 1.d-5'
write(21,*) ' visco     = 5.d-6'
write(21,*) ' central_density = 0.973'
write(21,*) ' central_mass = 2'
write(21,*) ' visco_T_dependent = .f.'
write(21,*) ' eta_T_dependent = .t.'
write(21,*) ' eta_num_T_dependent=.t.'
write(21,*) ' visco_num_T_dependent  = .f.'
write(21,*)
write(21,*) ' eta_num       = 1.d-15'
write(21,*) ' visco_num     = 1.d-14'
write(21,*) ' visco_par_num = 1.d-14'
write(21,*) ' d_perp_num    = 1.d-14'
write(21,*) ' zk_perp_num   = 1.d-14'
write(21,*) ' bc_natural_open = .true.'
write(21,*) ' gamma_stangeby = 8'
write(21,*) ' bcs(9)%dirichlet%T    = .t.'
write(21,*) ' bcs(9)%dirichlet%rho  = .t.'
write(21,*) ' bcs(9)%dirichlet%vpar = .t.'
write(21,*) ' bcs(9)%natural%T      = .f.'
write(21,*) ' bcs(9)%natural%vpar   = .f.'
write(21,*) ' bcs(9)%natural%rho    = .f.'
write(21,*) ' bcs(9)%mach1          = .f.'
write(21,*) ' bcs(11)%dirichlet%T   = .f.'
write(21,*) ' bcs(11)%dirichlet%rho = .f.'
write(21,*) ' bcs(11)%dirichlet%vpar= .f.'
write(21,*) ' bcs(11)%natural%rho   = .t.'
write(21,*) ' bcs(11)%natural%vpar  = .t.'
write(21,*) ' bcs(11)%natural%T     = .t.'
write(21,*) ' bcs(11)%natural%rhon  = .t.'
write(21,*) ' bcs(11)%mach1         = .t.'
write(21,*) ' bcs(5)%dirichlet%T    = .t.'
write(21,*) ' bcs(5)%dirichlet%rho  = .t.'
write(21,*) ' bcs(5)%dirichlet%vpar = .t.'
write(21,*) ' bcs(5)%natural%rho    = .f.'
write(21,*) ' bcs(5)%natural%vpar   = .f.'
write(21,*) ' bcs(5)%natural%T      = .f.'
write(21,*) ' bcs(5)%mach1          = .f.'
write(21,*) ' bcs(2)%dirichlet%rho  = .t.'
write(21,*) ' bcs(2)%dirichlet%vpar = .t.'
write(21,*) ' bcs(2)%dirichlet%T    = .t.'
write(21,*) ' bcs(2)%natural%T      = .f.'
write(21,*) ' bcs(2)%natural%rhon   = .f.'
write(21,*) ' bcs(2)%natural%rho    = .f.'
write(21,*) ' bcs(2)%natural%vpar   = .f.'
write(21,*) ' bcs(2)%mach1          = .f.'
write(21,*) ' bcs(15)%dirichlet%T   = .t.'
write(21,*) ' bcs(15)%dirichlet%rho = .t.'
write(21,*) ' bcs(15)%dirichlet%vpar= .t.'
write(21,*) ' bcs(15)%natural%rho   = .f.'
write(21,*) ' bcs(15)%natural%vpar  = .f.'
write(21,*) ' bcs(15)%natural%T     = .f.'
write(21,*) ' bcs(15)%mach1         = .f.'
write(21,*) ' use_simple_bnd_types  = .t.'
write(21,*) ' rho_file     = "jorek_density"'
write(21,*) ' T_file       = "jorek_temperature"'
write(21,*) ' ffprime_file = "jorek_ffprime"'
write(21,*)
write(21,*) ' D_par     = 0.d0'
write(21,*) ' D_perp(1) = 1.d-6'
write(21,*) ' D_perp(2) = 0.9d0'
write(21,*) ' D_perp(3) = 0.1d0'
write(21,*) ' D_perp(4) = 0.01d0'
write(21,*) ' D_perp(5) = 0.9d0'
write(21,*) ' D_perp(6) = 1.d-6'
write(21,*)
write(21,*) ' ZK_par     = 3.2412E+06'
write(21,*) ' ZK_perp(1) = 1.d-6'
write(21,*) ' ZK_perp(2) = 0.9d0'
write(21,*) ' ZK_perp(3) = 0.1d0'
write(21,*) ' ZK_perp(4) = 0.01d0'
write(21,*) ' ZK_perp(5) = 0.9d0'
write(21,*) ' ZK_perp(6) = 1.d-6'
write(21,*)
write(21,*) ' tauIC      = -5.2682E-04'
write(21,*) ' T_min= 1.d-6'
write(21,*) ' T_min_neg= 1.d-6'
write(21,*) ' rho_min=1d-5'
write(21,*) ' rho_min_neg=1d-5'
write(21,*) ' heatsource     = 1.d-7'
write(21,*) ' particlesource = 5.d-6'
write(21,*) ' eqdsk_psi_fact = 1.0'
write(21,*) ' RZ_grid_jump_thres = 0.9'
write(21,*)
write(21,*) ' &end'
close(21)

open(22,file='jorek_limit')
do j=1,limitr
  write(22,'(A,i3,A,e16.8,A,i3,A,e16.8,A)'), &
           '  R_limiter(',j,') =',rlim(j), &
           ', Z_limiter(',j,') =',zlim(j), &
           ','
enddo
close(22)
open(23,file='CfetrPSI.dat')
do i=1,nz 
  do j=1,nr
    write(23,'(e16.8)'), &
           psirz(i,j)
          
  enddo
enddo
close(23)

call finplt

!-------------------------------------- eqdsk to vtk (careful VTK expects single precision)
lf   = char(10)
ivtk = 23

write(*,'(A)') ' writing VTK output'

open(unit=ivtk,file='eqdsk.vtk',form='binary',convert='BIG_ENDIAN')

buffer = '# vtk DataFile Version 3.0'//lf                        ; write(ivtk) trim(buffer)
buffer = 'eqdsk'//lf                                             ; write(ivtk) trim(buffer)
buffer = 'BINARY'//lf                                            ; write(ivtk) trim(buffer)
buffer = 'DATASET RECTILINEAR_GRID'//lf                          ; write(ivtk) trim(buffer)

write(str2(1:24),'(3i8)') nr,nz,1
buffer = 'DIMENSIONS '//str2//lf                             ; write(ivtk) trim(buffer)

write(str1(1:12),'(i12)') nr
buffer = 'X_COORDINATES '//str1//' FLOAT'//lf                    ; write(ivtk) trim(buffer)
write(ivtk) (real(xx(i),4),i=1,nr)

write(str1(1:12),'(i12)') nz
buffer = 'Y_COORDINATES '//str1//' FLOAT'//lf                    ; write(ivtk) trim(buffer)
write(ivtk) (real(yy(j),4),j=1,nz)

write(str1(1:12),'(i12)') 1
buffer = 'Z_COORDINATES '//str1//' FLOAT'//lf                    ; write(ivtk) trim(buffer)
write(ivtk) real(0.d0,4)

! POINT_DATA SECTION
write(str1(1:12),'(i12)') nr*nz
buffer = lf//lf//'POINT_DATA '//str1//lf                             ; write(ivtk) trim(buffer)

buffer = 'SCALARS psi float'//lf                                     ; write(ivtk) trim(buffer)
buffer = 'LOOKUP_TABLE default'//lf                                  ; write(ivtk) trim(buffer)
write(ivtk) ((real(psirz(i,j),4), i=1,nr), j=1,nz)

close(ivtk)

end



!=======================================================================
! Extend pressure into SOL using monotonic cubic Hermite interpolation.
! Guarantees C1 continuity at psi=1 and suppresses oscillations via
! Steffen's monotonic derivative formula.
!=======================================================================
subroutine extend_pressure_cubic(n_psi, psi, p, n_ext, psi_ext, p_ext, sig_sol, pres_sol)
  implicit none

  integer, intent(in)  :: n_psi, n_ext
  real*8,  intent(in)  :: psi(n_psi), p(n_psi)
  real*8,  intent(in)  :: psi_ext(n_ext)
  real*8,  intent(in)  :: sig_sol, pres_sol
  real*8,  intent(out) :: p_ext(n_ext)

  integer, parameter :: n_sol_ctrl = 3
  integer :: n_ctrl, i
  real*8, allocatable :: x(:), y(:), d(:)
  real*8 :: p_sep

  p_sep = p(n_psi)

  ! Build control points: interior (last few) + SOL extension
  ! Use last 2 interior points to capture the edge derivative
  n_ctrl = 2 + n_sol_ctrl
  allocate(x(n_ctrl), y(n_ctrl), d(n_ctrl))

  x(1) = psi(n_psi-1)
  x(2) = psi(n_psi)          ! = 1.0
  x(3) = psi(n_psi) + sig_sol
  x(4) = psi(n_psi) + 3.0d0 * sig_sol
  x(5) = psi(n_psi) + 6.0d0 * sig_sol

  y(1) = p(n_psi-1)
  y(2) = p_sep
  y(3) = p_sep * 0.1d0
  y(4) = max(pres_sol, p_sep * 0.01d0)
  y(5) = max(pres_sol, p_sep * 0.001d0)

  ! Compute monotonic derivatives (Steffen's method)
  call steffen_deriv(n_ctrl, x, y, d)

  ! Copy interior points directly
  p_ext(1:n_psi) = p(1:n_psi)

  ! Evaluate cubic Hermite on SOL output grid
  do i = n_psi+1, n_ext
    call hermite_eval(n_ctrl, x, y, d, psi_ext(i), p_ext(i))
  end do

  deallocate(x, y, d)
end subroutine extend_pressure_cubic



!=======================================================================
! Steffen's monotonic derivative formula.
! Guarantees the cubic Hermite interpolant is monotone on each interval,
! preventing oscillations that plague standard cubic splines.
! Ref: Steffen, M., "A simple method for monotonic interpolation",
!      Astronomy & Astrophysics, 239, 443 (1990).
!=======================================================================
subroutine steffen_deriv(n, x, y, d)
  implicit none

  integer, intent(in)  :: n
  real*8,  intent(in)  :: x(n), y(n)
  real*8,  intent(out) :: d(n)

  integer :: i
  real*8  :: h(n-1), s(n-1), p, q

  ! Slopes of each segment
  do i = 1, n-1
    h(i) = x(i+1) - x(i)
    s(i) = (y(i+1) - y(i)) / h(i)
  end do

  ! Endpoint derivatives: one-sided quadratic
  d(1) = ((2.d0*h(1)+h(2))*s(1) - h(1)*s(2)) / (h(1)+h(2))
  d(n) = ((2.d0*h(n-1)+h(n-2))*s(n-1) - h(n-1)*s(n-2)) / (h(n-1)+h(n-2))

  ! Monotonicity-preserving correction on d(1)
  if (d(1)*s(1) <= 0.d0) then
    d(1) = 0.d0
  else if (abs(d(1)) > 3.d0*abs(s(1))) then
    d(1) = 3.d0 * s(1)
  end if

  ! Monotonicity-preserving correction on d(n)
  if (d(n)*s(n-1) <= 0.d0) then
    d(n) = 0.d0
  else if (abs(d(n)) > 3.d0*abs(s(n-1))) then
    d(n) = 3.d0 * s(n-1)
  end if

  ! Interior points: Steffen's formula with harmonic mean
  do i = 2, n-1
    p = (s(i-1)*h(i) + s(i)*h(i-1)) / (h(i-1) + h(i))
    if (s(i-1)*s(i) <= 0.d0) then
      d(i) = 0.d0           ! local extremum → zero derivative
    else
      q = 2.d0 * s(i-1) * s(i) / (s(i-1) + s(i))
      ! Steffen: choose between p (Bessel) and q (harmonic mean)
      ! with a smoothness weight
      d(i) = (sign(1.d0, p) * min(abs(p), abs(q)))
    end if
  end do

end subroutine steffen_deriv



!=======================================================================
! Evaluate cubic Hermite piecewise polynomial at a given x0.
! Uses the standard Hermite basis functions:
!   h00(t) = 2t^3 - 3t^2 + 1
!   h10(t) = t^3 - 2t^2 + t
!   h01(t) = -2t^3 + 3t^2
!   h11(t) = t^3 - t^2
!=======================================================================
subroutine hermite_eval(n, x, y, d, x0, y0)
  implicit none

  integer, intent(in)  :: n
  real*8,  intent(in)  :: x(n), y(n), d(n), x0
  real*8,  intent(out) :: y0

  integer :: left, right, mid
  real*8  :: t, h, h00, h10, h01, h11

  ! Binary search to find the interval containing x0
  left  = 1
  right = n
  do
    if (right <= left + 1) exit
    mid = (left + right) / 2
    if (x(mid) >= x0) then
      right = mid
    else
      left = mid
    end if
  end do

  ! Cubic Hermite evaluation on [x(left), x(right)]
  h  = x(right) - x(left)
  t  = (x0 - x(left)) / h

  h00 =  2.d0*t**3 - 3.d0*t**2 + 1.d0
  h10 = (t**3 - 2.d0*t**2 + t) * h
  h01 = -2.d0*t**3 + 3.d0*t**2
  h11 = (t**3 - t**2) * h

  y0 = h00*y(left) + h10*d(left) + h01*y(right) + h11*d(right)

end subroutine hermite_eval
 
 
