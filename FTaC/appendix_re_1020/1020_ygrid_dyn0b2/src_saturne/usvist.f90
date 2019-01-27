!===============================================================================

!===============================================================================
! Purpose:
! -------

!> usvist
!> \brief Modify turbulent viscosity
!>
!> This subroutine is called at beginning of each time step
!> after the computation of the turbulent viscosity
!> (physical quantities have already been computed in \ref usphyv).
!>
!> Turbulent viscosity \f$ \mu_T \f$ (kg/(m s)) can be modified.
!>
!> A modification of the turbulent viscosity can lead to very
!> significant differences betwwen solutions and even give wrong
!> results.
!>
!> This subroutine is therefore reserved to expert users.
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss
!> \param[in]     ncesmp        number of cells with mass source term
!> \param[in]     icepdc        head loss cell numbering
!> \param[in]     icetsm        numbering of cells with mass source term
!> \param[in]     itypsm        kind of mass source for each variable
!>                               (cf. \ref cs_user_mass_source_terms)
!> \param[in]     dt            time step (per cell)
!> \param[in]     ckupdc        work array for head loss terms
!> \param[in]     smacel        values of variables related to mass source
!>                              term. If ivar=ipr, smacel=mass flux
!_______________________________________________________________________________

subroutine usvist &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   ckupdc , smacel )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use numvar
use optcal
use cstphy
use entsor
use parall
use period
use mesh
use field
use field_operator

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)

! Local variables

integer          iel, inc
integer          iprev

double precision coef, deux, delta
double precision s11, s22, s33
double precision dudy, dudz, dvdx, dvdz, dwdx, dwdy
double precision xfil, xa  , xb  , radeux

double precision, dimension(:,:,:), allocatable :: gradv
double precision, dimension(:,:), pointer :: coefau
double precision, dimension(:,:,:), pointer :: coefbu
double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: visct

double precision mytmpdble, mycsmago2
! Variables with save attribute
integer, save :: ny, nxnz
integer, allocatable, dimension(:), save :: map_3d_1d
double precision, allocatable, dimension(:), save :: ygrid, moyenne
! Other variables
double precision :: xmin, zmin, yy!, eps
character(len=64) :: buffer
integer :: nyloc
integer, allocatable, dimension(:) :: lstelt, ncell_3d_1d
double precision, allocatable, dimension(:) :: ygridloc, w1d

!===============================================================================

if (ntcabs.le.100000) return ! nothing done before 100000 time steps

! Set up spatial averaging
if (.not.allocated(map_3d_1d)) then

  ! Compute minimal x and z values
  xmin = minval(xyzcen(1,1:ncel))
  zmin = minval(xyzcen(3,1:ncel))
!  eps = (minval(volume(1:ncel)))**(1.d0/3.d0)
  if (irangp.ge.0) then
    call parmin(xmin)
    call parmin(zmin)
!    call parmin(eps)
  endif

  ! Compute number of y cells
  ! Assume domain start at x=0, z=0 and has x>0 and z>0
  write (buffer, "(A4,ES12.5,A9,ES12.5)") "x < ",1.5*xmin," and  z < ",1.5*zmin
  allocate(lstelt(ncel))
  call getcel(buffer, nyloc, lstelt)
  ny = nyloc
  if (irangp.ge.0) then
    call parcpt(ny)
  endif
  if (irangp.le.0) then
    write(nfecra,*) 'Number of cells in y-direction : ',ny
  endif
  if (ny.eq.0) then
    if (irangp.le.0) write(nfecra,*) 'Problem with 1d y-grid, ABORT simulation'
    if (irangp.le.0) write(nfecra,*) 'User must adapt subroutine ussmag'
    call csexit(12345678)
  endif

  ! Compute global y-grid
  allocate(ygridloc(nyloc))
  if (nyloc.ge.1) then
    do iel=1,nyloc
      ygridloc(iel) = xyzcen(2,lstelt(iel))
    enddo
  endif
  allocate(ygrid(ny))
  if (irangp.ge.0) then
    call paragv(nyloc,ny,ygridloc,ygrid)
  else
    ygrid = ygridloc
  endif

  ! Allocate 1d moyenne
  allocate(moyenne(ny))
  moyenne = 0.

  ! Compute local 3D-1D mapping
  allocate(map_3d_1d(ncel))
  do iel = 1,ncel
    yy = xyzcen(2,iel)
    map_3d_1d(iel) = minloc( abs(ygrid-yy), 1)
  enddo

  ! Check mesh is homogeneous
  allocate(ncell_3d_1d(ny))
  ncell_3d_1d = 0
  do iel = 1,ncel
    ncell_3d_1d( map_3d_1d(iel) ) = ncell_3d_1d( map_3d_1d(iel) ) + 1
  enddo
  if (irangp.ge.0) call parism(ny,ncell_3d_1d)
  nxnz = ncell_3d_1d(1)
  if (irangp.le.0) then
    write(nfecra,*) 'Number of cells in y-constant planes : ',nxnz
  endif
  if (maxval(ncell_3d_1d).ne.minval(ncell_3d_1d)) then
    if (irangp.le.0) write(nfecra,*) 'Problem with 1d y-grid, ABORT simulation'
    if (irangp.le.0) write(nfecra,*) 'User must adapt subroutine ussmag'
    call csexit(123456789)
  endif

endif

!===============================================================================
! 1.  Initialization
!===============================================================================

call field_get_coefa_v(ivarfl(iu), coefau)
call field_get_coefb_v(ivarfl(iu), coefbu)

! Allocate temporary arrays for gradients calculation
allocate(gradv(3, 3, ncelet))

call field_get_val_s(iprpfl(ivisct), visct)
call field_get_val_s(icrom, crom)

! --- For the calculation of viscosity on the sub-mesh
xfil   = xlesfl
xa     = ales
xb     = bles
deux   = 2.d0
radeux = sqrt(deux)

!===============================================================================
! 2.  Calculation of velocity gradient and of
!       S11**2+S22**2+S33**2+2*(S12**2+S13**2+S23**2)
!===============================================================================

inc = 1
iprev = 1

call field_gradient_vector(ivarfl(iu), iprev, imrgra, inc,    &
                           gradv)

do iel = 1, ncel

  s11  = gradv(1, 1, iel)
  s22  = gradv(2, 2, iel)
  s33  = gradv(3, 3, iel)
  dudy = gradv(2, 1, iel)
  dvdx = gradv(1, 2, iel)
  dudz = gradv(3, 1, iel)
  dwdx = gradv(1, 3, iel)
  dvdz = gradv(3, 2, iel)
  dwdy = gradv(2, 3, iel)

  mytmpdble = s11**2 + s22**2 + s33**2                &
                     + 0.5d0*((dudy+dvdx)**2          &
                     +        (dudz+dwdx)**2          &
                     +        (dvdz+dwdy)**2)

  delta  = xfil* (xa*volume(iel))**xb
!  coef = csmago**2 * radeux
!  delta  = coef * delta**2
!  visct(iel) = crom(iel) * delta * sqrt(mytmpdble)
  mycsmago2 = visct(iel) / (crom(iel) * sqrt(mytmpdble) * (delta**2) * radeux )
  moyenne( map_3d_1d(iel) ) = moyenne( map_3d_1d(iel) ) &
      + (1./dble(ntmabs-100000+1))*mycsmago2/dble(nxnz)
enddo

if (ntcabs.eq.ntmabs) then
  if (irangp.ge.0) then
    call parrsm(ny,moyenne)
  endif
  if (irangp.le.0) then
    do iel = 1,ny
      write(nfecra,*) 'CEDRIC MOYENNE CSMAGO : ',ygrid(iel),moyenne(iel)
    enddo
  endif
endif


! Free memory
deallocate(gradv)

return
end subroutine usvist
