!-------------------------------------------------------------------------------

!                      Code_Saturne version 5.0-alpha
!                      --------------------------
! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2017 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

!===============================================================================
! Purpose:
! -------

!> \file cs_user_extra_operations.f90
!>
!> \brief This function is called at the end of each time step, and has a very
!>  general purpose
!>  (i.e. anything that does not have another dedicated user subroutine)
!>
!> See \subpage cs_user_extra_operations_examples and
!> \subpage cs_user_extra_operations-nusselt_calculation for examples.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________

subroutine cs_f_user_extra_operations &
 ( nvar   , nscal  ,                                              &
   dt     )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use lagran
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh
use field
use field_operator
use turbomachinery
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables
integer          iel, ivar, ivar2, nceltot
integer          errout, impout
integer          f_id, ifcvsl, iflwgr
integer          npoint, ii, iel1, irang1, irangv

double precision x

double precision errl2, errl2p
double precision sumerrl2, sumerrl2p

double precision errl1, errl1p
double precision sumerrl1, sumerrl1p

double precision errli, errlip

double precision xyz(3), sc, scth, time_evol_m, time_evol_p, dtdx, dtth
double precision, dimension(:), pointer :: cvar_scal, cvar_exa, cvar_err
double precision, dimension(:), pointer :: cvar_scal2, cvar_err2
double precision, dimension(:), pointer :: cpro_wgrec_s
double precision, dimension(:,:), pointer :: grad, grad2, grade

type(var_cal_opt) vcopt

!===============================================================================
! 1.  Initialization
!===============================================================================
if (ntcabs.eq.ntmabs) then

  ! Error on field
  call field_get_id('exa', f_id)
  call field_get_val_s(f_id, cvar_exa)

  call field_get_id('error', f_id)
  call field_get_val_s(f_id, cvar_err)
  call field_get_id('error2', f_id)
  call field_get_val_s(f_id, cvar_err2)

  ivar = isca(1)
  call field_get_val_s(ivarfl(ivar), cvar_scal)

  ivar2 = isca(2)
  call field_get_val_s(ivarfl(ivar2), cvar_scal2)

  if(irangp.le.0) then
    errout = impusr(1)
    impout = impusr(2)
    open(errout,file='error.dat')
  endif

  errl2 = 0d0
  errl2p = 0d0

  errl1 = 0d0
  errl1p = 0d0

  errli = 0d0
  errlip = 0d0

  do iel = 1, ncel
    x = xyzcen(1,iel)

    errl2 = errl2 + volume(iel) * (cvar_scal(iel)-cvar_exa(iel))**2
    errl2p = errl2p + volume(iel) * cvar_exa(iel)**2

    errl1 = errl1 + volume(iel) * abs(cvar_scal(iel)-cvar_exa(iel))
    errl1p = errl1p + volume(iel) * abs(cvar_exa(iel))

    errli = max(errli, abs(cvar_scal(iel)-cvar_exa(iel)))
    errlip = max(errlip, abs(cvar_exa(iel)))

    cvar_err(iel) = abs(cvar_scal(iel)-cvar_exa(iel))
    cvar_err2(iel) = abs(cvar_scal2(iel)-cvar_exa(iel))
  enddo

  sumerrl2 = errl2
  sumerrl2p = errl2p

  sumerrl1 = errl1
  sumerrl1p = errl1p

  nceltot = ncel

  if (irangp.ge.0) then
    call parsom(sumerrl2)
    call parsom(sumerrl2p)
    call parsom(sumerrl1)
    call parsom(sumerrl1p)
    call parmax(errli)
    call parmax(errlip)
    call parcpt(nceltot)
  endif

  if(irangp.le.0) then
    errl2 = sqrt(sumerrl2/sumerrl2p)
    errl1 = sumerrl1/sumerrl1p
    errli = errli/errlip
    write(errout,*) '# ncel^(1/3) | erreur L^1 | erreur L^2 | erreur L^inf'
    write(errout,*) nceltot**(1d0/3d0), errl1 , errl2, errli
    close(errout)
  endif

  ! Error on gradient
  call field_get_id('err_grad', f_id)
  call field_get_val_s(f_id, cvar_err)

  call field_get_id('err_grad2', f_id)
  call field_get_val_s(f_id, cvar_err2)

  allocate(grad(3,ncelet))
  call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)
  if (vcopt%iwgrec.eq.1) then
    call field_get_key_int(ivarfl(ivar), kwgrec, iflwgr)
!    call field_get_val_s(iflwgr, cpro_wgrec_s)
!    call synsca(cpro_wgrec_s)
!    if (irangp.ge.0.or.iperio.eq.1) then
!      call synsca(cpro_wgrec_s)
!    endif
  endif
  call field_gradient_scalar &
    (ivarfl(ivar), 0, vcopt%imrgra, 1, 1, grad)

  allocate(grad2(3,ncelet))
  call field_get_key_struct_var_cal_opt(ivarfl(ivar2), vcopt)
  if (vcopt%iwgrec.eq.1) then
    call field_get_key_int(ivarfl(ivar2), kwgrec, iflwgr)
!    call field_get_val_s(iflwgr, cpro_wgrec_s)
!    call synsca(cpro_wgrec_s)
!    if (irangp.ge.0.or.iperio.eq.1) then
!      call synsca(cpro_wgrec_s)
!    endif
  endif
  call field_gradient_scalar &
    (ivarfl(ivar2), 0, vcopt%imrgra, 1, 1, grad2)

  allocate(grade(3,ncelet))
  call solgraddiff(grade)

  if(irangp.le.0) then
    errout = impusr(1)
    impout = impusr(2)
    open(errout,file='error_grad.dat')
  endif

  errl2 = 0d0
  errl2p = 0d0

  errl1 = 0d0
  errl1p = 0d0

  errli = 0.d0
  errlip = 0.d0

  do iel = 1, ncel
    errl2 = errl2 + volume(iel) * sum((grad(1:3,iel)-grade(1:3,iel))**2)
    errl2p = errl2p + volume(iel) * sum(grade(1:3,iel)**2)

    errl1 = errl1 + volume(iel) * maxval(abs(grad(1:3,iel)-grade(1:3,iel)))
    errl1p = errl1p + volume(iel) * maxval(abs(grade(1:3,iel)))

    errli = max(errli, maxval(abs(grad(1:3,iel)-grade(1:3,iel))))
    errlip = max(errlip, maxval(abs(grade(1:3,iel))))

    cvar_err(iel) = maxval(abs(grad(1:3,iel)-grade(1:3,iel)))

    cvar_err2(iel) = maxval(abs(grad2(1:3,iel)-grade(1:3,iel)))
  enddo

  sumerrl2 = errl2
  sumerrl2p = errl2p

  sumerrl1 = errl1
  sumerrl1p = errl1p

  nceltot = ncel

  if (irangp.ge.0) then
    call parsom(sumerrl2)
    call parsom(sumerrl2p)
    call parsom(sumerrl1)
    call parsom(sumerrl1p)
    call parmax(errli)
    call parmax(errlip)
    call parcpt(nceltot)
  endif

  if(irangp.le.0) then
    errl2 = sqrt(sumerrl2/sumerrl2p)
    errl1 = sumerrl1/sumerrl1p
    errli = errli/errlip
    write(errout,*) '# ncel^(1/3) | erreur L^1 | erreur L^2 | erreur L^inf'
    write(errout,*) nceltot**(1d0/3d0), errl1 , errl2, errli
    close(errout)
  endif

  ! Profile

  if(irangp.le.0) then
    open(impout,file='profile.dat')
    write(impout,*) &
      '# x | scalaire | scalaire th | dtdx | dtdx th'
  endif

  npoint = 300
  iel1 = -999
  irang1 = -999
  do ii = 1, npoint

    xyz(1) = float(ii-1)/float(npoint-1)*2d0-1d0
    xyz(2) = 0.5d0
    xyz(3) = 0.5d0

    call findpt(ncelet, ncel, xyzcen, xyz(1), xyz(2), xyz(3), iel, irangv)


    if ((iel.ne.iel1).or.(irangv.ne.irang1)) then
      iel1   = iel
      irang1 = irangv

      ! Set temporary variables for the process containing
      ! the point and then send it to other processes.
      if (irangp.eq.irangv) then
        x    = xyzcen(1,iel)
        sc   = cvar_scal(iel)
        scth = cvar_exa(iel)
        dtdx = grad(1,iel)
        dtth = grade(1,iel)
      endif
      ! Broadcast to other ranks in parallel
      if (irangp.ge.0) then
        call parall_bcast_r(irangv, x)
        call parall_bcast_r(irangv, sc)
        call parall_bcast_r(irangv, scth)
        call parall_bcast_r(irangv, dtdx)
        call parall_bcast_r(irangv, dtth)
      endif
      if(irangp.le.0) then
        write(impout,*) x, sc, scth, dtdx, dtth
      endif
    endif
  enddo


  if(irangp.le.0) close(impout)

endif
!----
! End
!----

return
end subroutine cs_f_user_extra_operations
