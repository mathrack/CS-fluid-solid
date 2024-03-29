!-------------------------------------------------------------------------------

!                      Code_Saturne version 5.0.7-patch
!                      --------------------------
! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

integer          ifac
integer          ii
integer          ilelt  , nlelt

double precision xfor(3)
double precision, dimension(:,:), pointer :: bfprp_for

integer, allocatable, dimension(:) :: lstelt

integer :: iel
double precision :: y, tbulk, myvoltot
double precision, dimension(nscal) :: tbulk_nsca
double precision, dimension(:), pointer :: mytemp
!< [loc_var_dec]

!===============================================================================

!===============================================================================
! Initialization
!===============================================================================

if (ipstdv(ipstfo).ge.1) call field_get_val_v(iforbr, bfprp_for)

! Allocate a temporary array for cells or interior/boundary faces selection
allocate(lstelt(nfabor))!max(ncel,nfac,nfabor)))

!===============================================================================
! Example: compute global efforts on a subset of faces
!===============================================================================

! If efforts have been calculated correctly:

if (ipstdv(ipstfo).ge.1) then

  ! Bottom wall
  do ii = 1, ndim
    xfor(ii) = 0.d0
  enddo
  call getfbr('plane[0,1,0,1,0.0001]', nlelt, lstelt)
  do ilelt = 1, nlelt
    ifac = lstelt(ilelt)
    do ii = 1, ndim
      xfor(ii) = xfor(ii) + bfprp_for(ii, ifac)
    enddo
  enddo
  if (irangp.ge.0) then
    call parrsm(ndim,xfor)
  endif
  if (irangp.le.0) then
    write(nfecra,*) 'cedric output : friction bot wall = ',xfor
  endif

  ! Top wall
  do ii = 1, ndim
    xfor(ii) = 0.d0
  enddo
  call getfbr('plane[0,1,0,-1,0.0001]', nlelt, lstelt)
  do ilelt = 1, nlelt
    ifac = lstelt(ilelt)
    do ii = 1, ndim
      xfor(ii) = xfor(ii) + bfprp_for(ii, ifac)
    enddo
  enddo
  if (irangp.ge.0) then
    call parrsm(ndim,xfor)
  endif
  if (irangp.le.0) then
    write(nfecra,*) 'cedric output : friction top wall = ',xfor
  endif

endif ! ipstdv(ipstfo).ge.1

! Deallocate the temporary array
deallocate(lstelt)

! Compute and print bulk temperatures
if (nscal.ge.1) then

  ! Compute bulk temperature for each scalar and total volume
  do ii = 1, nscal

    ! Total volume
    select case (ii)
      case (1) ! diric + neuma
        myvoltot = 0.d0
        do iel = 1,ncel
          y = xyzcen(2, iel)
          if (-1.lt.y .and. y.lt.1.) then
            myvoltot = myvoltot + volume(iel)
          endif
        enddo
        if (irangp.ge.0) then
          call parsom(myvoltot)
        endif

      case (3) ! fluid-solid coupling
        myvoltot = 0.d0
        do iel = 1,ncel
          myvoltot = myvoltot + volume(iel)
        enddo
        if (irangp.ge.0) then
          call parsom(myvoltot)
        endif
    end select

    ! Bulk temperature
    call field_get_val_s(ivarfl(isca(ii)),mytemp)
    tbulk = 0.d0
    select case (ii)
      case (1 : 2) ! diric + neuma
        do iel = 1,ncel
          y = xyzcen(2, iel)
          if (-1.lt.y .and. y.lt.1.) then
            tbulk = tbulk + volume(iel)*mytemp(iel)
          endif
        enddo

      case (3 : ) ! fluid-solid coupling
        do iel = 1,ncel
          tbulk = tbulk + volume(iel)*mytemp(iel)
        enddo
    end select
    if (irangp.ge.0) then
      call parsom(tbulk)
    endif
    tbulk_nsca(ii) = tbulk/myvoltot

  enddo ! ii = 1, nscal

  ! Remove bulk part except for nscal=1 => diric
  if (nscal.ge.2) then
  do ii = 2, nscal
    call field_get_val_s(ivarfl(isca(ii)),mytemp)
    select case (ii)
      case (2) ! neuma
        do iel = 1,ncel
          y = xyzcen(2, iel)
          if (-1.lt.y .and. y.lt.1.) then
            mytemp(iel) = mytemp(iel) - tbulk_nsca(ii)
          endif
        enddo

      case (3 : ) ! fluid-solid coupling
        do iel = 1,ncel
          mytemp(iel) = mytemp(iel) - tbulk_nsca(ii)
        enddo
    end select

  enddo ! ii = 2, nscal
  endif ! nscal.ge.2

  ! Listing
  if (mod(ntcabs,10).eq.0) then

    ! Print bulk temperatures in fluid + solid
    if (irangp.le.0) then
      write(nfecra,*) 'cedric output : bulk fluid temperature = ',ttcabs,tbulk_nsca
    endif

    ! Print bulk temperatures in solid
    if (nscal.ge.3) then

      ! Compute solid volume
      myvoltot = 0.
      do iel = 1,ncel
        y = xyzcen(2, iel)
        if (1.lt.y .or. y.lt.-1.) then
          myvoltot = myvoltot + volume(iel)
        endif ! 1.lt.y .or. y.lt.-1.
      enddo ! iel = 1,ncel
      if (irangp.ge.0) then
        call parsom(myvoltot)
      endif ! irangp.ge.0

      ! Compute bulk solid temp
      do ii = 3, nscal
        tbulk = 0.
        call field_get_val_s(ivarfl(isca(ii)),mytemp)
        do iel = 1,ncel
          y = xyzcen(2, iel)
          if (1.lt.y .or. y.lt.-1.) then
            tbulk = tbulk + volume(iel)*mytemp(iel)
          endif
        enddo
        if (irangp.ge.0) then
          call parsom(tbulk)
        endif
        tbulk_nsca(ii) = tbulk/myvoltot
      enddo ! ii = 3, nscal

      ! Print
      if (irangp.le.0) then
        write(nfecra,*) 'cedric output : bulk solid temperature = ',ttcabs,tbulk_nsca(3:)
      endif

    endif ! nscal.ge.3

  endif ! mod(ntcabs,10).eq.0

endif ! nscal.ge.1

return
end subroutine cs_f_user_extra_operations
