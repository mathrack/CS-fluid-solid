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

!> \file cs_user_modules.f90
!>
!> \brief User-defined module: it allows to create any user array.
!>
!> See \subpage cs_user_modules for examples.
!>
!> This file is compiled before all other user Fortran files.
!> To ensure this, it must not be renamed.
!>
!> The user may define an arbitrary number of modules here, even though
!> only one is defined in the example.
!
!> \cond DOXYGEN_SHOULD_SKIP_THIS

!-------------------------------------------------------------------------------

module case_setup

  ! Amplitude of the scalar
  double precision, parameter :: Tm = 1.d0 ! x < 0
  double precision, parameter :: Tp = 2.d0 ! x > 0

  ! Number of harmonics
  integer, parameter :: nm = 1 ! x < 0
  integer, parameter :: np = 1 ! x > 0

  ! Conductivity
  ! Enforce continuity of heat flux at x = 0 and t > 0
  double precision, parameter :: num = 1.d0 ! x < 0
  double precision, parameter :: nup = num * Tm * dble(nm) / (Tp * dble(np))

  ! Density
  double precision, parameter :: rhom = 1.d0
  double precision, parameter :: rhop = (rhom * Tm / dble(nm)) * dble(np) / Tp

  ! Pi ** 2
  double precision, parameter :: hipi = dacos(-1.d0)
  double precision, parameter :: pipi = hipi*hipi

end module case_setup

!-------------------------------------------------------------------------------

!> (DOXYGEN_SHOULD_SKIP_THIS) \endcond
