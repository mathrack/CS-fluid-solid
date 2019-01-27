!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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
! Function:
! ---------

!> \file viswal.f90
!>
!> \brief Compute the turbulent viscosity for the WALE LES model
!>
!> The turbulent viscosity is:
!> \f$ \mu_T = \rho (C_{wale} L)^2 * \dfrac{(\tens{S}:\tens{Sd})^{3/2}}
!>                                         {(\tens{S} :\tens{S})^(5/2)
!>                                         +(\tens{Sd}:\tens{Sd})^(5/4)} \f$
!> with \f$ \tens{S}  = \frac{1}{2}(\gradt \vect{u} + \transpose{\gradt \vect{u}})\f$
!> and  \f$ \tens{Sd} = \deviator{(\symmetric{(\tens{S}^2)})}\f$
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------

subroutine viswal

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use entsor
use parall
use mesh
use field
use field_operator

!===============================================================================

implicit none

! Arguments

! Local variables

return

end subroutine
