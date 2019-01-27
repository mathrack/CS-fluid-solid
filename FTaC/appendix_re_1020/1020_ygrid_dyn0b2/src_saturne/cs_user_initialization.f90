!-------------------------------------------------------------------------------

!                      Code_Saturne version 4.0.4
!                      --------------------------
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
! Purpose:
! -------

!> \file cs_user_initialization.f90
!>
!> \brief Initialize variables
!>
!> This subroutine is called at beginning of the computation
!> (restart or not) before the loop time step.
!>
!> This subroutine enables to initialize or modify (for restart)
!> unkown variables and time step values.
!>
!> \c rom and \c viscl values are equal to \c ro0 and \c viscl0 or initialize
!> by reading the restart file.
!> variable diffusivity and cp variables (when there are defined) have no value
!> excepted if they are read from a restart file.
!>
!> Modification of the behaviour law of physical quantities (rom, viscl,
!> viscls, cp) is not done here. It is the purpose of the user subroutine
!> \ref usphyv (in cs_user_physical_properties.f90)
!>
!> \section cs_user_initialization_cell_id Cells identification
!>
!> Cell value field ids
!>
!> - Density:                        \c iprpfl(irom)
!> - Dynamic molecular viscosity:    \c iprpfl(iviscl)
!> - Turbulent viscosity:            \c iprpfl(ivisct)
!> - Specific heat:                  \c iprpfl(icp)
!> - Diffusivity(lambda):            \c field_get_key_int(ivarfl(isca(iscal)),
!>                                      kivisl, ...)
!>
!> Cells may be identified using the \ref getcel subroutine.
!> The syntax of this subroutine is described in the
!> \ref cs_user_boundary_conditions subroutine,
!> but a more thorough description can be found in the user guide.
!
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

subroutine cs_user_initialization &
 ( nvar   , nscal  ,                                              &
   dt     )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use atincl
use ctincl
use elincl
use ppcpfu
use cs_coal_incl
use cs_fuel_incl
use mesh
use field
use turbomachinery

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables

integer          iel
double precision y, noise
double precision, dimension(3,ncel) :: my_rand_var

double precision, dimension(:,:), pointer :: cvar_iu

integer, allocatable, dimension(:) :: lstelt

!===============================================================================


!===============================================================================
! Initialization
!===============================================================================

allocate(lstelt(ncel)) ! temporary array for cells selection

if (isuite.eq.0) then

  ! Generate random numbers
  noise=0.1
  if (noise.ne.0.d0) then
    call random_number(my_rand_var)
  else
    my_rand_var = 0.d0
  endif

  ! Initialize velocity
  call field_get_val_v(ivarfl(iu), cvar_iu)
  do iel=1,ncel
    y=xyzcen(2,iel)
    cvar_iu(1,iel)=(1.5+noise*my_rand_var(1,iel))*(1.0-y**2)
    cvar_iu(2,iel)=(0.0+noise*my_rand_var(2,iel))*(1.0-y**2)
    cvar_iu(3,iel)=(0.0+noise*my_rand_var(3,iel))*(1.0-y**2)
  enddo
endif
!--------
! Formats
!--------

!----
! End
!----

deallocate(lstelt) ! temporary array for cells selection

return
end subroutine cs_user_initialization
