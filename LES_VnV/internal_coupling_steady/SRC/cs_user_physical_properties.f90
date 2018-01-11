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

!> \file cs_user_physical_properties.f90
!>
!> \brief Definition of physical variable laws.
!>
!> See \subpage physical_properties for examples.
!>

!===============================================================================
!> usphyv
!>
!> \brief Definition of physical variable laws.
!>
!> \section Warning
!>
!> It is \b forbidden to modify turbulent viscosity \c visct here
!> (a specific subroutine is dedicated to that: \ref usvist)
!>
!> - icp = 0 must <b> have been specified </b>
!>    in \ref usipsu if we wish to define a variable specific heat
!>    cpro_cp (otherwise: memory overwrite).
!>
!> - the kivisl field integer key (scalar_diffusivity_id)
!>    must <b> have been specified </b>
!>    in \ref usipsu if we wish to define a variable viscosity
!>    \c viscls.
!>
!>
!> \remarks
!>  - This routine is called at the beginning of each time step
!>    Thus, <b> AT THE FIRST TIME STEP </b> (non-restart case), the only
!>    values initialized before this call are those defined
!>      - in the GUI or  \ref usipsu (cs_user_parameters.f90)
!>             - density    (initialized at \c ro0)
!>             - viscosity  (initialized at \c viscl0)
!>      - in the GUI or \ref cs_user_initialization
!>             - calculation variables (initialized at 0 by defaut
!>             or to the value given in the GUI or in \ref cs_user_initialization)
!>
!>  - We may define here variation laws for cell properties, for:
!>     - density:                                    rom    kg/m3
!>     - density at boundary faces:                  romb   kg/m3)
!>     - molecular viscosity:                        cpro_viscl  kg/(m s)
!>     - specific heat:                              cpro_cp     J/(kg degrees)
!>     - diffusivities associated with scalars:      cpro_vscalt kg/(m s)
!>
!> \b Warning: if the scalar is the temperature, cpro_vscalt corresponds
!> to its conductivity (Lambda) in W/(m K)
!>
!>
!> The types of boundary faces at the previous time step are available
!>   (except at the first time step, where arrays \c itypfb and \c itrifb have
!>   not been initialized yet)
!>
!> It is recommended to keep only the minimum necessary in this file
!>   (i.e. remove all unused example code)
!>
!>
!> \section usphyv_cell_id Cells identification
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
!> \param[in]     mbrom         indicator of filling of romb array
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________

subroutine usphyv &
 ( nvar   , nscal  ,                                              &
   mbrom  ,                                                       &
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
use ppincl
use field
use mesh
use lagran
use cs_c_bindings

use case_setup, only : num, nup

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          mbrom

double precision dt(ncelet)

! Local variables

!< [loc_var_dec]
integer          iel, iscal, ifcvsl

double precision, dimension(:), pointer :: cvar_diff, cvar_diff2
double precision x


!< [loc_var_dec]

!===============================================================================

!===============================================================================
! 0. Initializations to keep
!===============================================================================

!< [init]
iscal = 1
! --- Stop if the cvar_diffsivity is not variable
call field_get_key_int(ivarfl(isca(iscal)), kivisl, ifcvsl)
if (ifcvsl.lt.0) then
   write(nfecra,1010) iscal
   call csexit (1)
endif

call field_get_val_s(ifcvsl, cvar_diff)
call field_get_key_int(ivarfl(isca(2)), kivisl, ifcvsl)
call field_get_val_s(ifcvsl, cvar_diff2)

do iel = 1, ncel
   x = xyzcen(1, iel)
   if(x.lt.0d0) then
      cvar_diff(iel) = num
   else
      cvar_diff(iel) = nup
   endif
   cvar_diff2(iel) = cvar_diff(iel)
enddo

!< [example_1]

!===============================================================================


!----
! End
!----

 1010 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:  stop when computing physical quantities',/,       &
'@    =======',/,                                                 &
'@    Inconsistent calculation data',/,                           &
'@',/,                                                            &
'@    For scalar', i10,/,                                         &
'@      the cvar_diffsivity is uniform while',/,                      &
'@      usphyv prescribes a variable diffusivity.',/,             &
'@',/,                                                            &
'@    The calculation will not be run.',/,                        &
'@',/,                                                            &
'@    Modify usipsu or usphyv.',/,                                &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

return
end subroutine usphyv