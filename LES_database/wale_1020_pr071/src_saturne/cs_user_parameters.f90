!-------------------------------------------------------------------------------

!                      Code_Saturne version 5.0.3
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

!> \file cs_user_parameters.f90
!>
!> \brief User subroutines for input of calculation parameters (Fortran modules).
!>        These subroutines are called in all cases.
!>
!>  See \subpage f_parameters for examples.
!>
!>   If the Code_Saturne GUI is used, this file is not required (but may be
!>   used to override parameters entered through the GUI, and to set
!>   parameters not accessible through the GUI).
!>
!>   Several routines are present in the file, each destined to defined
!>   specific parameters.
!>
!>   To modify the default value of parameters which do not appear in the
!>   examples provided, code should be placed as follows:
!>   - usipsu   for numerical and physical options
!>   - usipes   for input-output related options
!>
!>   As a convention, "specific physics" defers to the following modules only:
!>   pulverized coal, gas combustion, electric arcs.
!>
!>   In addition, specific routines are provided for the definition of some
!>   "specific physics" options.
!>   These routines are described at the end of this file and will be activated
!>   when the corresponding option is selected in the usppmo routine.
!-------------------------------------------------------------------------------


!===============================================================================

!> \brief User subroutine for input of model selection parameters.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]      ixmlpu       indicates if the XML file from the GUI is used
!>                              used (1: yes, 0: no
!> \param[in, out] iturb        turbulence model
!> \param[in, out] itherm       thermal model
!> \param[in, out] iale         ale module
!> \param[in, out] ivofmt       vof method
!> \param[in, out] icavit       cavitation model
!______________________________________________________________________________!

subroutine usipph &
 ( ixmlpu, iturb , itherm, iale , ivofmt, icavit )

!===============================================================================
! Module files
!===============================================================================

use entsor, only: nfecra ! No other module should appear here
use optcal, only: irijco ! No other module should appear here

!===============================================================================

implicit none

! Arguments

integer ixmlpu
integer iturb, itherm, iale, ivofmt, icavit

  iturb = 42

return
end subroutine usipph


!===============================================================================

!> \brief User subroutine for the input of additional user parameters.
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nmodpp         number of active specific physics models
!______________________________________________________________________________!

subroutine usipsu &
 ( nmodpp )

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use dimens
use numvar
use optcal
use cstphy
use entsor
use parall
use period
use ihmpre
use albase
use ppppar
use ppthch
use ppincl
use coincl
use cpincl
use field
use cavitation
use post
use rotation
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer nmodpp

! Local variables

integer ii, kscmin, kscmax, ifcvsl
type(var_cal_opt) :: vcopt

  ipstdv(ipstfo) = 1
  ileaux = 1
  idtvar = 0
  ntmabs = 3000000
  dtref  = 0.002d0
  imvisf = 1
  imrgra = 0

  call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt)
  vcopt%blencv = 1.0d0
  call field_set_key_struct_var_cal_opt(ivarfl(iu), vcopt)
  if (nscaus.ge.1) then
    do ii = 1, nscaus
      call field_get_key_struct_var_cal_opt(ivarfl(isca(ii)), vcopt)
      vcopt%blencv = 1.0d0
      call field_set_key_struct_var_cal_opt(ivarfl(isca(ii)), vcopt)
    enddo
  endif

  ro0    = 1d0
  viscl0 = 2.d0/14124.d0
  cp0    = 0.71d0

! --- Variable diffusivity field id (ifcvsl>=0) or constant
! --- Variable density field id (ifcvsl>=0)
!     With ifcvsl = 0, the field will be added automatically
  do ii = 1, 2
    visls0(ii) = viscl0/cp0
    ifcvsl = -1
    call field_set_key_int(ivarfl(isca(ii)), kivisl, ifcvsl)
    call field_set_key_int(ivarfl(isca(ii)), kromsl, ifcvsl)
  enddo
  do ii = 3, nscaus
    visls0(ii) = viscl0/cp0
    ifcvsl = 0
    call field_set_key_int(ivarfl(isca(ii)), kivisl, ifcvsl)
    call field_set_key_int(ivarfl(isca(ii)), kromsl, ifcvsl)
  enddo

return
end subroutine usipsu


!===============================================================================

!> \brief User subroutine for the input of additional user parameters for
!>        input/output.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nmodpp       number of active specific physics models
!______________________________________________________________________________!

subroutine usipes &
 ( nmodpp )

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use dimens
use numvar
use optcal
use cstphy
use entsor
use field
use parall
use period
use ihmpre
use post
use ppppar
use ppthch
use ppincl
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer nmodpp

  ntlist = 500

return
end subroutine usipes
