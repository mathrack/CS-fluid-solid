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


subroutine usppmo &
!================
 ( ixmlpu )


!===============================================================================
! Purpose:
! -------

!> \brief User subroutine.

!> Define the use of a specific physics amongst the following:
!>   - combustion with gas / coal / heavy fuel oil
!>   - compressible flows
!>   - electric arcs
!>   - atmospheric modelling
!>   - radiative transfer
!>   - cooling towers modelling
!>
!>    Only one specific physics module can be activated at once.


!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ixmlpu        indicates if the XML file from the GUI is used
!>                              (1 : yes, 0 : no)
!______________________________________________________________________________!

!===============================================================================
! Module files
!===============================================================================

use paramx
use entsor
use cstphy
use ppppar
use ppthch
use ppincl
use ppcpfu
use coincl
use radiat
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer ixmlpu

return
end subroutine usppmo


!===============================================================================


subroutine usipph &
!================
 ( ixmlpu, iturb , itherm, iale , icavit )

!===============================================================================
! Purpose:
! --------

!> \brief User subroutine for input of parameters.

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
!> \param[in, out] icavit       cavitation model
!______________________________________________________________________________!

!===============================================================================
! Module files
!===============================================================================

use entsor, only: nfecra ! No other module should appear here
use optcal, only: irijco ! No other module should appear here

!===============================================================================

implicit none

! Arguments

integer ixmlpu
integer iturb, itherm, iale, icavit

! Local variables

!===============================================================================

!===============================================================================


!>    In this subroutine, only the parameters which already appear may
!>    be set, to the exclusion of any other.
!>
!>    If we are not using the Code_Saturne GUI:
!>    All the parameters which appear in this subroutine must be set.
!>
!>    If we are using the Code_Saturne GUI:
!>    parameters protected by a test of the form:
!>
!>      if (ixmlpu.eq.0) then
!>         ...
!>      endif
!>
!>    should already have been defined using the GUI, so only
!>    experts should consider removing the test and adapting them here.

!===============================================================================

!----
! Formats
!----

return
end subroutine usipph


!===============================================================================


subroutine usipsu &
!================

 ( nmodpp )

!===============================================================================
! Purpose:
! -------

!> \brief User subroutine for the input of additional user parameters.
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nmodpp         number of active specific physics models
!______________________________________________________________________________!

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
use rotation
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer nmodpp

!===============================================================================

!>  This subroutine allows setting parameters
!>  which do not already appear in the other subroutines of this file.
!>
!>  It is possible to add or remove parameters.
!>  The number of physical properties and variables is known here.

!===============================================================================

integer f_id, iscal

iccvfg = 1 ! frozen velocity

ntmabs = 1 ! Number of time steps
dtref = 0.0001d0/dble(ntmabs) ! Time step

! Reconstruction gradients

imvisf = 1 ! harmonic interpolation for face viscosity
imrgra = 1 ! 1 <=> LSQ gradient reconstruction

call field_get_id('exa', f_id)
call field_set_key_int(f_id , keyvis, 1)
call field_set_key_int(f_id , keylog, 1)

call field_get_id('error', f_id)
call field_set_key_int(f_id , keyvis, 1)
call field_set_key_int(f_id , keylog, 1)

call field_get_id('error2', f_id)
call field_set_key_int(f_id , keyvis, 1)
call field_set_key_int(f_id , keylog, 1)

call field_get_id('err_grad', f_id)
call field_set_key_int(f_id , keyvis, 1)
call field_set_key_int(f_id , keylog, 1)

call field_get_id('err_grad2', f_id)
call field_set_key_int(f_id , keyvis, 1)
call field_set_key_int(f_id , keylog, 1)

! --- Variable density field id (ifcvsl>=0)
!     With ifcvsl = 0, the field will be added automatically, and later calls to
!       field_get_key_int(ivarfl(isca(iscal)), kromsl, ifcvsl)
!       will return its id.
ro0 = 1.
viscl0 = 1.
cp0 = 1.
do iscal = 1, nscaus
  visls0(iscal) = viscl0
  call field_set_key_int(ivarfl(isca(iscal)), kromsl, 0)
  call field_set_key_int(ivarfl(isca(iscal)), kivisl, 0)
enddo

!----
! Formats
!----

return
end subroutine usipsu


!===============================================================================


subroutine usipes &
!================

 ( nmodpp )


!===============================================================================
! Purpose:
! --------

!> \brief User subroutine for the input of additional user parameters for
!>        input/output.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nmodpp       number of active specific physics models
!______________________________________________________________________________!

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
use ppppar
use ppthch
use ppincl

!===============================================================================

implicit none

! Arguments

integer nmodpp

!===============================================================================

!>     This subroutine allows setting parameters
!>     which do not already appear in the other subroutines of this file.
!>
!>     It is possible to add or remove parameters.
!>     The number of physical properties and variables is known here.

!===============================================================================

return
end subroutine usipes


!===============================================================================


subroutine usati1
!================


!===============================================================================
! Purpose:
! --------

!> \brief Initialize non-standard calculation options for the atmospheric version.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use ppppar
use atincl
use atsoil
use atchem
use atimbr
use siream

!===============================================================================

implicit none

!===============================================================================

!----
! End
!----

return
end subroutine usati1


!===============================================================================
! Purpose:
! -------
!
!> 1. Additional Calculation Options
!>    a. Density Relaxation
!>
!> 2. Physical Constants
!>    a.Dynamic Diffusion Coefficient
!>    b.Constants of the chosen model (EBU, Libby-Williams, ...)
!
!> This routine is called:
!>
!>
!>  - Eddy Break Up pre-mixed flame
!>  - Diffusion flame in the framework of ``3 points'' rapid complete chemistry
!>  - Libby-Williams pre-mixed flame
!>  - Lagrangian module coupled with pulverized coal:
!>    Eulerian combustion of pulverized coal and
!>    Lagrangian transport of coal particles
!>  - Pulverised coal combustion
!>  - Fuel (oil) combustion
!
!===============================================================================

subroutine cs_user_combustion


!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use parall
use ihmpre
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use cs_coal_incl
use cs_fuel_incl
use radiat

!===============================================================================

implicit none

!----
! End
!----

return
end subroutine cs_user_combustion


!===============================================================================

!===============================================================================
!  Purpose  :
!  -------
!> \brief User subroutines for input of calculation parameters,
!>        and to initialize variables used for radiative transfer module
!
!-------------------------------------------------------------------------------
subroutine cs_f_user_radiative_transfer_param

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use entsor
use optcal
use cstphy
use parall
use period
use ppppar
use radiat

!===============================================================================

implicit none

return

end subroutine cs_f_user_radiative_transfer_param


!===============================================================================

subroutine uscfx1
!================


!===============================================================================
! Purpose:
! -------

!> \brief User subroutine.

!> Initialize non standard options for the compressible flow scheme such
!> as the variability of the thermal conductivity and the volume viscosity.
!> Their values can be given in the subroutine \ref uscfx2 .

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!

!===============================================================================
! Module files
!===============================================================================

use paramx
use ihmpre
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use parall
use period
use ppppar
use ppthch
use ppincl
use field

!===============================================================================

implicit none

! Arguments

!===============================================================================

!===============================================================================


!----
! End
!----

return
end subroutine uscfx1


!===============================================================================


subroutine uscfx2
!================


! Purpose:
! -------

!> \brief User subroutine.
!>
!> Set values for the reference volumic viscosity, the reference
!> conductivity and the molar mass for compressible flow.
!>
!> Initialize non standard options for the compressible flow scheme such
!> as the hydrostatic equilibrium.
!>
!> In addition to options set in the user subroutine \ref uscfx1 (or in
!> the GUI): this subroutine allows to set a switch to indicate if the
!> molecular viscosity is constant, its values being given in the user
!> subroutine \ref usipsu .


!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!

!===============================================================================
! Module files
!===============================================================================

use paramx
use ihmpre
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use ppppar
use ppthch
use ppincl

!===============================================================================

implicit none

! Arguments

!===============================================================================


!----
! End
!----

return
end subroutine uscfx2


!===============================================================================


subroutine cs_user_cooling_towers
!================


!===============================================================================
! Purpose:
! -------

!> \brief Definition of cooling tower model and exchange zones

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use parall
use period
use ppppar
use ppthch
use ppincl
use ctincl

!===============================================================================

implicit none

!===============================================================================

!===============================================================================


!----
! End
!----

return
end subroutine cs_user_cooling_towers


!===============================================================================

subroutine user_darcy_ini1
!========================

!===============================================================================
!  Purpose:
!  --------

!> \brief User routine for definition of computation parameters dealing with darcy module

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!

!===============================================================================
! Module files
!===============================================================================

use ihmpre, only: iihmpr
use entsor
use darcy_module

!===============================================================================

implicit none

!===============================================================================


!----
! End
!----

return

end subroutine user_darcy_ini1

!===============================================================================
