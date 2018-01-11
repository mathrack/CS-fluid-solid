!===============================================================================
! User source terms definition.
!
! 1) Momentum equation (coupled solver)
! 2) Species transport
! 3) Turbulence (k-epsilon, k-omega, Rij-epsilon, v2-f, Spalart-Allmaras)
!===============================================================================

!-------------------------------------------------------------------------------

!VERS

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

!>    User subroutine.
!> \brief    Additional right-hand side source terms for scalar equations (user
!>     scalars and specific physics scalars).
!>
!>
!> Usage
!> -----
!> The routine is called for each scalar, user or specific physisc. It is
!> therefore necessary to test the value of the scalar number iscal to separate
!> the treatments of the different scalars (if (iscal.eq.p) then ....).
!>
!> The additional source term is decomposed into an explicit part (crvexp) and
!> an implicit part (crvimp) that must be provided here.
!> The resulting equation solved by the code for a scalar f is:
!>
!>   \f[ \rho*volume*\frac{df}{dt} + .... = crvimp*f + crvexp \f]
!>
!>
!> Note that crvexp and crvimp are defined after the Finite Volume integration
!> over the cells, so they include the "volume" term. More precisely:
!>   - crvexp is expressed in kg.[scal]/s, where [scal] is the unit of the scalar
!>   - crvimp is expressed in kg/s
!>
!>
!> The crvexp and crvimp arrays are already initialized to 0 before entering the
!> the routine. It is not needed to do it in the routine (waste of CPU time).
!>
!> For stability reasons, Code_Saturne will not add -crvimp directly to the
!> diagonal of the matrix, but Max(-crvimp,0). This way, the crvimp term is
!> treated implicitely only if it strengthens the diagonal of the matrix.
!> However, when using the second-order in time scheme, this limitation cannot
!> be done anymore and -crvimp is added directly. The user should therefore test
!> the negativity of crvimp by himself.
!>
!> When using the second-order in time scheme, one should supply:
!>   - crvexp at time n
!>   - crvimp at time n+1/2
!>
!>
!> The selection of cells where to apply the source terms is based on a getcel
!> command. For more info on the syntax of the getcel command, refer to the
!> user manual or to the comments on the similar command \ref getfbr in the routine
!> \ref cs_user_boundary_conditions.

!> WARNING: If scalar is the temperature, the resulting equation
!>          solved by the code is:
!>
!>  rho*Cp*volume*dT/dt + .... = crvimp*T + crvexp
!>
!>
!> Note that crvexp and crvimp are defined after the Finite Volume integration
!> over the cells, so they include the "volume" term. More precisely:
!>   - crvexp is expressed in W
!>   - crvimp is expressed in W/K
!>

!>
!> STEP SOURCE TERMS
!>===================
!> In case of a complex, non-linear source term, say F(f), for scalar f, the
!> easiest method is to implement the source term explicitely.
!>
!>   df/dt = .... + F(f(n))
!>   where f(n) is the value of f at time tn, the beginning of the time step.
!>
!> This yields :
!>   crvexp = volume*F(f(n))
!>   crvimp = 0
!>
!> However, if the source term is potentially steep, this fully explicit
!> method will probably generate instabilities. It is therefore wiser to
!> partially implicit the term by writing:
!>
!>   df/dt = .... + dF/df*f(n+1) - dF/df*f(n) + F(f(n))
!>
!> This yields:
!>   crvexp = volume*( F(f(n)) - dF/df*f(n) )
!>   crvimp = volume*dF/df
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss terms
!> \param[in]     ncesmp        number of cells with mass source terms
!> \param[in]     iscal         index number of the current scalar
!> \param[in]     icepdc        index number of cells with head loss terms
!> \param[in]     icetsm        index number of cells with mass source terms
!> \param[in]     itypsm        type of mass source term for each variable
!> \param[in]                    (see \ref cs_user_mass_source_terms)
!> \param[in]     dt            time step (per cell)
!> \param[in]     ckupdc        head loss coefficient
!> \param[in]     smacel        value associated to each variable in the mass
!> \param[in]                    source terms or mass rate (see \ref
!>                               cs_user_mass_source_terms)
!> \param[out]    crvexp        explicit part of the source term
!> \param[out]    crvimp        implicit part of the source term
!______________________________________________________________________________!


subroutine ustssc &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   ckupdc , smacel ,                                              &
   crvexp , crvimp )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use parall
use period
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iscal

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision crvexp(ncelet), crvimp(ncelet)

! Local variables
integer :: iel
double precision :: x, y, z

!===============================================================================

do iel = 1, ncel

  x = xyzcen(1,iel)
  y = xyzcen(2,iel)
  z = xyzcen(3,iel)

  crvimp(iel) = 0.
  crvexp(iel) = volume(iel) * ( &
               3.d0*(dacos(-1.d0)**2)             &
             * sin (dacos(-1.d0) * x )                                  &  
             * sin (dacos(-1.d0) * y )                         &
             * sin (dacos(-1.d0) * z )                     &
             - 1.d0*(dacos(-1.d0)**2)                          &
             * cos (dacos(-1.d0) * x )                                  &  
             * cos (dacos(-1.d0) * y )                         &
             * sin (dacos(-1.d0) * z )                     &
             - 1.d0*(dacos(-1.d0)**2)                          &
             * sin (dacos(-1.d0) * x )                                  &
             * cos (dacos(-1.d0) * y )                         &
             * cos (dacos(-1.d0) * z ) )

enddo


return
end subroutine ustssc
