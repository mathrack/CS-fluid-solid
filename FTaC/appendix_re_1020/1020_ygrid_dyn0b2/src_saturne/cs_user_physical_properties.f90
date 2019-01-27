!-------------------------------------------------------------------------------

!VERS

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

!> \file cs_user_physical_properties.f90
!> \brief Definition of physical variable laws.
!>
!> usphyv
!> \brief Definition of physical variable laws.
!>
!> \section Warning
!>
!> It is \b forbidden to modify turbulent viscosity \c visct here
!> (a specific subroutine is dedicated to that: \ref usvist)
!>
!> - icp = 1 must <b> have been specified </b>
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
use entsor
use parall
use period
use ppppar
use ppthch
use ppincl
use field
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          mbrom

double precision dt(ncelet)

! Local variables

integer          ivart, iel, ifac
integer          ith, iscal, ii, ifcvsl
double precision vara, varb, varc, varam, varbm, varcm, vardm
double precision                   varal, varbl, varcl, vardl
double precision                   varac, varbc
double precision xvart

double precision, dimension(:), pointer :: coefap, coefbp
double precision, dimension(:), pointer :: bfpro_rom, cpro_rom
double precision, dimension(:), pointer :: cpro_viscl, cpro_vscalt, cpro_cp, cpro_beta
double precision, dimension(:), pointer :: cvar_scalt

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if (1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 0. Initializations to keep
!===============================================================================

!===============================================================================

!   The following examples should be adapted by the user
!   ====================================================

!  Each example is bounded by a test using .false. as a precaution.
!  Replace .false. by .true to activate the example.

!  It is recommended to keep only the minimum necessary in this file
!  (i.e. remove all unused example code)


!  example 1: variable density as a function of temperature
!  example 2: variable viscosity as a function of tempeprature
!  example 3: variable specific heat as a function of tempeprature
!  example 4: variable Lambda/CP as a function of temperature
!             for temperature or enthalpy
!  example 5: variable sclalars diffusivity as a function of temperature
!===============================================================================


!===============================================================================
!  Example 1: variable density as a function of temperature
!  =========
!    Below, we define the same density law
!    Values of this property must be defined at cell centers
!      (and optionally, at boundary faces).
!  ===================================================================

!    The test on .false. allows deactivating instructions (which are defined
!       only as a starting example)

if (.false.) then

  ! Position of variables, coefficients
  ! -----------------------------------

  ! --- Number of the thermal variable
  !       (and of its boundary conditions)
  !       To use user scalar 2 instead, write 'ivart = isca(2)'

  if (iscalt.gt.0) then
    ivart = isca(iscalt)
    call field_get_val_s(ivarfl(ivart), cvar_scalt)
  else
    write(nfecra,9010) iscalt
    call csexit (1)
  endif

  ! --- Pointers to density values

  call field_get_val_s(icrom, cpro_rom)
  call field_get_val_s(ibrom, bfpro_rom)

  ! --- Coefficients of laws chosen by the user
  !       Values given here are fictitious

  vara  = -4.0668d-3
  varb  = -5.0754d-2
  varc  =  1000.9d0

  ! Density at cell centers
  !------------------------
  ! law                    rho  = t  * ( a *  t +  b) +   c
  ! so      cpro_rom(iel) = xvart * (vara*xvart+varb) + varc

  ! Volumic thermal expansion coefficient
  !--------------------------------------
  ! law                     cpro_beta  = -1/rho * (d rho / d T)
  ! so cpro_beta(iel) = (-1.d0/cpro_rom(iel))*(2.d0*vara*xvart+varb)

  call field_get_val_s(iprpfl(ibeta), cpro_beta)

  do iel = 1, ncel
    xvart = cvar_scalt(iel)
    cpro_rom(iel) = xvart * (vara*xvart+varb) + varc
    cpro_beta(iel)= (-1.d0/cpro_rom(iel))*(2.d0*vara*xvart+varb)
  enddo


  ! Density at boundary faces
  !---------------------------

  ! By default, the value of rho at the boundary is the value taken
  !   at the center of adjacent cells. This is the recommended approach.
  ! To be in this case, nothing needs to be done:
  !   do not prescribe a value for bfpro_rom(ifac) and
  !   do not modify mbrom

  ! For users who do not wish to follow this recommendation, we
  !   note that the boundary temperature may be fictitious, simply
  !   defined so as to conserve a flux (this is especially the case
  !   at walls). The value of rho which is computed at the boundary
  !   when introducing this fictitious temperature in a physical law
  !   may thus be completely false (negative for example).

  ! If we wish to specify a law anyways:
  !                        rho  = t  * ( a *  t +  b) +   c
  ! so      bfpro_rom(ifac) = xvart * (vara*xvart+varb) + varc

  ! 't' being the temperature at boundary face centers, we may use the
  ! following lines of code (voluntarily deactived, as the must be used
  ! with caution):

  ! Note that when we prescribe the density at the boundary, it must be done
  ! at ALL boundary faces.
  !    ===

  if (.false.) then

    ! Boundary condition coefficients
    call field_get_coefa_s(ivarfl(ivart), coefap)
    call field_get_coefb_s(ivarfl(ivart), coefbp)

    ! Caution: mbrom = 1 is necessary for the law to be taken
    !                           into account.
    mbrom = 1

    do ifac = 1, nfabor

      ! ifabor(ifac) is the cell adjacent to the boundary face
      iel = ifabor(ifac)
      xvart = coefap(ifac) + cvar_scalt(iel)*coefbp(ifac)
      bfpro_rom(ifac) = xvart * (vara*xvart+varb) + varc
    enddo

  endif ! --- Test on .false.

endif ! --- Test on .false.


!===============================================================================
!  Example 2: variable viscosity as a function of temperature
!  =========
!    Below, we define the same viscosity law
!    Values of this property must be defined at cell centers
!  ===================================================================

!    The test on .false. allows deactivating instructions (which are defined
!       only as a starting example)

if (.false.) then

  ! Position of variables, coefficients
  ! -----------------------------------

  ! --- Number of the thermal variable
  !       To use user scalar 2 instead, write 'ivart = isca(2)'

  if (iscalt.gt.0) then
    ivart = isca(iscalt)
    call field_get_val_s(ivarfl(ivart), cvar_scalt)
  else
    write(nfecra,9010) iscalt
    call csexit(1)
  endif

  ! --- Molecular dynamic viscosity

  call field_get_val_s(iprpfl(iviscl), cpro_viscl)

  ! --- Coefficients of laws chosen by the user
  !       Values given here are fictitious

  varam = -3.4016d-9
  varbm =  6.2332d-7
  varcm = -4.5577d-5
  vardm =  1.6935d-3

  ! Molecular dynamic viscosity in kg/(m.s) at cell centers
  !--------------------------------------------------------
  ! law                    mu   = t * (t * (am * t + bm) + cm) + dm
  ! so      cpro_viscl(iel) = xvart*(xvart*(varam*xvart+varbm)+varcm)+vardm

  do iel = 1, ncel
    xvart = cvar_scalt(iel)
    cpro_viscl(iel) = xvart*(xvart*(varam*xvart+varbm)+varcm)+vardm
  enddo

endif ! --- Test on .false.


!===============================================================================
!  Example 3: specific heat as a function of temperature
!  =========
!    Below, we define the same viscosity law
!    Values of this property must be defined at cell centers
!  ===================================================================

!    The test on .false. allows deactivating instructions (which are defined
!       only as a starting example)

if (.false.) then

  ! Position of variables, coefficients
  ! -----------------------------------

  ! --- Number of the thermal variable
  !       To use user scalar 2 instead, write 'ivart = isca(2)'

  if (iscalt.gt.0) then
    ivart = isca(iscalt)
    call field_get_val_s(ivarfl(ivart), cvar_scalt)
  else
    write(nfecra,9010) iscalt
    call csexit (1)
  endif

  ! --- Specific heat

  if (icp.gt.0) call field_get_val_s(iprpfl(icp), cpro_cp)

  ! --- Stop if Cp is not variable

  if (icp.le.0) then
    write(nfecra,1000) icp
    call csexit (1)
  endif

  ! --- Coefficients of laws chosen by the user
  !       Values given here are fictitious

  varac = 0.00001d0
  varbc = 1000.0d0

  ! Specific heat in J/(kg.degrees) at cell centers
  !------------------------------------------------
  ! law                    cpro_cp  = ac * t + bm
  ! so          cpro_cp(iel) = varac*xvart + varbc

  do iel = 1, ncel
    xvart = cvar_scalt(iel)
    cpro_cp(iel) = varac*xvart + varbc
  enddo

endif ! --- Test on .false.


!======================================================================================
!  Example 4: Lambda/Cp a function of temperature for enthalpy or
!             Lambda    a function of temperature for temperature because Cp is put
!                       outside the divergence term
!  =========
!  Below, we define the same lambda/Cp ratio law (or lambda law if temperature is used)
!  Values of this property must be defined at cell centers
!  ====================================================================================

!    The test on .false. allows deactivating instructions (which are defined
!       only as a starting example)

if (.false.) then

  ! Position of variables, coefficients
  ! -----------------------------------

  ! --- Number of the thermal variable
  !       To use user scalar 2 instead, write 'ivart = isca(2)'

  if (iscalt.gt.0) then
    ivart = isca(iscalt)
    call field_get_val_s(ivarfl(ivart), cvar_scalt)
  else
    write(nfecra,9010) iscalt
    call csexit (1)
  endif

  ! --- Lambda/Cp of the thermal (or Lambda if temperature is used)

  call field_get_key_int(ivarfl(isca(iscalt)), kivisl, ifcvsl)

  if (ifcvsl.ge.0) then
    call field_get_val_s(ifcvsl, cpro_vscalt)
  else
    cpro_vscalt => NULL()
  endif

  ! --- Stop if Lambda/CP (or Lambda if temperature is used) is not variable

  if (ifcvsl.lt.0) then
    write(nfecra,1010) iscalt
    call csexit (1)
  endif

  ! if thermal variable is not temperature
  if (iscacp(iscal).le.0) then

    ! --- Specific heat

    if (icp.gt.0) call field_get_val_s(iprpfl(icp), cpro_cp)

    ! --- Coefficients of laws chosen by the user
    !       Values given here are fictitious

    varal = -3.3283d-7
    varbl =  3.6021d-5
    varcl =  1.2527d-4
    vardl =  0.58923d0

    ! Lambda/Cp in kg/(m.s) at cell centers
    !--------------------------------------
    ! law    Lambda/Cp = {t * (t * (al * t +  bl) + cl) + dl} / Cp
    ! so     cpro_vscalt(iel) &
    !             = (xvart*(xvart*(varal*xvart+varbl)+varcl)+vardl)/cp0

    ! We assume Cp has been defined previously.

    if (icp.le.0) then

      ! --- If Cp is uniform, we use cp0
      do iel = 1, ncel
        xvart = cvar_scalt(iel)
        cpro_vscalt(iel) =                                           &
             (xvart*(xvart*(varal*xvart+varbl)+varcl)+vardl)         &
             /cp0
      enddo

    else

      ! --- If Cp is not uniform, we use cpro_vscalt above
      do iel = 1, ncel
        xvart = cvar_scalt(iel)
        cpro_vscalt(iel) =                                           &
             (xvart*(xvart*(varal*xvart+varbl)+varcl)+vardl)         &
             /cpro_cp(iel)
      enddo

    endif

  ! if iscalt is temperature, the Cp division is not needed
  else

    ! --- Coefficients of laws chosen by the user
    !       Values given here are fictitious

    varal = -3.3283d-7
    varbl =  3.6021d-5
    varcl =  1.2527d-4
    vardl =  0.58923d0

    ! Lambda in W/(m.K) at cell centers
    !--------------------------------------
    ! law    Lambda = {t * (t * (al * t +  bl) + cl) + dl}
    ! so     cpro_vscalt(iel) &
    !             = (xvart*(xvart*(varal*xvart+varbl)+varcl)+vardl)

    do iel = 1, ncel
      xvart = cvar_scalt(iel)
      cpro_vscalt(iel) = (xvart*(xvart*(varal*xvart+varbl)+varcl)+vardl)
    enddo

  endif

endif ! --- Test on .false.


!===============================================================================
!  Example 5: Diffusivity as a function of temperature for user scalars
!  =========
!    Excluding:
!      - temperature, enthalpy (handled above)
!      - fluctuation variances (property equal to that of the associated scalar)
!
!    Below, we define the same diffusivity law for all scalars (except the
!      ones excluded above).
!    Values of this property must be defined at cell centers
!  ===================================================================

!    The test on .false. allows deactivating instructions (which are defined
!       only as a starting example)

if (.false.) then

  do ii = 1, nscaus ! Loop on scalars

    ! --- Number of user scalar 'ii' in the lsit of scalars
    iscal = ii

    ! --- If it is a thermal variable, it has already been handled above
    ith = 0
    if (iscal.eq.iscalt) ith = 1

    ! --- If the variable is a fluctuation, its diffusivity is the same
    !       as that of the scalar to which it is attached:
    !       there is nothing to do here, we move on to the next variable
    !       without setting cpro_vscalt(iel).

    ! We only handle here non-thermal variables which are not fluctuations
    if (ith.eq.0.and.iscavr(iscal).le.0) then

      ! Position of variables, coefficients
      ! -----------------------------------

      ! --- Number of the thermal variable
      !       To use user scalar 2 instead, write 'ivart = isca(2)'

      if (iscalt.gt.0) then
        ivart = isca(iscalt)
        call field_get_val_s(ivarfl(ivart), cvar_scalt)
      else
        write(nfecra,9010) iscalt
        call csexit (1)
      endif

      ! --- Scalar's Lambda

      call field_get_key_int(ivarfl(isca(iscal)), kivisl, ifcvsl)
      if (ifcvsl.ge.0) then
        call field_get_val_s(ifcvsl, cpro_vscalt)
      else
        cpro_vscalt => NULL()
      endif

      ! --- Stop if Lambda is not variable

      if (ifcvsl.lt.0) then
        write(nfecra,1010) iscal
        call csexit (1)
      endif

      ! --- Coefficients of laws chosen by the user
      !       Values given here are fictitious

      varal = -3.3283d-7
      varbl =  3.6021d-5
      varcl =  1.2527d-4
      vardl =  0.58923d0

      ! Lambda in kg/(m.s) at cell centers
      !--------------------------------------
      ! law    Lambda = {t * (t * (al * t +  bl) + cl) + dl}
      ! so     cpro_vscalt(iel) &
      !             = (xvart*(xvart*(varal*xvart+varbl)+varcl)+vardl)

      do iel = 1, ncel
        xvart = cvar_scalt(iel)
        cpro_vscalt(iel) = (xvart*(xvart*(varal*xvart+varbl)+varcl)+vardl)
      enddo

    endif ! --- Tests on 'ith' and 'iscavr'

  enddo ! --- Loop on scalars
endif ! --- Test on .false.


!===============================================================================

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    DONNEES DE CALCUL INCOHERENTES                          ',/,&
'@                                                            ',/,&
'@      usipsu indique que la chaleur specifique est uniforme ',/,&
'@        ICP = ',I10   ,' alors que                          ',/,&
'@      usphyv impose une chaleur specifique variable.        ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Modifier usipsu ou usphyv.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1010 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    DONNEES DE CALCUL INCOHERENTES                          ',/,&
'@                                                            ',/,&
'@    Pour le scalaire ', i10                                  ,/,&
'@      la diffusivite est uniforme alors que                 ',/,&
'@      usphyv impose une diffusivite variable.               ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Modifier usipsu ou usphyv.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9010 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    APPEL A csexit DANS LE SOUS PROGRAMME usphyv            ',/,&
'@                                                            ',/,&
'@    La variable dont dependent les proprietes physiques ne  ',/,&
'@      semble pas etre une variable de calcul.               ',/,&
'@    En effet, on cherche a utiliser la temperature alors que',/,&
'@      ISCALT = ',I10                                         ,/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Verifier le codage de usphyv (et le test lors de la     ',/,&
'@      definition de IVART).                                 ',/,&
'@    Verifier la definition des variables de calcul dans     ',/,&
'@      usipsu. Si un scalaire doit jouer le role de la       ',/,&
'@      temperature, verifier que ISCALT a ete renseigne.     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:  stop when computing physical quantities',/,       &
'@    =======',/,                                                 &
'@    Inconsistent calculation data',/,                           &
'@',/,                                                            &
'@      usipsu specifies that the specific heat is uniform',/,    &
'@        icp = ',i10   ,' while',/,                              &
'@      usphyv prescribes a variable specific heat.',/,           &
'@',/,                                                            &
'@    The calculation will not be run.',/,                        &
'@',/,                                                            &
'@    Modify usipsu or usphyv.',/,                                &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 1010 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:  stop when computing physical quantities',/,       &
'@    =======',/,                                                 &
'@    Inconsistent calculation data',/,                           &
'@',/,                                                            &
'@    For scalar', i10,/,                                         &
'@      the diffusivity is uniform while',/,                      &
'@      usphyv prescribes a variable diffusivity.',/,             &
'@',/,                                                            &
'@    The calculation will not be run.',/,                        &
'@',/,                                                            &
'@    Modify usipsu or usphyv.',/,                                &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 9010 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:  stop when computing physical quantities',/,       &
'@    =======',/,                                                 &
'@',/,                                                            &
'@    The variable on which physical properties depend does',/,   &
'@      seem to be a calculation variable.',/,                    &
'@    Indeed, we are trying to use the temperature while',/,      &
'@      iscalt = ',i10                                  ,/,&
'@',/,                                                            &
'@    The calculation will not be run.',/,                        &
'@',/,                                                            &
'@    Check the programming in usphyv (and the test when',/,      &
'@      defining ivart).',/,                                      &
'@    Check the definition of calculation variables in',/,        &
'@      usipsu. If a scalar should represent the,',/,             &
'@      temperature, check that iscalt has been defined',/,       &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

#endif

!----
! End
!----

return
end subroutine usphyv


!===============================================================================


subroutine uselph &
!================

 ( nvar   , nscal  ,                                              &
   mbrom  , izfppp ,                                              &
   dt     )

!===============================================================================
! FONCTION :
! --------

!   REMPLISSAGE DES VARIABLES PHYSIQUES POUR LE MODULE ELECTRIQUE

!     ----> Effet Joule
!     ----> Arc Electrique
!     ----> Conduction Ionique

!      1) Masse Volumique
!      2) Viscosite moleculaire
!      3) Chaleur massique Cp
!      4) Lambda/Cp moleculaire
!      4) Diffusivite moleculaire



! ATTENTION :
! =========


! Il est INTERDIT de modifier la viscosite turbulente VISCT ici
!        ========
!  (une routine specifique est dediee a cela : usvist)

! Pour le module electrique, toutes les proprietes physiques sont
!   supposees variables et contenues dans le tableau PROPCE
!   (meme si elles sont physiquement constantes)


! Remarques :
! ---------

! Cette routine est appelee au debut de chaque pas de temps

!    Ainsi, AU PREMIER PAS DE TEMPS (calcul non suite), les seules
!    grandeurs initialisees avant appel sont celles donnees
!      - dans usipsu :
!             . la masse volumique (initialisee a RO0)
!             . la viscosite       (initialisee a VISCL0)
!      - dans usiniv/useliv :
!             . les variables de calcul  (initialisees a 0 par defaut
!             ou a l
!             ou a la valeur donnee dans usiniv)

! On peut donner ici les lois de variation aux cellules
!     - de la masse volumique                      ROM    kg/m3
!         (et eventuellememt aux faces de bord     ROMB   kg/m3)
!     - de la viscosite moleculaire                VISCL  kg/(m s)
!     - de la chaleur specifique associee          CP     J/(kg degres)
!     - des "diffusivites" associees aux scalaires VISCLS kg/(m s)


! On dispose des types de faces de bord au pas de temps
!   precedent (sauf au premier pas de temps, ou les tableaux
!   ITYPFB et ITRIFB n'ont pas ete renseignes)


! Il est conseille de ne garder dans ce sous programme que
!    le strict necessaire.


! Cells identification
! ====================

! Cells may be identified using the 'getcel' subroutine.
! The syntax of this subroutine is described in the
! 'cs_user_boundary_conditions' subroutine,
! but a more thorough description can be found in the user guide.


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! mbrom            ! te ! <-- ! indicateur de remplissage de romb              !
! izfppp           ! te ! <-- ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module phys. part.                    !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use entsor
use ppppar
use ppthch
use ppincl
use elincl
use field
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          mbrom
integer          izfppp(nfabor)

double precision dt(ncelet)

! Local variables

integer          iel, ifcvsl
integer          mode

double precision tp
double precision xkr   , xbr
double precision rom0  , temp0 , dilar , aa    , bb    , cc
double precision srrom1, rhonp1
double precision, dimension(:), pointer :: cpro_rom
double precision, dimension(:), pointer :: cpro_viscl, cpro_vscalt
double precision, dimension(:), pointer :: cpro_vpoti, cpro_vpotr
double precision, dimension(:), pointer :: cpro_cp, cpro_temp
double precision, dimension(:), pointer :: cvar_scalt

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================

!===============================================================================
! 0 - INITIALISATIONS A CONSERVER
!===============================================================================

! --- Initialisation memoire


ipass = ipass + 1

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
!       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
!       ET NON CELLE DE LA BIBLIOTHEQUE
!===============================================================================


!     En Joule, on s'arrete : il faut que l'utilisateur
!       donne les proprietes physiques
if ( ippmod(ieljou).ge.1 ) then

  write(nfecra,9010)
  call csexit (1)

!     En Arc on continue car on a un fichier de donnees
!       Un message indique que l'utilisateur n'a rien fourni
elseif (ippmod(ielarc).ge.1) then

  if (ipass.eq.1) then
    write(nfecra,9011)
  endif

  return

endif

 9010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA DEFINITION DES PROP. PHYSIQUES   ',/,&
'@    =========                                               ',/,&
'@                      MODULE ELECTRIQUE                     ',/,&
'@                                                            ',/,&
'@     LE SOUS-PROGRAMME UTILISATEUR uselph DOIT ETRE COMPLETE',/,&
'@                                                            ',/,&
'@     Ce sous-programme utilisateur permet de definir les    ',/,&
'@       proprietes physiques. Il est indispensable.          ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9011 format(/,                                                   &
' Module arc electrique: pas d''intervention utilisateur pour ',/,&
'                          le calcul des proprietes physiques.',/)

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!     Message au premier passage pour indiquer que l'utilisateur a
!       rapatrie le sous-programme.
if (ipass.eq.1) then
  write(nfecra,1000)
endif

!===============================================================================
! 1 - EFFET JOULE
!===============================================================================

if (ippmod(ieljou).ge.1) then


!     Attention, dans les modules electriques, la chaleur massique, la
!       conductivite thermique et la conductivite electriques sont
!       toujours dans le tableau PROPCE
!       qu'elles soient physiquement variables ou non.

!       On n'utilisera donc PAS les variables
!          =====================
!                                cp0, visls0(iscalt)
!                                visls0(ipotr) et visls0(ipoti)

!       Informatiquement, ceci se traduit par le fait que
!                                icp > 0, et le numero associe au mot cle
!                                scalar_diffusivity_id est >= 0 pour
!                                les scalaires iscalt, ipotr, ipoti



!       Calcul de la temperature a partir de l'enthalpie
!       ------------------------------------------------

!       Ceci depend largement des choix utilisateur en
!         matiere de loi H-T (T en Kelvin)

!       On demande de fournir cette loi dans le sous programme usthht
!          (users/usthht.f90)
!           usthht fournit en particulier un exemple d'interpolation
!            a partir d'une tabulation utilisateur
!           usthht en mode T->H sera utilise pour l'initialisation
!            de l'enthalpie dans useliv.

!       mode = 1 : h => ivarfl(isca(ihm)) -> t => iprpfl(itemp)
  mode = 1

  call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt)
  call field_get_val_s(iprpfl(itemp), cpro_temp)

  do iel = 1, ncel
    call usthht (mode,cvar_scalt(iel),cpro_temp(iel))
  enddo


!       Masse volumique au centre des cellules
!       --------------------------------------

!     ATTENTION :
!     =========
!       Dans le module electrique effet Joule, on fournira
!       OBLIGATOIREMENT la loi de variation de la masse volumique ici
!       en renseignant cpro_rom(iel)
!       (meme si elle est uniforme ou constante).


!     Masse Vol : RO = ROM0 / (1+DILAR*(T-T0)
!         (Choudhary) semblable a Plard (HE-25/94/017)

!          avec sous-relaxation (sauf au premier pas de temps)

  temp0  = 300.d0
  rom0   = 2500.d0
  dilar  = 7.5d-5
  if (ntcabs.gt.1) then
    srrom1 = srrom
  else
    srrom1 = 0.d0
  endif

  call field_get_val_s(icrom, cpro_rom)
  do iel = 1, ncel
    rhonp1 = rom0 /                                               &
            (1.d0+ dilar * (cpro_temp(iel)-temp0) )
    cpro_rom(iel) = srrom1*cpro_rom(iel)+(1.d0-srrom1)*rhonp1
  enddo


!       Viscosite moleculaire dynamique en kg/(m s)
!        ------------------------------------------

!     ATTENTION :
!     =========
!       Dans le module electrique effet Joule, on fournira
!       OBLIGATOIREMENT la loi de variation de la viscosite ici
!       en renseignant cpro_viscl(iel)
!       (meme si elle est uniforme ou constante).


!     Viscosite : MU = EXP((AA/T-BB)-CC)
!          (Choudhary)
!      Plard (HE-25/94/017) ; limite a 1173K par C Delalondre

  call field_get_val_s(iprpfl(iviscl), cpro_viscl)
  aa     = 10425.d0
  bb     =   500.d0
  cc     =-6.0917d0

  do iel = 1, ncel
    if ( cpro_temp(iel) .gt. 1173.d0 ) then
      tp = cpro_temp(iel)
    else
      tp= 1173.d0
    endif
    cpro_viscl(iel) = exp( (aa/(tp-bb))+cc )
  enddo


!       Chaleur specifique J/(kg degres)
!       --------------------------------

!     ATTENTION :
!     =========
!       Dans le module electrique effet Joule, on fournira
!       OBLIGATOIREMENT la loi de variation de la chaleur massique ici
!       en renseignant PROPCE(IEL,IPCPP)
!       (meme si elle est uniforme ou constante).


!        CP = 1381 (Choudhary)
!          coherent avec Plard (HE-25/94/017)

  if (icp.gt.0) call field_get_val_s(iprpfl(icp), cpro_cp)
  do iel = 1, ncel
    cpro_cp(iel) = 1381.d0
  enddo


!       Lambda/Cp en kg/(m s)
!       ---------------------

!     ATTENTION :
!     =========
!       Dans le module electrique effet Joule, on fournira
!       OBLIGATOIREMENT la loi de variation de la conductivite ici
!       en renseignant PROPCE(IEL,IPCVSL)
!       (meme si elle est uniforme ou constante).


!         Lambda
!          On suppose Cp renseigne au prealable.

!          Plard (HE-25/94/017)

  call field_get_key_int(ivarfl(isca(iscalt)), kivisl, ifcvsl)
  call field_get_val_s(ifcvsl, cpro_vscalt)

  do iel = 1, ncel
    xbr = 85.25d0                                                      &
         -5.93d-2*(cpro_temp(iel)-tkelvi)                              &
         +2.39d-5*(cpro_temp(iel)-tkelvi)**2
    xkr = 16.d0*stephn*(1.4d0)**2*(cpro_temp(iel))**3                  &
         /(3.d0*xbr)

    cpro_vscalt(iel) = 1.73d0 + xkr
  enddo

! --- On utilise CP calcule dans cpro_vscalt ci dessus
  do iel = 1, ncel
    cpro_vscalt(iel) = cpro_vscalt(iel)/cpro_cp(iel)
  enddo


!       Conductivite electrique en S/m
!       ==============================

!     ATTENTION :
!     =========
!       Dans le module electrique effet Joule, on fournira
!       OBLIGATOIREMENT la loi de variation de la conductivite ici
!       en renseignant PROPCE(IEL,IPCSIG)
!       (meme si elle est uniforme ou constante).


!         SIGMA  (Plard HE-25/94/017)

  call field_get_key_int(ivarfl(isca(ipotr)), kivisl, ifcvsl)
  call field_get_val_s(ifcvsl, cpro_vpotr)
  do iel = 1, ncel
    cpro_vpotr(iel) =                                                  &
         exp(7.605d0-7200.d0/cpro_temp(iel))
  enddo

!     La conductivite electrique pour le potentiel imaginaire est
!       toujours implicitement prise egale a la conductivite
!       utilisee pour le potentiel reel.
!       IL NE FAUT PAS la renseigner.

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
!     Conductivite electrique imaginaire :
!     La conductivite reelle et imaginaire sont dans le meme tableau.
!       Ce choix est fait en dur dans varpos.
!       Les pointeurs pour les deux existent quand meme.
!     Sinon, on pourrait faire ceci :
  if (1.eq.0) then
    if ( ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4 ) then
      call field_get_key_int(ivarfl(isca(ipoti)), kivisl, ifcvsl)
      call field_get_val_s(ifcvsl, cpro_vpoti)
      do iel = 1, ncel
        cpro_vpoti(iel) = cpro_vpotr(iel)
      enddo
    endif
  endif
!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!       Diffusivite variable a l'exclusion de l'enthalpie et du potentiel
!       -----------------------------------------------------------------
!     Pour le moment, il n'y a pas d'autres scalaires et
!                                                  on ne fait donc rien

endif

!===============================================================================
! 2 - ARC ELECTRIQUE
!===============================================================================

!     Les proprietes physiques sont a priori fournies par fichier
!       de donnees. IL n'y a donc rien a faire ici.

!      IF ( IPPMOD(IELARC).GE.1 ) THEN
!      ENDIF


!===============================================================================
! 3 - CONDUCTION IONIQUE
!===============================================================================

!     CETTE OPTION N'EST PAS ACTIVABLE

!--------
! Formats
!--------

 1000 format(/,                                                   &
' Module electrique: intervention utilisateur pour        ',/,    &
'                      le calcul des proprietes physiques.',/)

!----
! End
!----

return
end subroutine uselph


!===============================================================================


subroutine ussmag &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   ckupdc , smacel ,                                              &
   smagor , mijlij , mijmij )

!===============================================================================
! FONCTION :
! --------

! MODIFICATION UTILISATEUR DE LA CONSTANTE DE SMAGORINSKY
! DANS LE CAS DE L'UTILISATION D'UN MODELE DYNAMIQUE

!              SMAGOR = Mij.Lij / Mij.Mij

! EN FAIT, DES MOYENNES LOCALES DU NUMERATEUR ET DU DENOMINATEUR
! SONT REALISEES AVANT L'APPEL A USSMAG, SOIT

!              SMAGOR = < Mij.Lij > / < Mij.Mij >

! DANS CET ROUTINE, Mij.Lij ET Mij.Mij SONT PASSES EN ARGUMENT
! AVANT LA MOYENNE LOCALE.
! DANS L'EXEMPLE CI-DESSOUS ON REALISE UNE MOYENNE LOCALE DU
! RAPPORT.
!              SMAGOR = < Mij.Lij / Mij.Mij >

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. cs_user_mass_source_terms)     !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
! smagor(ncelet)   ! tr ! <-- ! constante de smagorinsky dans le cas           !
!                  !    !     ! d'un modlele dynamique                         !
! mijlij(ncelet    ! tr ! <-- ! mij.lij avant moyenne locale                   !
! mijmij(ncelet    ! tr ! <-- ! mij.mij avant moyenne locale                   !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use cstnum
use optcal
use cstphy
use entsor
use parall
use field
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision smagor(ncelet), mijlij(ncelet), mijmij(ncelet)

! Local variables
integer          iel
double precision, allocatable, dimension(:) :: w1

! Local variables used for spatial averaging
! 
! Variables with parameter attribute

! Variables with save attribute
integer, save :: ny, nxnz
integer, allocatable, dimension(:), save :: map_3d_1d
double precision, allocatable, dimension(:), save :: ygrid
! Other variables
double precision :: xmin, zmin, yy!, eps
character(len=64) :: buffer
integer :: nyloc
integer, allocatable, dimension(:) :: lstelt, ncell_3d_1d
double precision, allocatable, dimension(:) :: ygridloc, w1d


!===============================================================================
! 1.  INITIALISATION
!===============================================================================

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

! Allocate work arrays
allocate(w1(ncel))
allocate(w1d(ny))

!===============================================================================
! 2.  MOYENNE SPATIALE SUR DIRECTIONS HOMOGENES
!
!     Dans le cas ou l'utilisateur souhaite utilise le voisinage
!       etendu, il est fortement conseille de passer le mode de
!       calcul des gradients en IMRGRA = 2, afin de conserver
!       la totalite du voisinage etendu. En effet, le calcul de
!       moyenne locale est generalement degradee en voisinage reduit
!       (IMRGRA = 3).

!===============================================================================

! Sum numerator in homogeneous directions
! Put sum in temporary 3d array w1
w1d = 0.
do iel = 1, ncel
  w1d( map_3d_1d(iel) ) = w1d( map_3d_1d(iel) ) + mijlij(iel)
enddo
if (irangp.ge.0) call parrsm(ny,w1d)
w1d = w1d / nxnz
do iel = 1, ncel
  w1(iel) = w1d( map_3d_1d(iel) )
enddo

! Sum denominator in homogeneous directions
! Put sum in temporary 1d array w1d
w1d = 0.
do iel = 1, ncel
  w1d( map_3d_1d(iel) ) = w1d( map_3d_1d(iel) ) + mijmij(iel)
enddo
if (irangp.ge.0) call parrsm(ny,w1d)
w1d = w1d / nxnz

!     On calcule le rapport
do iel = 1, ncel
  if (abs(w1d( map_3d_1d(iel) )).le.epzero) then
    w1(iel) = smagmx**2
  else
    w1(iel) = w1(iel)/w1d( map_3d_1d(iel) )
  endif
enddo

! Free memory
deallocate(w1)
deallocate(w1d)

!----
! End
!----

return
end subroutine ussmag


!===============================================================================

!===============================================================================
! Purpose:
! -------

!> usvima
!> \brief User subroutine dedicated the use of ALE
!>  (Arbitrary Lagrangian Eulerian Method): fills mesh viscosity arrays.
!>
!> This subroutine is called at the beginning of each time step.
!>
!> Here one can modify mesh viscosity value to prevent cells and nodes
!> from huge displacements in awkward areas, such as boundary layer for example.
!>
!> IF variable IORTVM = 0, mesh viscosity modeling is isotropic therefore VISCMX
!> array only needs to be filled.
!> IF variable IORTVM = 1, mesh viscosity modeling is orthotropic therefore
!> all arrays VISCMX, VISCMY and VISCMZ need to be filled.
!>
!> Note that VISCMX, VISCMY and VISCMZ arrays are initialized at the first time step
!> to the value of 1.
!>
!> \section usvima_cell_id Cells identification
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
!> \param[out]    viscmx        mesh viscosity in X direction
!> \param[out]    viscmy        mesh viscosity in Y direction
!> \param[out]    viscmz        mesh viscosity in Z direction
!_______________________________________________________________________________

subroutine usvima &
 ( nvar   , nscal  ,                                              &
   dt     ,                                                       &
   viscmx , viscmy , viscmz )

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
use entsor
use parall
use period
use albase
use field
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet)
double precision viscmx(ncelet), viscmy(ncelet), viscmz(ncelet)

! Local variables

integer          iel
double precision rad, xr2, xcen, ycen, zcen

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if (1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
!  1. Example :
!       One gives a huge mesh viscosity value to the cells included in a circle
!       with a radius 'rad' and a center '(xcen, ycen, zcen)'

!     In general it appears quite much easier to fill mesh viscosity arrays at
!     the beginning of the calculations basing on the initial geometry.

!  The test on .false. allows deactivating instructions (which are defined
!     only as a starting example)

if (.false.) then

  if (ntcabs.eq.0) then
    rad = (1.d-3)**2
    xcen  = 1.d0
    ycen  = 0.d0
    zcen  = 0.d0

    do iel = 1, ncel
      xr2 = (xyzcen(1,iel)-xcen)**2 + (xyzcen(2,iel)-ycen)**2       &
           + (xyzcen(3,iel)-zcen)**2
      if (xr2.lt.rad) viscmx(iel) = 1.d10
    enddo

    ! 2. In case of orthotropic mesh viscosity modeling one can choose
    !    to submit nodes to a lower stress in Z direction

    if (iortvm.eq.1) then
      do iel = 1, ncel
        viscmy(iel) = viscmx(iel)
        viscmz(iel) = 1.d0
      enddo
    endif

  endif

endif ! --- Test on .false.

!----
! End
!----

return
end subroutine usvima

!===============================================================================
! Purpose:
! -------

!> usatph
!> \brief User subroutine dedicated to modifie physical properties of the
!>        atmospheric module
!>
!> This subroutine is called at beginning of each time step at the end of
!> atphyv.
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!_______________________________________________________________________________

subroutine usatph

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use pointe
use numvar
use optcal
use cstphy
use entsor
use parall
use period
use albase
use field
use mesh

!===============================================================================

implicit none

! Arguments

! Local variables

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if (1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!----
! Formats
!----

!----
! End
!----

return
end subroutine usatph
