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

!> \file iniini.f90
!> \brief Commons default initialization before handing over the user.
!>
!------------------------------------------------------------------------------

subroutine iniini

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use dimens
use numvar
use optcal
use cstphy
use rotation
use entsor
use pointe
use albase
use alstru
use alaste
use parall
use period
use ihmpre
use cplsat
use ppincl
use ppcpfu
use mesh
use field
use cavitation
use darcy_module
use radiat
use cs_nz_condensation, only: nzones

!===============================================================================

implicit none

! Local variables

integer          ii, jj, iscal, iprop, iest
integer          istr

!===============================================================================

!===============================================================================
! 0. STOCKAGE DES ARGUMENTS ET IMPRESSIONS INITIALES
!===============================================================================

write(nfecra, 900)

#if defined(_CS_LANG_FR)

 900  format(/,                                                   &
'===============================================================',&
/,/,                                                        &
'                   PREPARATION DU CALCUL                     ',/,&
'                   =====================                     ',/,&
                                                                /,&
                                                                /,&
' =========================================================== ',/,&
                                                                /,&
                                                                /)

#else

 900  format(/,                                                   &
'===============================================================',&
/,/,                                                        &
'                   CALCULATION PREPARATION'                   ,/,&
'                   ======================='                   ,/,&
                                                                /,&
                                                                /,&
' ===========================================================' ,/,&
                                                                /,&
                                                                /)

#endif



!===============================================================================
! 0. Global field keys
!===============================================================================

call field_get_key_id("label", keylbl)
call field_get_key_id('log', keylog)
call field_get_key_id('post_vis', keyvis)
call field_get_key_id("post_id", keyipp)

call field_get_key_id("inner_mass_flux_id", kimasf)
call field_get_key_id("boundary_mass_flux_id", kbmasf)

call field_get_key_id("scalar_diffusivity_id", kivisl)
call field_get_key_id("scalar_diffusivity_ref", kvisl0)

call field_get_key_id("gradient_weighting_id", kwgrec)

call field_get_key_id("source_term_prev_id", kstprv)

icrom = -1
ibrom = -1

ipori = -1
iporf = -1

!===============================================================================
! 1. Map Fortran pointers to C global data
!===============================================================================

call time_step_init
call time_step_options_init
call thermal_model_init
call turb_model_init
call turb_rans_model_init
call turb_les_model_init
call wall_functions_init
call stokes_options_init
call physical_constants_init
call fluid_properties_init
call space_disc_options_init
call piso_options_init
call wall_reference_values_init
call turb_reference_values_init

!===============================================================================
! 2. ENTREES SORTIES entsor.f90
!===============================================================================

! ---> Impressions standard
!         -10000 : non initialise : voir modini
do ii = 1, nvarmx
  iwarni(ii) = -10000
enddo

! ---> NFECRA vaut 6 par defaut ou 9 en parallele (CSINIT)

!    Methode des vortex : on utilise la meme unite
!                  Pour le fichier de donnees, on specifie l'unite ici, mais le
!                  nom est laisse dans usvort. On utilise 20 et pas 11
!                  car en cas d'entree multiple, les deux fichiers doivent etre
!                  ouverts en meme temps dans VORINI
impmvo = 11
impdvo = 20

! ---> Fichier aval

!     NTSUIT : Periode de sauvegarde du fichier suite
!            (il est de toutes facons ecrit en fin de calcul)
!              -1   : a la fin seulement
!              0    : par defaut (4 fois par calcul)
!              > 0  : periode

ntsuit = 0

!    Methode des vortex : on utilise la meme unite
!                  et le format ascii obligatoirement
!                  (pas de detection automatique, fichiers de taille faible)
impvvo = 20

!    Fichier listing Lagrangien

implal = 80
ficlal = 'listla'
ntlal  = 1

! ---> Fichier thermochinie
!        FPP : utilisateur
!        JNF : Janaf
!        Les deux fichiers peuvent partager la meme unite
!          puisqu'ils sont lus l'un a pres l'autre.
!      En prime, INDJON (janaf=1 ou non=0)

impfpp = 25
ficfpp = 'define_ficfpp_in_usppmo'

indjon = 1

! ---> Fichiers module atmospherique
impmet = 26
ficmet = 'meteo'

! ---> Fichiers historiques

!     EMPHIS : EMPlacement
!     PREHIS : PREfixe
!     EXTHIS : EXTension
!     IMPUSH : Unite fichiers specifiques ushist
!     FICUSH : Nom   fichiers specifiques ushist
!     IMPSTH : fichier stock + unite d'ecriture des variables
!              des structures mobiles

impsth(1) = 30
impsth(2) = 31

emphis = 'monitoring/'
prehis = 'probes_'

do ii = 1, nushmx
  impush(ii) = 32+ii
  if (irangp .le. 0) then
    write(ficush(ii),'(1a3,i3.3)')'ush',ii
  else
    write(ficush(ii),'(1a3,i3.3,1a3,i4.4)')'ush',ii,'.n_',irangp+1
  endif
enddo

! tplfmt : time plot format (1: .dat, 2: .csv, 3: both)
! ncapt  : nombre de sondes total (limite a ncaptm)
! nthist : periode de sortie (> 0 ou -1 (jamais))
! frhist : frequence de sortie, en secondes (prioritaire sur nthist si > 0)
! nthsav : periode de sauvegarde (> 0 (fichiers ouverts et refermes) ou -1 )
! ihisvr : nb de sonde et numero par variable (-999 non initialise)
! ihistr : indicateur d'ecriture des historiques des structures
!          mobiles internes (=0 ou 1)
! ncapt  : nombre de sondes total (limite a ncaptm)
! nodcap : element correspondant aux sondes
! ndrcap : rang du processus contenant nodcap (parallelisme)
! xyzcap : position demandee des sondes
! tplflw : time plot flush wall-time interval (none if <= 0)

tplfmt = 1
ncapt = 0

nthist = 1
frhist = -1.d0
nthsav = -1

do ii = 1, nvppmx
  do jj = 1, ncaptm+1
    ihisvr(ii ,jj) = -999
  enddo
enddo

ihistr = 0

do ii = 1, ncaptm
  nodcap(ii) = 1
  ndrcap(ii) = -1
enddo

do ii = 1, ncaptm
  xyzcap(1,ii) = 0.d0
  xyzcap(2,ii) = 0.d0
  xyzcap(3,ii) = 0.d0
enddo

tplflw = -1.d0

! ---> Fichiers Lagrangiens

!     IMPLA1 : Unite fichier specifique Lagrangien
!     IMPLA2 : Unite fichier specifique Lagrangien
!     IMPLA3 : Unite fichier SCRATCH pour stockage temporaire
!     IMPLA4 : Unite fichier SCRATCH pour stockage temporaire
!     IMPLA5 : Unite d'ecriture des variables associees

impla1 = 50
impla2 = 51
impla3 = 52
impla4 = 53

impla5(1)  = 54
impla5(2)  = 55
impla5(3)  = 56
impla5(4)  = 57
impla5(5)  = 58
impla5(6)  = 59
impla5(7)  = 60
impla5(8)  = 61
impla5(9)  = 62
impla5(10) = 63
impla5(11) = 64
impla5(12) = 65
impla5(13) = 66
impla5(14) = 67
impla5(15) = 68

! ---> Fichiers utilisateurs

do ii = 1, nusrmx
  impusr(ii) = 69+ii
  if (irangp .le. 0) then
    write(ficusr(ii),'(1a4,i2.2)')'usrf',ii
  else
    write(ficusr(ii),'(1a4,i2.2,1a3,i4.4)')'usrf',ii,'.n_',irangp+1
  endif
enddo

! ---> Sorties listing

!   COMMUNES
!     IPP*   : Pointeurs de reperage des variables pour les sorties
!              1 pointe sur une case poubelle et constitue donc une
!              bonne initialisation
!     NOMVAR : Nom des variables
!     ITRSVR : numero de variable si IPP correspond a une variable resolue (p,u,k...)
!              0 si IPP correspond a une variable annexe (cp, mut...)ou a rien
!     NTLIST : periode d'ecriture
!       ( -1 : dernier pas de temps : > 0 : periode)

do ii = 1, npromx
  ipppro(ii) = 1
enddo
ippdt        = 1
ipptx        = 1
ippty        = 1
ipptz        = 1

ntlist = 1

! ---> Post traitement automatique (bord)

do ii = 1, 5
  ipstdv(ii) = 0
enddo

! ---> CPU
!      TMARUS : marge (Arret du calcul avant limite CPU)
!        Si TMARUS negatif, le code calcule une marge seul
!        Sinon, il utilise TMARUS (donnee en secondes)

tmarus = -1.d0


! Ici entsor.f90 est completement initialise

!===============================================================================
! 3. DIMENSIONS DE dimens.f90 (GEOMETRIE, sauf NCELBR)
!===============================================================================

! Geometry

ncel   = 0
ncelet = 0
nfac   = 0
nfabor = 0
nfml   = 0
nnod   = 0
lndfac = 0
lndfbr = 0

ncelgb = 0
nfacgb = 0
nfbrgb = 0
nsomgb = 0

! By default, assume no periodicity
iperio = 0
iperot = 0

! Get mesh metadata.

call ledevi(iperio, iperot)
!==========

call tstjpe(iperio, iperot)
!==========

!===============================================================================
! 4. DIMENSIONS de dimens.f90 (PHYSIQUE)
!===============================================================================

! --- Nombre de scalaires, de scalaires a diffusivite
!            variable, de variables, de proprietes

nscal  = 0
nscaus = 0
nscapp = 0
nscasp = 0
nvar   = 0
nproce = 0

!===============================================================================
! 5. POSITION DES VARIABLES DE numvar.f90
!===============================================================================

! --- Initialize mappings of field ids

do ii = 1, nvarmx
  ivarfl(ii) = -1
enddo

do ii = 1, npromx
  iprpfl(ii) = -1
enddo

idtten = -1

! --- Variables de calcul resolues (RTP, RTPA)

ipr    = 0
iu     = 0
iv     = 0
iw     = 0
ik     = 0
iep    = 0
ir11   = 0
ir22   = 0
ir33   = 0
ir12   = 0
ir13   = 0
ir23   = 0
iphi   = 0
ifb    = 0
ial    = 0
inusa  = 0

do iscal = 1, nscamx
  isca  (iscal) = 0
  iscapp(iscal) = 0
  iscasp(iscal) = 0
enddo

! --- Initialisation par defaut des commons pour la physique particuliere

call ppinii
!==========

! --- Proprietes physiques au sens large

do iprop  = 1, npromx
  ipproc(iprop) = iprop
enddo

irom   = 0
iroma  = 0
iviscl = 0
ivisct = 0
icour  = 0
ifour  = 0
icp    = 0
iprtot = 0

ibeta  = 0

! --- Ici tout numvar est initialise.

!===============================================================================
! 6. OPTIONS DU CALCUL : TABLEAUX DE optcal.f90
!===============================================================================

! --- Definition des equations
!       (convection-diffusion instationnaire, avec dirichlet
!        sauf pour la pression, diffusion instationnaire)
!        IDIFFT en particulier multiplie la diffusion turbulente
!           (quand elle est activee par le modele)

do ii = 1, nvarmx
  iconv (ii) = 1
  istat (ii) = 1
  idiff (ii) = 1
  idifft(ii) = 1
  idften(ii) = 1
  iswdyn(ii) = 0
enddo

! --- Schema en temps

!     NOTER BIEN que les valeurs de THETA pour
!       rho, visc, cp, flux de masse, termes sources
!       sont renseignees a partir des indicateurs I..EXT et ISTMPF
!       Elles ne sont pas completees directement par l'utilisateur
!       Les tests dans l'algo portent sur ces indicateurs


!   -- Schema en temps (regroupera les options suivantes)
!     = 1 Standard ordre 1
!     = 2 Standard ordre 2
!     si on veut, on peut en rajouter.
!     Cette variable conditionne toutes les autres de maniere automatique,
!     pour definir des schemas coherents dans modini.
ischtp = -999

!   -- Variables
!     Pour un schema centre (en n+1/2) il faut prendre theta = 0.5
!     On devrait pouvoir faire de l'explicite pur avec theta = 0 mais
!       ce point reste a voir
!     Ce theta sert aux termes de convection diffusion (partie implicitee
!       d'ordinaire)
!     Il est applique sous la forme (1-theta) ancien + theta nouveau
!       (ce n'est pas une extrapolation, contrairement aux termes sources)
do ii = 1, nvarmx
  thetav(ii) =-999.d0
enddo

!   -- Flux de masse (-999 = non initialise)
!     = 1 Standard d'ordre 1 (THETFL = -999 inutile)
!     = 0 Explicite (THETFL = 0)
!     = 2 Ordre deux (THETFL = 0.5)
istmpf = -999
thetfl =-999.d0

!   -- Termes sources Navier Stokes
!     Pour les termes sources explicites en std, I..EXT definit
!       l'extrapolation -theta ancien + (1+theta) nouveau
!     = 0 explicite
!     = 1 extrapolation avec theta = 1/2
!     = 2 extrapolation avec theta = 1
!       0 implique pas de reservation de tableaux
!       1 et 2 sont deux options equivalentes, la difference etant faite
!       uniquement au moment de fixer theta
!     Pour les termes sources implicites en std, I..EXT definit
!       la mise a l'ordre 2 ou non avec le thetav de la variable associee
!     = 0 implicite (std)
!     > 0 utilisation du thetav
!     Noter cpdt que le TS d'acc. masse n'est pas regi par I..EXT
!       (il suit bilsc2)
isno2t = -999
thetsn =-999.d0
!   -- Termes sources Grandeurs turbulentes
!     Pour les termes sources explicites en std, I..EXT definit
!       l'extrapolation -theta ancien + (1+theta) nouveau
!     = 0 explicite
!     = 1 extrapolation avec theta = 1/2
!     = 2 extrapolation avec theta = 1
!       0 implique pas de reservation de tableaux
!       1 et 2 sont deux options equivalentes, la difference etant faite
!       uniquement au moment de fixer theta
!     Pour les termes sources implicites en std, I..EXT definit
!       la mise a l'ordre 2 ou non avec le thetav de la variable associee
!     = 0 implicite (std)
!     > 0 utilisation du thetav
!     Noter cpdt que le TS d'acc. masse n'est pas regi par I..EXT
!       (il suit bilsc2)
isto2t = -999
thetst = -999.d0

! Backward differential time scheme order
do ii = 1, nvarmx
  ibdtso(ii) = 1
enddo

!    -- Proprietes physiques
!     I..EXT definit l'extrapolation -theta ancien + (1+theta) nouveau
!     = 0 explicite
!     = 1 extrapolation avec theta = 1/2
!     = 2 extrapolation avec theta = 1
!       0 implique pas de reservation de tableaux
!       1 et 2 sont deux options equivalentes, la difference etant faite
!       uniquement au moment de fixer theta
!     INIT.. =1 indique que la variable a ete proprement initialisee (dans un
!       fichier suite portant les valeurs adaptees)

!     Masse volumique
iroext = -999
thetro = -999.d0
initro = 0
!     Viscosite totale
iviext = -999
thetvi = -999.d0
initvi = 0
!     Chaleur specifique
icpext = -999
thetcp = -999.d0
initcp = 0

!   -- Convergence point fixe vitesse pression
epsup  = 1.d-5

!   -- Tab de travail pour normes de navsto
xnrmu0 = 0.d0
xnrmu  = 0.d0

!   -- Nb d'iter point fixe vitesse pression
nterup = 1

do iscal = 1, nscamx

!   -- Termes sources des scalaires
!     Pour les termes sources explicites en std, I..EXT definit
!       l'extrapolation -theta ancien + (1+theta) nouveau
!     = 0 explicite
!     = 1 extrapolation avec theta = 1/2
!     = 2 extrapolation avec theta = 1
!       0 implique pas de reservation de tableaux
!       1 et 2 sont deux options equivalentes, la difference etant faite
!       uniquement au moment de fixer theta
!     Pour les termes sources implicites en std, I..EXT definit
!       la mise a l'ordre 2 ou non avec le thetav de la variable associee
!     = 0 implicite (std)
!     > 0 utilisation du thetav
!     Noter cpdt que le TS d'acc. masse n'est pas regi par I..EXT
!       (il suit bilsc2)
  isso2t(iscal) = -999
  thetss(iscal) =-999.d0

!    -- Proprietes physiques
!     I..EXT definit l'extrapolation -theta ancien + (1+theta) nouveau
!     = 0 explicite
!     = 1 extrapolation avec theta = 1/2
!     = 2 extrapolation avec theta = 1
!       0 implique pas de reservation de tableaux
!       1 et 2 sont deux options equivalentes, la difference etant faite
!       uniquement au moment de fixer theta
!     INIT.. =1 indique que la variable a ete proprement initialisee (dans un
!       fichier suite portant les valeurs adaptees)

  ivsext(iscal) = -999
  thetvs(iscal) = -999.d0
  initvs(iscal) = 0

enddo



! --- Schema convectif
!       (a decider, centre-upwind si ordre 2,  test de pente a decider)

do ii = 1, nvarmx
  blencv(ii) = -999.d0
enddo
do ii = 1, nvarmx
  ischcv(ii) = 1
  isstpc(ii) = -999
enddo

! Method to compute interior mass flux due to ALE mesh velocity
! default: based on cell center mesh velocity
iflxmw = 0

! --- Reconstruction des gradients
!       On donne les valeurs par defaut
!       Pour la methode de limitation, on decidera plus tard
!         selon les choix de l'utilisateur
!       On n'active pas l'extrapolation des gradients par defaut
!         meme pour la pression (par securite : moins stable sur certains cas)

imrgra = 0
anomax = -grand*10.d0

do ii = 1, nvarmx
  nswrgr(ii) = 100
  nswrsm(ii) = -999
  imligr(ii) = -999
  ircflu(ii) = 1
enddo

do ii = 1, nvarmx
  epsrgr(ii) = 1.d-5
  climgr(ii) = 1.5d0
  extrag(ii) = 0.d0
enddo

! --- Solveurs iteratifs
!       La methode de resolution sera choisie selon les equations
!       La valeur de epsilon relatif est tres faible
!         (1.D-5 pourrait suffire)
!       On met IDIRCL a 1 pour toutes les variables. Pour toutes les
!       variables sauf pression et fb (en v2f) on a ISTAT=1, le
!       decalage de diagonale ne sera pas active. Pour la pression
!       on decalera la diagonale si necessaire. Pour fb, on sait qu'il
!       y a un autre terme diagonal (meme si ISTAT=0), donc IDIRCL
!       sera mis a 0 dans varpos.

do ii = 1, nvarmx
  idircl(ii) = 1
  ndircl(ii) = 0
  epsilo(ii) = -999.d0
  epsrsm(ii) = -999.d0
enddo

! --- Restarted calculation
!       By default, non-restarted calculation
!       Write auxiliary restart file by default
!       Read auxiliary restart file by default (in case of restarted calculation)
!       The match between new scalars and old scalars will be established later
!         (GUI, cs_user_parameters.f90, and lecamo)
!       The restart indicator of the 1D wall thermal model is initialized by default
!         to -1, to force the user to set it in uspt1d.
!       The same goes for the vortex method restart indicator.
!       The same goes for the synthetic turbulence method restart indicator.

isuite = 0
iecaux = 1
ileaux = 1
isuit1 = -1
isuivo = -1
isuisy = -1

! --- Reperage du temps

ntpabs = 0
ntcabs = ntpabs
ntmabs = 10
ntinit = 2

inpdt0 = 0

ttpabs = 0.d0
ttcabs = ttpabs

ttpmob = 0.d0
ttcmob = ttpmob

! --- Marche en temps
!       Par defaut pas de temps uniforme et constant,
!         sans coef multiplicatif
!       Dans les cas a pas de temps variable, par defaut
!         COUMAX = FOUMAX = 0.5 et variation max 10%
!       Les autres grandeurs sont a renseigner par l'utilisateur
!       Pour DTMIN et DTMAX, si on impose DTMIN > DTMAX,
!         les bornes sont prises egales a +/-1D12

idtvar = 0

coumax = 1.d0
foumax = 10.d0
cflmmx = 0.99d0
dtmin  = -grand*10.d0
dtmax  = -grand*10.d0
varrdt = 0.1d0

dtref  = -grand*10.d0

do ii = 1, nvarmx
  cdtvar(ii) = 1.d0
  relaxv(ii) =-999.d0
enddo
relxst = 0.7d0


!     Par defaut, pas de limitation du pas de temps liee aux
!     effets de densite

iptlro = 0

! --- Thermique

! No thermal scalar by default
itherm = 0
itpscl = 0
iscalt =-1

! No enthalpy for the gas (combustion) by default
ihgas = -1

! --- Turbulence
!     Le modele de turbulence devra etre choisi par l'utilisateur
!     En fait on n'a pas besoin d'initialiser ITYTUR (cf. varpos)
!     On suppose qu'il n'y a pas de temperature

iturb  =-999
itytur =-999

! Parfois, IGRHOK=1 donne des vecteurs non physiques en paroi
!        IGRHOK = 1
igrhok = 0
igrake = 1
ideuch =-999
iwallt = 0
ilogpo = 1
iclkep = 0
ikecou =-999
irijnu = 0
irijrb = 0
irijec = 0
igrari = 1
idifre = 1
iclsyr = 1
iclptr = 0
idries =-1

! --- Rotation/curvature correction of turbulence models
!     Unactivated by default
!     Correction type (itycor) is set in fldvar
irccor = 0
itycor = -999

! --- Turbulent diffusion model for second moment closure
!     Set to 1 (Daly and Harlow model) by default
idirsm = 1

! --- Stokes


! --- Take porosity into account

iporos = 0

! --- Algorithm to take into account the density variation in time

!     by default:
!     ----------
!      - the thermodynamic pressure (pther) is initialized with p0 = p_atmos
!      - the maximum thermodynamic pressure (pthermax) is initialized with -1
!        (no maximum by default, this term is used to model a venting effect when
!         a positive value is given by the user)

pther  = -1.d0
pthermax= -1.d0


! --- Cavitation module (not activated by default)
!       -1: module not activated
!        0: no vaporization/condensation model
!        1: Merkle's model
icavit = -1

! --- Radiative Transfert (not activated by default)
iirayo = 0

! ---  Choice the way to compute the condensation exchange coefficient (hcond)
!      associated to the condensation source term
!      ( not activated by default icophc =0)

icophc = 0

! ---  Choice the way to compute the thermal exchange coefficient (hpcond)
!      associated to the heat transfer of the condensation to the cooling wall
!      ( not activated by default icophg =0)

icophg = 0

! ---  Choice the way to compute to compute the wall temperature at
!      the solid/fluid interface coupled with condensation to the wall
!      ( not activated by default itag1d =0)

itag1d = 0

! --- Initialize the zones number for the condensation modelling
nzones = -1

! ---  Choice the way to compute to compute the wall temperature at
!      the solid/fluid interface coupled with condensation to the
!      metalli mass strucutres wall
!      ( not activated by default itagms =0)

itagms = 0

! --- Interpolation face des viscosites
!     = 0 ARITHMETIQUE
!     = 1 HARMONIQUE

imvisf = 0

! --- Gradient calculation
!     = 0 Standard
!     = 1 Weighted (shoulf be used with imvisf = 1)

do ii = 1, nvarmx
  iwgrec(ii) = 0
enddo

! --- Type des CL, tables de tri
!       Sera calcule apres cs_user_boundary_conditions.

do ii = 1, ntypmx
  idebty(ii) = 0
  ifinty(ii) = 0
enddo

! --- Traitement de la temperature pour couplage SYRTHES

do iscal = 1, nscamx
  icpsyr(iscal) = -999
enddo


! --- Estimateurs d'erreur pour Navier-Stokes
!       En attendant un retour d'experience et pour l'avoir,
!       on active les estimateurs par defaut.

!     Le numero d'estimateur IEST prend les valeurs suivantes
!        IESPRE : prediction
!                 L'estimateur est base sur la grandeur
!                 I = rho_n (u*-u_n)/dt + rho_n u_n grad u*
!                   - rho_n div (mu+mu_t)_n grad u* + grad P_n
!                   - reste du smb(u_n, P_n, autres variables_n)
!                 Idealement nul quand les methodes de reconstruction
!                   sont parfaites et le systeme est resolu exactement
!        IESDER : derive
!                 L'estimateur est base sur la grandeur
!                 I = div (flux de masse corrige apres etape pression)
!                 Idealement nul quand l'equation de Poisson est resolue
!                   exactement
!        IESCOR : correction
!                 L'estimateur est base sur la grandeur
!                 I = div (rho_n u_(n+1))
!                 Idealement nul quand IESDER est nul et que le passage
!                   des flux de masse aux faces vers les vitesses au centre
!                   se fait dans un espace de fonctions a divergence nulle.
!        IESTOT : total
!                 L'estimateur est base sur la grandeur
!                 I = rho_n (u_(n+1)-u_n)/dt + rho_n u_(n+1) grad u_(n+1)
!                   - rho_n div (mu+mu_t)_n grad u_(n+1) + gradP_(n_+1)
!                   - reste du smb(u_(n+1), P_(n+1), autres variables_n)
!                 Le flux du terme convectif est calcule a partir de u_(n+1)
!                   pris au centre des cellules (et non pas a partir du flux
!                   de masse aux faces actualise)

!     On evalue l'estimateur IEST selon les valeurs de IESCAL

!        iescal(iest) = 0 : l'estimateur IEST n'est pas calcule
!        iescal(iest) = 1 : l'estimateur IEST   est     calcule,
!                         sans contribution du volume  (on prend abs(I))
!        iescal(iest) = 2 : l'estimateur IEST   est     calcule,
!                         avec contribution du volume ("norme L2")
!                         soit abs(I)*SQRT(Volume_cellule),
!                         sauf pour IESCOR : on calcule abs(I)*Volume_cellule
!                         pour mesurer l'ecart en kg/s


do iest = 1, nestmx
  iescal(iest) = 0
enddo

! --- Somme de NCEPDC (pour les calculs paralleles)

ncpdct = 0

! --- Somme de NCETSM (pour les calculs paralleles)

nctsmt = 0

! --- Somme de nfbpcd (pour les calculs paralleles)

nftcdt= 0

! --- Somme de NFPT1D (pour les calculs paralleles)

nfpt1t = 0

!     Non utilisateur

! --- Calcul de la distance a la paroi
!     Seules variables utilisateur : ICDPAR, IWARNY

ineedy = 0
imajdy = 0
icdpar = -999
nitmay = 10000
nswrsy = 1
nswrgy = 100
imligy = -999
ircfly = 1
ischcy = 1
isstpy = 0
iwarny = -999
ntcmxy = 1000


blency = 0.0d0
epsily = 1.0d-8
epsrsy = 1.0d-5
epsrgy = 1.0d-5
climgy = 1.5d0
extray = 0.0d0
coumxy = 5000.d0
epscvy = 1.0d-8
yplmxy = 200.d0

! --- Methode des vortex
ivrtex = 0

! --- Calcul des efforts aux parois
ineedf = 0

! --- Ici tout optcal.f90 est initialise

!===============================================================================
! 7. TABLEAUX DE cstphy.f90
!===============================================================================

! --- Gravite

gx = 0.d0
gy = 0.d0
gz = 0.d0

! --- Vecteur rotation

icorio = 0

! --- Constantes physiques de chaque phase
!       RO0,VISCL0 et CP0 devront etre initialises par l'utilisateur
!       P0 est donne par phase, mais seul P0(1) est utilise, idem
!        pour PRED0 et XYZREF
!       T0 ne sert a rien, sauf a etre dispo pour l'utilisateur

!       IROVAR : en attendant mieux, IROVAR indique si rho est constant
!         (et c'est donc RO0) : utilise uniquement pour les suites de
!         calcul, il indique qu'il ne faut pas relire la valeur de rho
!         qui a ete stockee dans le fichier suite. Ceci permet de ne pas
!         ecraser la valeur de rho0 (eventuellement modifiee par
!         l'utilisateur) par la valeur de l'ancien calcul.
!         Le pb ne se pose pas pour la viscosite car elle n'est relue
!         que si elle est extrapolee et elle est alors forcement variable
!         (du moins, on l'espere...) : on adopte la meme methode pour la
!         symetrie.

irovar = 0
ivivar = 0
ro0    = 1.17862d0
viscl0 = 1.83337d-5
p0     = 1.01325d5
! Reset pther to p0
pther = p0
pred0  = 0.d0
xyzp0(1)= -rinfin
xyzp0(2)= -rinfin
xyzp0(3)= -rinfin
ixyzp0 = -1
t0 = 20.d0 + 273.15d0
cp0    = 1017.24d0

! --- Turbulence
!     YPLULI est mis a -GRAND*10. Si l'utilisateur ne l'a pas specifie dans usipsu, on
!     modifie sa valeur dans modini (10.88 avec les lois de paroi invariantes,
!     1/kappa sinon)
ypluli = -grand*10.d0
xkappa  = 0.42d0
cstlog  = 5.2d0

apow    = 8.3d0
bpow    = 1.d0/7.d0
cpow    = apow**(2.d0/(1.d0-bpow))
dpow    = 1.d0/(1.d0+bpow)

cmu     = 0.09d0
cmu025  = cmu**0.25d0

!   pour le k-epsilon
ce1     = 1.44d0
ce2     = 1.92d0
ce4     = 1.20d0
sigmak  = 1.00d0
! sigmae is set to 1.30d0 in modini

!   pour le Rij-epsilon standard (et SSG pour CRIJ3)
crij1  = 1.80d0
crij2  = 0.60d0
crij3  = 0.55d0
! csrij is set to 0.22d0 in modini
! sigmae is fixed to csrij/0.18d0 in modini
crijp1 = 0.50d0
crijp2 = 0.30d0

!   pour le Rij-epsilon SSG
cssgs1  = 1.70d0
cssgs2  =-1.05d0
cssgr1  = 0.90d0
cssgr2  = 0.80d0
cssgr3  = 0.65d0
cssgr4  = 0.625d0
cssgr5  = 0.20d0
cssge2  = 1.83d0

!   Rij EB-RSM
cebms1  = 1.70d0
cebms2  = 0.d0
cebmr1  = 0.90d0
cebmr2  = 0.80d0
cebmr3  = 0.65d0
cebmr4  = 0.625d0
cebmr5  = 0.20d0
! cebmr6  is used in the buoyant term
cebmr6  = 0.6d0
! csrij is set to 0.21d0 in modini
cebme2  = 1.83d0
cebmmu  = 0.22d0
xcl     = 0.122d0
! sigmae is fixed to 1.15d0 in modini
xa1     = 0.1d0
xceta   = 80.d0
xct     = 6.d0

!   pour la LES
xlesfl = 2.d0
ales   = 1.d0
bles   = 1.d0/3.d0
csmago = 0.065d0
cwale  = 0.325d0/(xlesfl*(2.d0**0.25d0))!0.25d0 ! constant defined in cs_turbulence_model.c
xlesfd = 1.5d0
smagmx = 10.d0*csmago
cdries = 26.d0

!   pour le v2f phi-model
cv2fa1 = 0.05d0
cv2fe2 = 1.85d0
cv2fmu = 0.22d0
cv2fc1 = 1.4d0
cv2fc2 = 0.3d0
cv2fct = 6.d0
cv2fcl = 0.25d0
cv2fet = 110.d0

!   pour le v2f BL-v2/k (previously known as v2-f phi-alpha, hence the cpal**)
cpale1 = 1.44d0
cpale2 = 1.83d0
cpale3 = 2.3d0
cpale4 = 0.4d0
cpalse = 1.5d0
cpalmu = 0.22d0
cpalct = 4.d0
cpalcl = 0.164d0
cpalet = 75.d0
cpalc1 = 1.7d0
cpalc2 = 0.9d0

!   pour le modele k-omega sst
ckwsk1 = 1.d0/0.85d0
ckwsk2 = 1.d0
ckwsw1 = 2.d0
ckwsw2 = 1.d0/0.856d0
ckwbt1 = 0.075d0
ckwbt2 = 0.0828d0
ckwgm1 = ckwbt1/cmu - xkappa**2/(ckwsw1*sqrt(cmu))
ckwgm2 = ckwbt2/cmu - xkappa**2/(ckwsw2*sqrt(cmu))
ckwa1  = 0.31d0
ckwc1  = 10.d0

!   pour le modele de Spalart Allmaras
csab1    = 0.1355d0
csab2    = 0.622d0
csav1    = 7.1d0
csasig   = 2.d0/3.d0
csaw1    = csab1/xkappa**2 + 1.d0/csasig*(1.d0 + csab2)
csaw2    = 0.3d0
csaw3    = 2.d0

! for the Spalart-Shur rotation/curvature correction
cssr1 = 1.d0
cssr2 = 12.d0
cssr3 = 1.d0

! for the Cazalbou rotation/curvature correction
ccaze2 = 1.83d0
ccaza  = 4.3d0
ccazsc = 0.119d0
ccazb  = 5.130d0
ccazc  = 0.453d0
ccazd  = 0.682d0

!   echelle de longueur negative, recalculee par la suite
!    ou entree par l'utilisateur
almax   = -999.d0

!   vitesse de reference pour l'initialisation de la turbulence
!    doit etre entree par l'utilisateur, sauf s'il initialise lui-meme
!    la turbulence.
uref    = -grand*10.d0

!   longueur caracteristique pour le modele de longueur de melange
!    doit etre entree par l'utilisateur
  xlomlg    = -grand*10.d0

! --- Scalaires
!       L'utilisateur devra remplir VISLS0
!       On remplira plus tard, selon les modifs utilisateur,
!         ISCSTH
!       On donne la valeur par defaut pour les autres
!       En particulier, on suppose qu'on n'a pas de variance
!         (field first_moment_id < 0)
!         qu'on clippe les variances a zero seulement,
!         qu'on ne clippe pas les scalaires (sauf a +/-GRAND)

do iscal = 1, nscamx
  iscacp(iscal) =-10
  iclvfl(iscal) = -1
  visls0(iscal) =-grand*10.d0
  sigmas(iscal) = 1.0d0
  rvarfl(iscal) = 0.8d0
enddo

! --- Turbulent flux for a scalar (Default: SGDH)
do iscal = 1, nscamx
  iturt(iscal) = 0
enddo

! For the turbulent fluxes of the scalar
c1trit = 4.15d0
c2trit = 0.55d0
c3trit = 0.5d0
c4trit = 0.d0

! For the AFM model (Algebraic flux model)
xiafm  = 0.7d0
etaafm = 0.4d0
cthafm = 0.236d0

! For the DFM (tranport equation on the turbulent flux)
cthdfm = 0.31d0

! --- Ici tout cstphy a ete initialise

!===============================================================================
! 9. INITIALISATION DES PARAMETRES DE IHM de ihmpre.f90
!===============================================================================

!     Par defaut, pas de fichier IHM consulte (on regarde ensuite si on
!       en trouve un a partir de la ligne de commande)

iihmpr = 0

!===============================================================================
! 10. INITIALISATION DES PARAMETRES ALE de albase.f90 et alstru.f90
!===============================================================================

! --- Methode ALE
iale = 0

! --- Iterations d'initialisation fluide seul
nalinf = 0

! --- Type de viscosite de maillage (isotrope par defaut)
iortvm = 0

! --- Nombre de structures internes
!     (sera relu ou recalcule)
nbstru = -999

! --- Nombre de structures externes
!     (sera relu ou recalcule)
nbaste = -999

! --- Numero d'iteration de couplage externe
ntcast = 0

! --- Parametres du couplage implicite
nalimx = 1
epalim = 1.d-5

! --- Iteration d'initialisation de l'ALE
italin = -999

! --- Tableaux des structures
do istr = 1, nstrmx
  dtstr(istr) = dtref
  do ii = 1, 3
    xstr(ii,istr)   = 0.d0
    xpstr(ii,istr)  = 0.d0
    xppstr(ii,istr) = 0.d0
    xsta(ii,istr)   = 0.d0
    xpsta(ii,istr)  = 0.d0
    xppsta(ii,istr) = 0.d0
    xstp(ii,istr)   = 0.d0
    forstr(ii,istr) = 0.d0
    forsta(ii,istr) = 0.d0
    forstp(ii,istr) = 0.d0
    xstreq(ii,istr) = 0.d0
    do jj = 1, 3
      xmstru(ii,jj,istr) = 0.d0
      xcstru(ii,jj,istr) = 0.d0
      xkstru(ii,jj,istr) = 0.d0
    enddo
  enddo
enddo

! --- Schema de couplage des structures
aexxst = -grand
bexxst = -grand
cfopre = -grand

! --- Methode de Newmark HHT
alpnmk = 0.d0
betnmk = -grand
gamnmk = -grand

!===============================================================================
! 11. INITIALISATION DES PARAMETRES DE COUPLAGE CS/CS
!===============================================================================

! --- Nombre de couplage
nbrcpl = 0

! --- Couplage uniquement par les faces
ifaccp = 0

!===============================================================================
! 12. Lagrangian arrays
!===============================================================================

statis => rvoid2
nullify(stativ)
parbor => rvoid2
tslagr => rvoid2
nullify(dlgeo)

!===============================================================================
! 13. Cavitation module
!===============================================================================

call init_cavitation
!===================

!===============================================================================
! 14. Exit
!===============================================================================

return
end subroutine
