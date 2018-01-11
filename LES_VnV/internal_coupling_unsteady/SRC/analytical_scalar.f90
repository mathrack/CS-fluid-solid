subroutine soldiff &
!=====================
 (u)

use mesh
use paramx
use cstphy
use cstnum
use albase
use cplsat
use parall
use period
use case_setup, only : Tm, Tp, nm, np

!===============================================================================

implicit none

! Arguments

double precision u(ncel)

! local variables

integer          iel
double precision x

!-------------------------------------------------------------------------------

do iel = 1,ncel
   x = xyzcen(1,iel)
   if(x.le..0d0) then
      u(iel) = Tm * sin(dacos(-1.d0)*x*nm)
   else
      u(iel) = Tp * sin(dacos(-1.d0)*x*np)
   endif
enddo

endsubroutine
