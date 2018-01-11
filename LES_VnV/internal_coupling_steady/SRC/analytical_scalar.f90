subroutine soldiff &
!=====================
 ( u , num, nup)

use mesh
use paramx
use cstphy
use cstnum
use albase
use cplsat
use parall
use period

!===============================================================================

implicit none

! Arguments

double precision u(ncel)
double precision num, nup

! local variables

integer          iel
double precision x

!-------------------------------------------------------------------------------

do iel = 1,ncel
   x = xyzcen(1,iel)
   if(x.le.0d0) then
      u(iel) = (2*nup*x+nup-num)/(num+nup)
   else
      u(iel) = (2*num*x+nup-num)/(num+nup)
   endif
enddo

endsubroutine
