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

!===============================================================================

implicit none

! Arguments

double precision u(ncel)

! local variables

integer          iel
double precision x, y, z

!-------------------------------------------------------------------------------

do iel = 1,ncel
   x = xyzcen(1,iel)
   y = xyzcen(2,iel)
   z = xyzcen(3,iel)
   u(iel) = 1.d0 + sin(dacos(-1.d0) * x ) &
                 * sin(dacos(-1.d0) * y ) &
                 * sin(dacos(-1.d0) * z )
enddo

endsubroutine
