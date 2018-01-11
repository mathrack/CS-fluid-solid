subroutine solgraddiff(u)

use mesh
use paramx
use cstphy
use cstnum
use albase
use cplsat
use parall
use period

implicit none

! Arguments

double precision u(3,ncel)

! local variables

integer          iel
double precision x, y, z

!-------------------------------------------------------------------------------

do iel = 1,ncel
   x = xyzcen(1,iel)
   y = xyzcen(2,iel)
   z = xyzcen(3,iel)
   u(1,iel) = dacos(-1.d0) * cos(dacos(-1.d0) * x ) &
                           * sin(dacos(-1.d0) * y )       &
                           * sin(dacos(-1.d0) * z )
   u(2,iel) = dacos(-1.d0) * sin(dacos(-1.d0) * x ) &
                           * cos(dacos(-1.d0) * y )       &
                           * sin(dacos(-1.d0) * z )
   u(3,iel) = dacos(-1.d0) * sin(dacos(-1.d0) * x ) &
                           * sin(dacos(-1.d0) * y )       &
                           * cos(dacos(-1.d0) * z )

enddo

end subroutine
