#include "scan_geo.f90"

#include "split_tensor6.f90"

program import_geo_split_tensor

  implicit none

  ! Variables for command-line arguments
  integer :: narg, iarg
  character*80 :: argument

  ! Number of elements hexa8 in the mesh
  integer :: nelt

  ! Check command-line arguments
  narg=command_argument_count()
  if (narg.le.1) then
    print *,'Invalid command-line arguments, abort program'
    stop -1
  else
    print *,''
    print *,'Program designed to split a symetric tensor written in the format Ensight Gold'
    print *,''
    print *,'Step 1 : Import .geo file'
    print *,''
  endif

  iarg = 1
  call get_command_argument(iarg,argument)
  print *,'Reading file : ',argument
  call scan_geo(argument, nelt)
  print *,'Reading complete, nelt(hexa8) = ',nelt

  print *,''
  print *,'Step 2 : splitting tensors'
  print *,''

  do iarg=2,narg
    call get_command_argument(iarg,argument)
    print *,'Splitting file : ',argument
    call split_tensor6(argument, nelt)
    print *,'Completed'
  enddo

end program import_geo_split_tensor
