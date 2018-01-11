#include "get_np.f90"

subroutine scan_geo(filename, geone)

  implicit none

  character*80, intent(in) :: filename
  integer, intent(out) :: geone

  integer :: ierr
  integer :: fileunit=40

  character*80 :: buffer1
  integer :: ii, ipart, igeo, geonn, geonnb, geoneb, np
  real*4 :: xgeo, ygeo, zgeo
  real*4 :: xgeob, ygeob, zgeob

  logical :: is_here
  character*5 :: extension = '.3d1d'

  ! Look after .geo.3d1d file before .geo file
  inquire(file=trim(filename)//extension,exist=is_here)
  if (is_here) then
    print *, 'Using available .3d1d file for ', trim(filename)
    open(unit=fileunit,file=trim(filename)//extension,form='UNFORMATTED', &
         access='STREAM',action='READ',status='OLD',iostat=ierr)
    if (ierr.ne.0) then
      print *, 'Error in open statement, abort program',ierr
      stop
    endif
    read(unit=fileunit) ii
    read(unit=fileunit) geone
    close(unit=fileunit)
    return
  endif

  ! .3d1d file not found, open .geo file
  open(unit=fileunit,file=filename,form='UNFORMATTED', &
       access='STREAM',action='READ',status='OLD',iostat=ierr)
  if (ierr.ne.0) then
    print *, 'Error in open (read) statement, abort program',ierr
    stop -2
  else
    print *,''
  endif

  ! Read header
  do ii=1,6
    read(unit=fileunit) buffer1
    print *, buffer1
  enddo

  ! Read ipart
  read(unit=fileunit) ipart ! Part 1 : fluid domain
  print *, ipart
  do ii=1,2
    read(unit=fileunit) buffer1
    print *, buffer1
  enddo

  ! Read nn and X(1:nn) Y(1:nn) Z(1:nn)
  read(unit=fileunit) geonn
  print *, geonn
  do ii=1,geonn
    read(unit=fileunit) xgeo
  enddo
  do ii=1,geonn
    read(unit=fileunit) ygeo
  enddo
  do ii=1,geonn
    read(unit=fileunit) zgeo
  enddo

  ! Read element-type, ne and connectivity
  read(unit=fileunit) buffer1
  print *, buffer1
  if (buffer1(1:5).ne.'hexa8') then
    print *,'WARNING : Expected hexa8 for ipart = 1'
  endif
  call get_np(buffer1, np)
  read(unit=fileunit) geone
  print *, geone
  do ii=1,geone*np
    read(unit=fileunit) igeo
  enddo

  ! Part 2, ...
  ! Part 2 = boundary but there may be other parts...

  do
    ! Test that EOF is correctly reached
    read(unit=fileunit,iostat=ierr) buffer1
    if (ierr.lt.0) then
      print *,'EOF correctly reached'
      close(unit=fileunit)
      exit
    elseif (ierr.gt.0) then
      print *,'EOF not reached.'
      print *,'.geo file wrongly handled.'
      print *,'You must adapt the source code'
      print *,'Abort'
      close(unit=fileunit)
      stop
    endif
    print *, buffer1

    ! EOF not reached, keep reading more parts
    ! Read ipart
    read(unit=fileunit) ipart
    print *, ipart
    do ii=1,2
      read(unit=fileunit) buffer1
      print *, buffer1
    enddo

    ! Read nn and X(1:nn) Y(1:nn) Z(1:nn)
    read(unit=fileunit) geonnb
    print *, geonnb
    do ii=1,geonnb
      read(unit=fileunit) xgeob
    enddo
    do ii=1,geonnb
      read(unit=fileunit) ygeob
    enddo
    do ii=1,geonnb
      read(unit=fileunit) zgeob
    enddo

    ! Read element-type, ne and connectivity
    read(unit=fileunit) buffer1
    print *, buffer1
    if ((buffer1(1:5).ne.'quad4').and.(ipart.eq.2)) then
      print *,'WARNING : Expected quad4 for ipart = 2'
    endif
    call get_np(buffer1, np)
    read(unit=fileunit) geoneb
    print *, geoneb
    do ii=1,geoneb*np
      read(unit=fileunit) igeo
    enddo

  enddo

end subroutine scan_geo
