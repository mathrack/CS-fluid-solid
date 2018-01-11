subroutine split_tensor6(filename, nelt)

  implicit none

  character*80, intent(in) :: filename
  integer, intent(in) :: nelt

  integer :: fileunit, ierr
  character*80 :: buffer1, buffer2, buffer3
  character*4 :: extension
  integer :: ipart, ii, myint
  real*4 :: myreal4

  ! Open tensor file to split
  fileunit = 40
  open(unit=fileunit,file=filename,form='UNFORMATTED', &
       access='STREAM',action='READ',status='OLD',iostat=ierr)
  if (ierr.ne.0) then
    print *, 'Error in open (read) statement, abort program',ierr
    stop -2
  endif

  ! Read header of tensor file
  read(unit=fileunit) buffer1
  print *, buffer1
  read(unit=fileunit) buffer2 !print *, buffer2
  read(unit=fileunit) ipart !print *, ipart
  read(unit=fileunit) buffer3 !print *, buffer3

  ! Process all 6 components of the tensor
  do ii=1,6

    print *,ii,' / 6'

    ! Select extension for the name of the splitted file
    select case (ii)
      case (1)
        extension='.rxx'
      case (2)
        extension='.ryy'
      case (3)
        extension='.rzz'
      case (4)
        extension='.rxy'
      case (5)
        extension='.rxz'
      case (6)
        extension='.ryz'
    end select

    ! Open splitted file to write
    open(unit=fileunit+ii,file=trim(filename)//extension,form='UNFORMATTED', &
         access='STREAM',action='WRITE',status='NEW',iostat=ierr)
    if (ierr.ne.0) then
      print *, 'Error in open (write) statement, abort program',ierr
      stop -3
    endif

    ! Write splitted file while reading tensor file
    write(unit=fileunit+ii) buffer1
    write(unit=fileunit+ii) buffer2
    write(unit=fileunit+ii) ipart
    write(unit=fileunit+ii) buffer3
    do myint=1,nelt
      read(unit=fileunit) myreal4
      write(unit=fileunit+ii) myreal4
    enddo

    ! Write complete, close splitted file
    close(unit=fileunit+ii)

  enddo

  ! Test that EOF is correctly reached
  ! and close tensor file
  read(unit=fileunit,iostat=ierr) myreal4
  close(unit=fileunit)
  if (ierr.lt.0) then
    print *,'EOF correctly reached'
  else
    print *,'EOF not yet reached.'
    print *,'.geo file wrongly handled.'
    print *,'You must adapt the source code'
    print *,'Abort'
    stop
  endif

end subroutine split_tensor6
