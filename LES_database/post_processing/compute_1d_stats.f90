subroutine get_size(argument, ny, geone)

  implicit none

  character*160, intent(in) :: argument
  integer, intent(out) :: ny, geone

  integer :: fileunit, ierr

  fileunit = 40

  ! Open .3d1d file
  open(unit=fileunit,file=trim(argument),form='UNFORMATTED', &
       access='STREAM',action='READ',status='OLD',iostat=ierr)
  if (ierr.ne.0) then
    print *, 'Error in open (read) statement, abort program',ierr
    stop
  else
    print *,''
  endif

  ! Read ny
  read(unit=fileunit) ny

  ! Read geone
  read(unit=fileunit) geone

  ! Close .3d1d file
  close(unit=fileunit)

end subroutine get_size

subroutine get_mapping(filename, ny, geone, n1d, y1d, grid_3d1d)

  implicit none

  character*160, intent(in) :: filename
  integer, intent(in) :: ny, geone
  integer, intent(out) :: n1d(ny)
  real*4, intent(out) :: y1d(ny)
  integer, intent(out) :: grid_3d1d(geone)

  integer :: fileunit, ierr, ii
  character*1 :: buffer1

  fileunit = 40

  ! Open .3d1d file
  open(unit=fileunit,file=trim(filename),form='UNFORMATTED', &
       access='STREAM',action='READ',status='OLD',iostat=ierr)
  if (ierr.ne.0) then
    print *, 'Error in open (read) statement, abort program',ierr
    stop
  else
    print *,''
  endif

  ! Read ny
  read(unit=fileunit) ii
  ! Read geone
  read(unit=fileunit) ii

  ! 1D y-grid
  do ii=1,ny
    read(unit=fileunit) n1d(ii)
  enddo
  do ii=1,ny
    read(unit=fileunit) y1d(ii)
  enddo

  ! 3D->1D mapping
  do ii=1,geone
    read(unit=fileunit) grid_3d1d(ii)
  enddo

  ! Test that EOF is correctly reached
  read(unit=fileunit,iostat=ierr) buffer1
  if (ierr.lt.0) then
    print *,'EOF correctly reached'
    close(unit=fileunit)
  elseif (ierr.ge.0) then
    print *,'EOF not reached.'
    print *,'.3d1d file wrongly handled.'
    print *,'You must adapt the source code.'
    print *,'Abort'
    close(unit=fileunit)
    stop
  endif

end subroutine get_mapping

subroutine stats_3d_1d(filename, ny, geone, n1d, y1d, grid_3d1d)

  implicit none

  character*160, intent(in) :: filename
  integer, intent(in) :: ny, geone
  integer, intent(in) :: n1d(ny)
  real*4, intent(in) :: y1d(ny)
  integer, intent(in) :: grid_3d1d(geone)

  logical :: is_here
  integer :: fileunit, ierr, ii
  character*1 :: buffer1
  character*80 :: buffer
  real*4 :: myreal4
  real*4, dimension(ny) :: mystat1d
  real*8, dimension(ny) :: stattmp
  character*7 extension

  ! Init
  fileunit = 40
  mystat1d = 0.0
  stattmp = 0.d0
  extension = '.1d.dat'

  ! Open 3d stats file
  open(unit=fileunit,file=trim(filename),form='UNFORMATTED', &
       access='STREAM',action='READ',status='OLD',iostat=ierr)
  if (ierr.ne.0) then
    print *, 'Error in open (read) statement, abort program',ierr
    stop
  else
    print *,''
  endif

  ! Read 3D scalar file and compute 1D stat
  read(unit=fileunit) buffer
  read(unit=fileunit) buffer
  read(unit=fileunit) ii
  read(unit=fileunit) buffer
  do ii=1,geone
    read(unit=fileunit) myreal4
    stattmp(grid_3d1d(ii)) = stattmp(grid_3d1d(ii)) + dble(myreal4)
  enddo
  mystat1d = real( stattmp / dble(n1d) )

  ! Test that EOF is correctly reached
  read(unit=fileunit,iostat=ierr) buffer1
  if (ierr.lt.0) then
    print *,'EOF correctly reached'
    close(unit=fileunit)
  elseif (ierr.ge.0) then
    print *,'EOF not reached.'
    print *,'.3d1d file wrongly handled.'
    print *,'You must adapt the source code.'
    print *,'Abort'
    close(unit=fileunit)
    stop
  endif

  ! If file is present, delete it
  inquire(file=trim(filename)//extension,exist=is_here)
  if (is_here) then
    print *, 'Remove old 1d stat file for ', trim(filename)
    open(unit=fileunit,file=trim(filename)//extension, &
         status='OLD',iostat=ierr)
    if (ierr.ne.0) then
      print *, 'Error in open (delete) statement, abort program',ierr
      stop
    endif
    close(unit=fileunit, status='delete')
  endif

  ! Open 1d stat file
  open(unit=fileunit,file=trim(filename)//extension, &
       action='WRITE',status='NEW',iostat=ierr)
  if (ierr.ne.0) then
    print *, 'Error in open (write) statement, abort program',ierr
    stop
  else
    print *,''
  endif

  ! Write 1D data
  do ii=1,ny
    write(unit=fileunit,fmt=*) y1d(ii), mystat1d(ii)
  enddo

  ! Close 1D stat file
  close(unit=fileunit)

end subroutine stats_3d_1d

program main

  implicit none

  ! Variables for command-line arguments
  integer :: narg, iarg
  character*160 :: argument

  ! Number of elements hexa8 in the mesh
  integer :: geone
  ! Number of elements in 1D grid
  integer :: ny
  ! Number of 3D nodes for each 1D point
  integer, allocatable, dimension(:) :: n1d
  ! 1D y-grid
  real*4, allocatable, dimension(:) :: y1d
  ! 3D->1D mapping
  integer, allocatable, dimension(:) :: grid_3d1d

  ! Check command-line arguments
  narg=command_argument_count()
  if (narg.le.1) then
    print *,'Invalid command-line arguments, abort program'
    stop -1
  else
    print *,''
    print *,'Program designed to compute 1D statistics'
    print *,'  based on 3D->1D mapping'
    print *,''
  endif

  ! File must be the .3d1d mapping
  iarg = 1
  call get_command_argument(iarg,argument)
  print *,'Reading .3d1d file : ',trim(argument)
  call get_size(argument, ny, geone)
  allocate( n1d(ny) )
  allocate( y1d(ny) )
  allocate( grid_3d1d(geone) )
  call get_mapping(argument, ny, geone, n1d, y1d, grid_3d1d)
  print *,'Completed'

  do iarg=2,narg
    call get_command_argument(iarg,argument)
    print *,'Reading 3D stats file : ',trim(argument)
    call stats_3d_1d(argument, ny, geone, n1d, y1d, grid_3d1d)
    print *,'Completed'
  enddo

end program main
