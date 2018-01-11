#include "get_np.f90"

subroutine compute_xdx(geonn, geone, npcell, cgeo, xgeo, xcell, dxcell)

  implicit none

  integer, intent(in) :: geonn, geone, npcell
  integer, intent(in), dimension(npcell,geone) :: cgeo
  real*4, intent(in), dimension(geonn) :: xgeo
  real*4, intent(out), dimension(geone) :: xcell
  real*4, intent(out), dimension(geone) :: dxcell

  integer :: ii, jj
  real*4 xmin, xmax

  do ii=1,geone
    xcell(ii) = 0.
    xmax = xgeo(cgeo(1,ii))
    xmin = xmax
    do jj=1,npcell
      xcell(ii) = xcell(ii) + xgeo(cgeo(jj,ii))
      xmax = max( xmax,xgeo(cgeo(jj,ii)) )
      xmin = min( xmin,xgeo(cgeo(jj,ii)) )
    enddo
    dxcell(ii) = xmax - xmin
  enddo
  xcell = xcell / real(npcell)

end subroutine compute_xdx

subroutine import_geo(filename, geone)

  implicit none

  character*160, intent(in) :: filename
  integer, intent(out) :: geone

  integer :: fileunit, ierr

  character*80 :: buffer1
  character*5 :: extension
  integer :: ii, jj, ipart, igeo, geonn, geonnb, geoneb, np

  ! Coordinates for parts
  real*4, allocatable, dimension(:) :: xgeo, ygeo, zgeo
  real*4, allocatable, dimension(:) :: xgeob, ygeob, zgeob
  ! Connectivity for parts
  integer, allocatable, dimension(:,:) :: cgeo
  ! 1D grid
  real*4, allocatable, dimension(:) :: y1d
  ! 3D->1D mapping
  integer, allocatable, dimension(:) :: grid_3d1d

  ! Temporary user variables
  logical :: is_here, ipart2
  integer :: npcell, nx, ny, nz
  real*4 :: yb, dy
  real*4 :: xmin, dxmin, ymin, dymin, zmin, dzmin
  real*4, allocatable, dimension(:) :: xcell, ycell, zcell
  real*4, allocatable, dimension(:) :: dxcell, dycell, dzcell
  integer, allocatable, dimension(:) :: n1d
  logical, allocatable, dimension(:) :: this_ones, remaining

  ! Open .geo file
  fileunit = 40
  open(unit=fileunit,file=trim(filename),form='UNFORMATTED', &
       access='STREAM',action='READ',status='OLD',iostat=ierr)
  if (ierr.ne.0) then
    print *, 'Error in open (read) statement, abort program',ierr
    stop
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
  allocate(xgeo(geonn))
  allocate(ygeo(geonn))
  allocate(zgeo(geonn))
  do ii=1,geonn
    read(unit=fileunit) xgeo(ii)
  enddo
  do ii=1,geonn
    read(unit=fileunit) ygeo(ii)
  enddo
  do ii=1,geonn
    read(unit=fileunit) zgeo(ii)
  enddo

  ! Read element-type, ne and connectivity
  read(unit=fileunit) buffer1
  print *, buffer1
  if (buffer1(1:5).ne.'hexa8') then
    ! ipart = 1 is not the fluid domain, we assume ipart=2 will be
    ipart2=.true.
    print *,'WARNING : Expected hexa8 for ipart = 1.'
  else
    ! ipart = 1 is the fluid domain
    ipart2=.false.
  endif

  if (ipart2) then
    if (buffer1(1:5).ne.'quad4') then
      print *,'ERROR : Expected quad4 or hexa8 for ipart = 1. Abort.'
      stop
    endif
    call get_np(buffer1, np)
    read(unit=fileunit) geoneb
    print *, geoneb
    do ii=1,geoneb*np
      read(unit=fileunit) igeo ! We do not store connectivity on boundary
    enddo
  else
    call get_np(buffer1, np)
    read(unit=fileunit) geone
    print *, geone
    allocate(cgeo(np,geone))
    do ii=1,geone
      do jj=1,np
        read(unit=fileunit) cgeo(jj,ii) ! We store connectivity in the fluid domain
      enddo
    enddo
    npcell = np
  endif

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
      print *,'You must adapt the source code.'
      print *,'Abort.'
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
    if (allocated(xgeob)) deallocate(xgeob)
    if (allocated(ygeob)) deallocate(ygeob)
    if (allocated(zgeob)) deallocate(zgeob)
    allocate(xgeob(geonnb))
    allocate(ygeob(geonnb))
    allocate(zgeob(geonnb))
    do ii=1,geonnb
      read(unit=fileunit) xgeob(ii)
    enddo
    do ii=1,geonnb
      read(unit=fileunit) ygeob(ii)
    enddo
    do ii=1,geonnb
      read(unit=fileunit) zgeob(ii)
    enddo

    ! Read element-type, ne and connectivity
    read(unit=fileunit) buffer1
    print *, buffer1
    if (ipart2.and.(buffer1(1:5).ne.'hexa8').and.(ipart.eq.2)) then
      ! The fluid domain is not full hexa8
      ! OR
      ! The fluid domain is not ipart <= 2
      print *,'We are not ready for that yet. Abort.'
      stop
    endif
    if ((buffer1(1:5).ne.'quad4').and.(ipart.eq.2)) then
      print *,'WARNING : Expected quad4 for ipart = 2. Looking for hexa8.'
      call get_np(buffer1, np)
      read(unit=fileunit) geone
      print *, geone
      allocate(cgeo(np,geone))
      do ii=1,geone
        do jj=1,np
          read(unit=fileunit) cgeo(jj,ii) ! We store connectivity in the fluid domain
        enddo
      enddo
      npcell = np
      geonn = geonnb
      deallocate(xgeo)
      deallocate(ygeo)
      deallocate(zgeo)
      allocate(xgeo(geonn))
      allocate(ygeo(geonn))
      allocate(zgeo(geonn))
      xgeo = xgeob
      ygeo = ygeob
      zgeo = zgeob
    else
      call get_np(buffer1, np)
      read(unit=fileunit) geoneb
      print *, geoneb
      do ii=1,geoneb*np
        read(unit=fileunit) igeo ! We do not store connectivity on boundary
      enddo
    endif

  enddo

  ! Compute x, y, z and dx, dy, dz for each cell
  print *,'Compute cell position and size'
  allocate(xcell(geone))
  allocate(dxcell(geone))
  call compute_xdx(geonn, geone, npcell, cgeo, xgeo, xcell, dxcell)
  allocate(ycell(geone))
  allocate(dycell(geone))
  call compute_xdx(geonn, geone, npcell, cgeo, ygeo, ycell, dycell)
  allocate(zcell(geone))
  allocate(dzcell(geone))
  call compute_xdx(geonn, geone, npcell, cgeo, zgeo, zcell, dzcell)

  ! Extract cells at x=xmin and z=zmin to determine ny
  xmin = minval(xcell)
  dxmin = minval(dxcell)
  ymin = minval(ycell)
  dymin = minval(dycell)
  zmin = minval(zcell)
  dzmin = minval(dzcell)
  nx = count( (ycell.lt.(ymin+dymin/2.)) .and. (zcell.lt.(zmin+dzmin/2.)) )
  ny = count( (xcell.lt.(xmin+dxmin/2.)) .and. (zcell.lt.(zmin+dzmin/2.)) )
  nz = count( (xcell.lt.(xmin+dxmin/2.)) .and. (ycell.lt.(ymin+dymin/2.)) )
  print *,'Number of cells : ', nx, ny, nz

  ! Compute the 1D y-grid and the 3D->1D mapping
  print *, 'Compute 1D y-grid and 3D->1D mapping'
  allocate(y1d(ny))
  allocate(n1d(ny))
  allocate(this_ones(geone))
  allocate(remaining(geone))
  allocate(grid_3d1d(geone))
  y1d = 0.
  n1d = 0
  this_ones = .false.
  remaining = .true.
  grid_3d1d = 0
  yb = ymin
  dy = dymin
  do ii=1,ny
    this_ones = (ycell.ge.(yb-0.49*dy)).and.(ycell.le.yb+dy*0.49)
    n1d(ii) = count(this_ones)
    y1d(ii) = real(sum(dble(ycell), mask=this_ones) / dble(n1d(ii)))
    if ( n1d(ii).ne.(nx*nz) ) then
      print *, 'y-grid : seems to be wrong'
    endif
    remaining = remaining .and. ( .not.(this_ones) )
    yb = minval(ycell, mask=remaining)
    dy = minval(dycell, mask=remaining)
    where (this_ones)
      grid_3d1d = ii
    end where
  enddo

  ! If file is present, delete it
  extension = '.3d1d'
  inquire(file=trim(filename)//extension,exist=is_here)
  if (is_here) then
    print *, 'Remove old 3d1d mapping for ', trim(filename)
    open(unit=fileunit,file=trim(filename)//extension, &
         status='OLD',iostat=ierr)
    if (ierr.ne.0) then
      print *, 'Error in open (delete) statement, abort program',ierr
      stop
    endif
    close(unit=fileunit, status='delete')
  endif

  ! Export 1D grid and 3D->1D mapping
  print *, 'Write 1D y-grid and 3D->1D mapping'
  open(unit=fileunit,file=trim(filename)//extension,form='UNFORMATTED', &
       access='STREAM',action='WRITE',status='NEW',iostat=ierr)
  if (ierr.ne.0) then
    print *, 'Error in open (write) statement, abort program',ierr
    stop
  endif

  ! Number of 1D / 3D nodes
  write(unit=fileunit) ny
  write(unit=fileunit) geone

  ! 1D y-grid
  do ii=1,ny
    write(unit=fileunit) n1d(ii)
  enddo
  do ii=1,ny
    write(unit=fileunit) y1d(ii)
  enddo

  ! 3D->1D mapping
  do ii=1,geone
    write(unit=fileunit) grid_3d1d(ii)
  enddo

  ! Close fileunit
  close(unit=fileunit)

end subroutine import_geo

program main

  implicit none

  ! Variables for command-line arguments
  integer :: narg, iarg
  character*160 :: argument

  ! Number of elements hexa8 in the mesh
  integer :: geone

  ! Check command-line arguments
  narg=command_argument_count()
  if (narg.eq.0) then
    print *,'Invalid command-line arguments, abort program'
    stop -1
  else
    print *,''
    print *,'Program designed to extract 1D y-grid and 3D->1D mapping'
    print *,''
  endif

  do iarg=1,narg
    call get_command_argument(iarg,argument)
    print *,'Reading .geo file : ',trim(argument)
    call import_geo(argument, geone)
    print *,'Completed'
  enddo

end program main
