module mod_utils

use mod_precision
use mod_cartesian_grid

contains

subroutine read_data(tini, tend, nt, xmin, xmax, ymin, ymax, nxib, nyib, method, tol)
  implicit none
  real(kind=dp) :: tini, tend
  integer :: nt
  real(kind=dp) :: xmin, xmax
  real(kind=dp) :: ymin, ymax
  integer :: nxib
  integer :: nyib
  character(len=20) :: method
  real(kind=dp) :: tol

  open(unit=10, file="input.dat")
  read(10, *) tini
  read(10, *) tend
  read(10, *) nt
  read(10, *) xmin
  read(10, *) xmax
  read(10, *) ymin
  read(10, *) ymax
  read(10, *) nxib
  read(10, *) nyib
  read(10, *) method
  read(10, *) tol
  method = trim(method)
  close(10)

end subroutine read_data


subroutine save_sol(file_name, grid, u)
  implicit none
  character(len=*) :: file_name
  type(cartesian_grid_type) :: grid
  real(kind=dp), dimension(:) :: u
  ! local variables
  real(kind=dp) :: xi, yi
  real(kind=dp) :: xmin, ymin
  real(kind=dp) :: dx, dy
  integer :: nx, ny
  integer :: inx, iny, irow

  nx = grid%nx
  ny = grid%ny

  dx = grid%dx
  dy = grid%dy

  xmin = grid%xmin
  ymin = grid%ymin

  open(unit=10, file=file_name)

  write(10,*) "#", nx*ny

  do iny = 1, ny
    yi = (ymin + dy)  + (iny-1)*dy
    do inx = 1, nx
      irow = inx + (iny-1)*nx
      xi = (xmin + dx)  + (inx-1)*dx
      write(10,*) xi, yi, u(irow)
    end do
    write(10,*) 
  end do

  close(10)

end subroutine save_sol

subroutine compute_error(grid, u)
  implicit none
  type(cartesian_grid_type) :: grid
  real(kind=dp), dimension(:) :: u
  real(kind=dp) :: norm_err
  ! local variables
  real(kind=dp), allocatable, dimension(:) :: uexa
  integer :: nx, ny
  real(kind=dp) :: dx, dy
  real(kind=dp) :: xi, yi
  integer :: irow
  character(len=10) :: str
  integer :: ntot, read_ntot
  integer :: ierr

  nx = grid%nx
  ny = grid%ny
  ntot = nx * ny

  dx = grid%dx
  dy = grid%dy

  ! read quasi exact solution
  open(unit=10, file="sol_ref.dat", action="read", iostat=ierr)
  if (ierr /= 0) then
    print *
    print *, "Cannot open sol_ref.dat to compute error"
    return
  else
    read(10,*) str, read_ntot
    if (ntot /= read_ntot) then
      print *, "  Bad dimension of vector read "
      print *, "  Cannot compute norm of error "
      return 
    else
      allocate(uexa(ntot)) 
      do irow = 1, ntot
        read(10,*) xi, yi, uexa(irow)
      end do

      norm_err = 0.d0
      do irow = 1, ntot
        norm_err = norm_err + ((u(irow) - uexa(irow)) * (u(irow) - uexa(irow)))
      end do 

      norm_err = sqrt((dx*dy)*norm_err)
      print *
      print *, "|| unum - uexa || = ", norm_err
    end if
  end if

  deallocate(uexa)


end subroutine compute_error

end module mod_utils
