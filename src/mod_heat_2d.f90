module mod_heat_2d

use mod_precision
use mod_cartesian_grid
use mod_csr_matrix

real(kind=dp), private, parameter :: d = 1.d0
real(kind=dp), private :: doverdxdx
real(kind=dp), private :: doverdydy
integer, private :: nx
integer, private :: ny

contains

subroutine init_heat(grid)
  implicit none
  type(cartesian_grid_type) :: grid

  doverdxdx = d / (grid%dx*grid%dx)
  doverdydy = d / (grid%dy*grid%dy)

  nx = grid%nx
  ny = grid%ny

end subroutine init_heat

function f_heat(n, u)
  implicit none
  integer :: n
  real(kind=dp), dimension(n) :: u
  real(kind=dp), dimension(n) :: f_heat
  ! local 
  real(kind=dp) :: diag_x_coef, diag_y_coef
  integer :: inx, iny, irow

  f_heat = 0.d0

  do iny = 1, ny

    if (iny==1 .or. iny==ny) then
      diag_y_coef = -1.d0
    else
      diag_y_coef = -2.d0             
    end if

    do inx = 1, nx

      if (inx==1 .or. inx==nx) then
        diag_x_coef = -1.d0
      else
        diag_x_coef = -2.d0             
      end if

      irow = inx + (iny-1)*nx

      if (iny > 1) then
        f_heat(irow) = f_heat(irow) + (doverdydy * u(irow-nx))
      end if

      if (inx > 1) then
        f_heat(irow) = f_heat(irow) + (doverdxdx * u(irow-1))
      end if

      f_heat(irow) = f_heat(irow) + ((doverdxdx * (diag_x_coef*u(irow))) + (doverdydy * (diag_y_coef*u(irow))))

      if (inx < nx) then
        f_heat(irow) = f_heat(irow) + (doverdxdx * u(irow+1))
      end if

      if (iny < ny) then
        f_heat(irow) = f_heat(irow) + (doverdydy * u(irow+nx))
      end if

    end do
  end do

end function f_heat


subroutine build_laplacian_matrix(n, a)
  implicit none
  integer :: n
  type(csr_matrix) :: a
  ! local variables
  real(kind=dp) :: diag_x_coef, diag_y_coef
  integer :: nz
  integer :: inx, iny, irow, inz


  nz = 5*(nx-2)*(ny-2) + 2*(4*(nx-2)) + 2*(4*(ny-2)) + 4*3
  call create_csr_matrix(nx*ny, nz, a)

  a%row_ptr(1) = 1
  inz = 1

  do iny = 1, ny

    if (iny==1 .or. iny==ny) then
      diag_y_coef = -1.d0
    else
      diag_y_coef = -2.d0             
    end if

    do inx = 1, nx

      if (inx==1 .or. inx==nx) then
        diag_x_coef = -1.d0
      else
        diag_x_coef = -2.d0             
      end if

      irow = inx + (iny-1)*nx

      if (iny > 1) then
        a%col_ind(inz) = irow - nx
        a%val(inz) = doverdydy
        inz = inz + 1
      end if

      if (inx > 1) then
        a%col_ind(inz) = irow - 1
        a%val(inz) = doverdxdx
        inz = inz + 1
      end if

      a%col_ind(inz) = irow
      a%val(inz) = diag_x_coef*doverdxdx + diag_y_coef*doverdydy
      inz = inz + 1

      if (inx < nx) then
        a%col_ind(inz) = irow + 1
        a%val(inz) = doverdxdx
        inz = inz + 1
      end if

      if (iny < ny) then
        a%col_ind(inz) = irow + nx
        a%val(inz) = doverdydy
        inz = inz + 1
      end if

      a%row_ptr(irow+1) = inz

    end do
  end do

end subroutine build_laplacian_matrix

subroutine heat_init_sol(t, grid, u)
  implicit none
  real(kind=dp) :: t
  type(cartesian_grid_type) :: grid
  real(kind=dp), dimension(grid%nx) :: u
  ! local
  real(kind=dp), parameter :: pi = 4.d0 * atan(1.d0)
  real(kind=dp) :: xi, yi
  real(kind=dp) :: xmin, ymin
  real(kind=dp) :: dx, dy
  integer :: inx, iny, irow

  dx = grid%dx
  dy = grid%dy
  xmin = grid%xmin
  ymin = grid%ymin

  do iny = 1, ny
    yi = (ymin + dy) + (iny-1)*dy
    do inx = 1, nx
      xi = (xmin + dx) + (inx-1)*dx
      irow = inx + (iny-1)*nx
      u(irow) = (1.d0 / (4.d0*pi*d*t)) * exp(-((xi*xi)+(yi*yi)) / (4.d0*d*t))
    end do
  end do

end subroutine heat_init_sol

end module mod_heat_2d
