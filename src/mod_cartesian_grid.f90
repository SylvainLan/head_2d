module mod_cartesian_grid

type cartesian_grid_type
  ! minimal and maximal coordinates of domain
  real(kind=kind(0.d0)) :: xmin, xmax
  real(kind=kind(0.d0)) :: ymin, ymax
  ! number of points including boundary points
  integer :: nxib
  integer :: nyib
  ! number of interior points 
  integer :: nx
  integer :: ny
  ! step size 
  real(kind=kind(0.d0)) :: dx
  real(kind=kind(0.d0)) :: dy
end type cartesian_grid_type

contains

subroutine init_cartesian_grid(xmin, xmax, ymin, ymax, nxib, nyib, grid)
  implicit none
  real(kind=kind(0.d0)) :: xmin, xmax
  real(kind=kind(0.d0)) :: ymin, ymax
  integer :: nxib
  integer :: nyib
  type(cartesian_grid_type) :: grid

  grid%xmin = xmin
  grid%xmax = xmax
  grid%ymin = ymin
  grid%ymax = ymax


  grid%nxib = nxib
  grid%nyib = nyib

  grid%nx = nxib - 2
  grid%ny = nyib - 2

  grid%dx = (xmax-xmin) / (nxib-1)
  grid%dy = (ymax-ymin) / (nyib-1)

end subroutine init_cartesian_grid

end module mod_cartesian_grid
