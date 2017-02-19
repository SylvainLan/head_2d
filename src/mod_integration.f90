module mod_integration

use mpi
use mod_precision
use mod_forward_euler
use mod_backward_euler
use mod_radau
use mod_rock
use mod_heat_2d

contains

subroutine integrate(method, tol, nx, ny, u, nt, tini, tend)
  implicit none
  ! arguments
  character(len=20) :: method
  real(kind=dp) :: tol
  integer :: nx, ny
  real(kind=dp), dimension(:) :: u
  integer :: nt
  real(kind=dp) :: tini, tend
  ! local variables
  integer :: n 
  integer :: proc_nb
  integer :: ierr

  n = nx * ny

  call mpi_comm_rank(mpi_comm_world, proc_nb, ierr) 
 
  select case (method) 
    case ("forward_euler")
      if (proc_nb == 0) then
        print *
        print *, "Forward Euler integration"
        call forward_euler_integration(f_heat, n, u, nt, tini, tend)
      endif
    case ("backward_euler")
      if (proc_nb == 0) then
        print *
        print *, "Forward Euler integration"
      endif
      call backward_euler_integration(build_laplacian_matrix, n, u, nt, tini, tend)
    case ("radau5")
      if (proc_nb == 0) then
        print *
        print *, "Radau5 integration"
        call radau5_integration(f_heat, tol, nx, ny, u, nt, tini, tend)
      endif
    case ("rock4")
      if (proc_nb == 0) then
        print *
        print *, "Rock4 integration"
        call rock4_integration(f_heat, tol, n, u, nt, tini, tend)
      endif
    case default
      print *
      print *, method, "Unknown integration method"
      call exit(0)
  end select

end subroutine integrate

end module mod_integration
