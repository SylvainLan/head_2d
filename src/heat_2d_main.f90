!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program heat_2d_main

  use mpi
  use mod_precision
  use mod_cartesian_grid
  use mod_utils
  use mod_heat_2d
  use mod_integration

  implicit none

  real(kind=dp) :: tini, tend
  integer :: nt
  real(kind=dp) :: xmin, xmax
  real(kind=dp) :: ymin, ymax
  integer :: nxib
  integer :: nyib
  type(cartesian_grid_type) :: grid
  real(kind=dp), allocatable, dimension(:) :: unum
  character(len=20) :: method
  real(kind=dp) :: tol
  real(kind=dp) :: norm_err
  real(kind=dp) :: loc_elapsed, elapsed
  integer :: proc_nb, nb_procs
  integer :: ierr

  ! init mpi environnement for mumps
  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world, nb_procs, ierr)
  call mpi_comm_rank(mpi_comm_world, proc_nb, ierr)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! initialisation

  call read_data(tini, tend, nt, xmin, xmax, ymin, ymax, nxib, nyib, method, tol)

  if (proc_nb == 0) then

    print *, "Resolution of 2d heat equation" 
 
    print "(a, f12.3)", "   tini = ", tini
    print "(a, f12.3)", "   tend = ", tend
    print *, "  nt   =", nt
    print "(a, f12.3)", "   xmin =", xmin
    print "(a, f12.3)", "   xmax =", xmax
    print "(a, f12.3)", "   xmin =", ymin
    print "(a, f12.3)", "   xmax =", ymax
    print *, "  nxib =", nxib
    print *, "  nyib =", nyib
    print *, "  integration method : ", method
    print "(a, es12.3)", "   tolerance (for Radau5 and Rock4) =", tol
 
    call init_cartesian_grid(xmin, xmax, ymin, ymax, nxib, nyib, grid)
  
    ! init heat 
    call init_heat(grid)
 
    ! allocation
    allocate(unum(grid%nx*grid%ny)) 
  
    ! compute and save initial solution 
    call heat_init_sol(tini, grid, unum)
    call save_sol("sol_ini.dat", grid, unum)

  endif
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! integration
  loc_elapsed = mpi_wtime()
  call integrate(method, tol, grid%nx, grid%ny, unum, nt, tini, tend)
  loc_elapsed = mpi_wtime() - loc_elapsed
  call mpi_reduce(loc_elapsed, elapsed, 1, mpi_double, mpi_max, 0, mpi_comm_world, ierr)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! print sol and compute error
  if (proc_nb == 0) then

    print *
    print *, "Time (s) to integrate :", elapsed
 
    ! write numeric solution
    call save_sol("sol_num.dat", grid, unum)
 
    ! compute error
    call compute_error(grid, unum)
   
    deallocate(unum) 
  end if
 
  ! terminate mpi environment
  call mpi_finalize(ierr)

end program heat_2d_main
