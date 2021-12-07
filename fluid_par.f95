program fluid_par

use grid
use mpi_interface
use netcdf_write
use therm
use init
use step

implicit none

integer :: i,j
real*8 :: start,finish

call mpi_begin

call define_name_vars

if(myid==0)then
	call init_netcdf_vars
endif

call init_grids

call sponge_grid(N_vert,x_dim)

call calc_write_dim_freq_vec

call init_qvs

call calc_dtheta_e_vec

call init_conds

call neighbors

call get_ranks


if(myid==0)then
	allocate(write_array(N_vert,N_x,6))
endif

call cpu_time(start)

call time_stepper

call cpu_time(finish)

if(myid==0)then
	write(*,*) finish-start
endif


!call get_full_array(var(3:N_vert+2,3:x_dim+2,:))

!call write_netcdf_vars(write_array(:,:,1),write_array(:,:,2),write_array(:,:,3),&
!									write_array(:,:,4),write_array(:,:,5),write_array(:,:,6),1)

if(myid==0)then
	deallocate(write_array)
endif

if(myid==0)then
	call end_netcdf
endif

call mpi_end


end program fluid_par
