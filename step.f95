module step

use grid
use mpi_interface
use netcdf_write
use therm
use rk_time_step

contains

subroutine time_stepper

implicit none

real*8 :: time=0.0
integer :: m=1
integer :: i,j


do while(time<Tend)
	if((write_freq_vec(m)-write_freq_dt <= time) .and. (time <= write_freq_vec(m)+write_freq_dt))then
		call calc_ql(var(3:N_vert+2,3:x_dim+2,6),var(3:N_vert+2,3:x_dim+2,4),qvs,N_vert,x_dim)
		call get_full_array(var(3:N_vert+2,3:x_dim+2,:))
		call write_netcdf_vars(write_array(:,:,1),write_array(:,:,2),write_array(:,:,3),&
									write_array(:,:,4),write_array(:,:,5),write_array(:,:,6),m)
		m=m+1
	endif
	write(*,*) "Time =", time,myid
	call var_exchange(var)

	var(1,:,1)=var(4,:,1)
	var(2,:,1)=var(3,:,1)
	var(1,:,3)=var(4,:,3)
	var(2,:,3)=var(3,:,3)
	var(1,:,4)=var(4,:,4)
	var(2,:,4)=var(3,:,4)
	var(1,:,6)=var(4,:,6)
	var(2,:,6)=var(3,:,6)
	var(1,:,2)=-var(4,:,2)
	var(2,:,2)=0.0

	var(N_vert+3,:,1)=var(N_vert+2,:,1)
	var(N_vert+4,:,1)=var(N_vert+1,:,1)
	var(N_vert+3,:,3)=var(N_vert+2,:,3)
	var(N_vert+4,:,3)=var(N_vert+1,:,3)
	var(N_vert+3,:,4)=var(N_vert+2,:,4)
	var(N_vert+4,:,4)=var(N_vert+1,:,4)
	var(N_vert+3,:,6)=var(N_vert+2,:,6)
	var(N_vert+4,:,6)=var(N_vert+1,:,6)
	var(N_vert+3,:,2)=0.0
	var(N_vert+4,:,2)=0.0

	call rk_stepper(var(:,:,1),var(1:N_vert+3,:,2),var(:,:,3),var(:,:,4),var(2:N_vert+3,2:x_dim+3,5),N_vert,x_dim,dt,hx,hz)

	time=time+dt
enddo


end subroutine time_stepper

end module step
