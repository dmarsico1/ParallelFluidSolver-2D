module init

use grid
use mpi_interface

contains

subroutine init_grids

implicit none

integer :: i,j

do i=1,N_vert
	zt(i)=hz*(dble(i)-0.5)
	zm(i)=hz*dble(i)
enddo

if(myid > 0)then
	do i=1,x_dim
		xm(i)=(N_x-(num_proc-1)*(N_x/num_proc))*hx+x_dim*(myid-1)*hx+hx*(dble(i)-1.0)
		xt(i)=(N_x-(num_proc-1)*(N_x/num_proc))*hx+x_dim*(myid-1)*hx+hx*(dble(i)-0.5)
	enddo
else
	do i=1,x_dim
		xm(i)=hx*(dble(i)-1.0)
		xt(i)=hx*(dble(i)-0.5)
	enddo
endif

end subroutine init_grids

subroutine init_conds

implicit none

integer :: i,j

do i=1,N_vert
	do j=1,x_dim
		var(2+i,2+j,1)=0.0
		var(2+i,2+j,2)=0.0
		var(2+i,2+j,3)=3.0*exp((-1.0/(2.0*450.0**2))*((xt(j)- (N_x*hx/2.0))**2 + (zt(i)- (N_vert*hz/2.0))**2))
		var(2+i,2+j,4)=0.0!0.9*qvs(k)%p_var(z_dim_vec_ext(k)-z_dim_vec(k)+i)
		var(2+i,2+j,5)=0.0
	enddo
enddo

end subroutine init_conds

subroutine calc_write_dim_freq_vec

implicit none

integer :: dim_freq,i

dim_freq=dint(Tend/write_freq)

allocate(write_freq_vec(dim_freq))

do i=1,dim_freq
	write_freq_vec(i)=(i-1.0)*write_freq
enddo

end subroutine calc_write_dim_freq_vec


end module init
