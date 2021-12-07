module mpi_interface

use grid
use mpi

integer :: ierr,total_sum,comm,myid,num_proc
real*8 :: rold,rnew,rnew_tot,rold_tot,alpha_tot,gam,gam_tot,rhs_sum,rhs_sum_tot
integer :: status(mpi_status_size)
integer :: left,right
integer, dimension(:), allocatable :: counts,disp
integer, dimension(:), allocatable :: reqs
real*8, dimension(:,:), allocatable :: p

contains 

subroutine mpi_begin

comm=mpi_comm_world

call mpi_init(ierr)

call mpi_comm_rank(comm,myid,ierr)

call mpi_comm_size(comm,num_proc,ierr)

if(myid==0)then
	x_dim=N_x-(num_proc-1)*(N_x/num_proc)
else
	x_dim=N_x/num_proc
endif

allocate(var(2*Ngz+N_vert,2*NGx+x_dim,6),qvs(N_vert),ql(x_dim,N_vert),qt_vec(N_vert),dqt_vec(N_vert),&
		dtheta_vec(N_vert),dtheta_e_vec(N_Vert),xm(x_dim),xt(x_dim),zm(N_vert),zt(N_vert),spge(N_vert,x_dim),sf(N_vert,x_dim))

allocate(p(2+N_vert,2+x_dim))

allocate(reqs(4*num_proc))

end subroutine mpi_begin


subroutine neighbors

implicit none

!DAVE ADDS this thing here to fix the case with one processor

if(num_proc==1)then
	left=0
	right=0
else
	if(myid==0)then
		left=num_proc-1
		right=myid+1
	elseif(myid==num_proc-1)then
		left=myid-1
		right=0
	else
		left=myid-1
		right=myid+1
	endif
endif

end subroutine neighbors


subroutine var_exchange(var)

implicit none
real*8, dimension(N_vert+4,x_dim+4,6) :: var

integer :: i

do i=1,6
	call mpi_sendrecv(var(3,3,i),N_vert,mpi_real8,left,0,var(3,x_dim+3,i),N_vert,mpi_real8,right,0,comm,status,ierr)
	call mpi_sendrecv(var(3,x_dim+2,i),N_vert,mpi_real8,right,1,var(3,2,i),N_vert,mpi_real8,left,1,comm,status,ierr)
enddo

do i=1,6
	call mpi_sendrecv(var(3,4,i),N_vert,mpi_real8,left,0,var(3,x_dim+4,i),N_vert,mpi_real8,right,0,comm,status,ierr)
	call mpi_sendrecv(var(3,x_dim+1,i),N_vert,mpi_real8,right,1,var(3,1,i),N_vert,mpi_real8,left,1,comm,status,ierr)
enddo

end subroutine var_exchange


subroutine conj_exchange(p,k)

implicit none

real*8, dimension(N_vert+2,x_dim+2) :: p
integer :: k
integer :: local_reqs(4)

call mpi_sendrecv(p(2,2),N_vert,mpi_real8,left,0,p(2,x_dim+2),N_vert,mpi_real8,right,0,comm,status,ierr)
call mpi_sendrecv(p(2,x_dim+1),N_vert,mpi_real8,right,1,p(2,1),N_vert,mpi_real8,left,1,comm,status,ierr)

!call mpi_Irecv(p(2,x_dim+2),N_vert,mpi_real8,right,1,comm,local_reqs(1),ierr)
!call mpi_Irecv(p(2,1),N_vert,mpi_real8,left,0,comm,local_reqs(2),ierr)

!call mpi_Isend(p(2,2),N_Vert,mpi_real8,left,1,comm,local_reqs(3),ierr)
!call mpi_Isend(p(2,x_dim+1),N_vert,mpi_real8,right,0,comm,local_reqs(4),ierr)

!call mpi_allgather(local_reqs(1),4,mpi_int,reqs,4,mpi_int,comm,ierr)


end subroutine conj_exchange


subroutine pressure_exchange(p,k)

implicit none

real*8, dimension(N_vert+2,x_dim+2) :: p
integer :: k

call mpi_sendrecv(p(2,2),N_vert,mpi_real8,left,0,p(2,x_dim+2),N_vert,mpi_real8,right,0,comm,status,ierr)
call mpi_sendrecv(p(2,x_dim+1),N_vert,mpi_real8,right,1,p(2,1),N_vert,mpi_real8,left,1,comm,status,ierr)

end subroutine pressure_exchange

subroutine get_ranks

implicit none

integer :: i

if(myid==0)then
	allocate(counts(num_proc),disp(num_proc))
endif

call mpi_gather(x_dim,1,mpi_int,counts,1,mpi_int,0,comm,ierr)

if(myid==0)then
	counts=counts*N_vert
	disp(1)=0
	do i=2,num_proc
		disp(i)=disp(i-1)+counts(i-1)
	enddo
endif

end subroutine get_ranks


subroutine get_full_array(var)

implicit none
real*8, dimension(N_vert,x_dim,6) :: var
integer :: i

do i=1,6
	call mpi_gatherv(var(:,:,i),x_dim*N_vert,mpi_real8,write_array(:,:,i),counts,disp,mpi_real8,0,comm,ierr)
enddo

end subroutine get_full_array

subroutine mpi_end

implicit none

call mpi_finalize(ierr)

end subroutine mpi_end


end module mpi_interface
