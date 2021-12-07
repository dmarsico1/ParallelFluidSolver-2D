module matrix

use grid
use mpi_interface

contains

function SumMat(M)

implicit none

real*8 :: M(:,:)
integer :: numrows,numcols,i,j
real*8 :: sum
real*8 :: SumMat
	
numrows=size(M,1)
numcols=size(M,2)
sum=0

do i=1,numrows
	do j=1,numcols
		sum=sum+M(i,j)
	enddo
enddo

SumMat=sum

end function SumMat

subroutine MatMult(p,q,M,N,hx,hz)

implicit none
real*8, dimension(M+2,N+2), intent(in) :: q !this is what gets INPUT.  Layer of ghost cells around q.
real*8, dimension(M,N), intent(out) :: p !this is what gets OUTPUT.  No ghost cells around p.
real*8, intent(in) :: hx,hz
integer, intent(in) :: M,N
integer :: i,j

do j=2,N+1
	do i=2,M+1
		p(i-1,j-1)=(1.0/hx**2)*(q(i,j+1)-2.0*q(i,j)+q(i,j-1))+(1.0/hz**2)*(q(i+1,j)-2.0*q(i,j)+q(i-1,j))
	enddo
enddo


end subroutine MatMult

subroutine ConjGrad(x,f,tol,N,M,hx,hz,dt)

implicit none
real*8, dimension(M+2,N+2), intent(out) :: x !This is what gets output
real*8, dimension(M,N), intent(in) :: f !This is the right hand side
real*8, intent(in) :: tol
integer :: i,j,k
real*8, dimension(M,N) :: r,q,Ax
integer :: N,M
real*8 :: hx,hz,dt
integer :: req


r=f

rhs_sum=sum(r)

call mpi_allreduce(rhs_sum,rhs_sum_tot,1,mpi_real8,mpi_sum,comm,ierr)

r=r-(rhs_sum_tot/(N_x*N_vert))

p(2:M+1,2:N+1)=r
rold=sum(r*r)

call mpi_allreduce(rold,rold_tot,1,mpi_real8,mpi_sum,comm,ierr)


do k=1,N_x**2
	call conj_exchange(p,k)
	p(1,:)=p(2,:)
	p(M+2,:)=p(M+1,:)
	call matmult(q,p,M,N,hx,hz)
	gam=sum(q*p(2:M+1,2:N+1))
	call mpi_allreduce(gam,gam_tot,1,mpi_real8,mpi_sum,comm,ierr)
	alpha_tot=rold_tot/gam_tot
	x(2:M+1,2:N+1) = x(2:M+1,2:N+1)+alpha_tot*p(2:M+1,2:N+1)
	r=r-alpha_tot*q
	rnew=sum(r*r)
	call mpi_allreduce(rnew,rnew_tot,1,mpi_real8,mpi_sum,comm,ierr)
	if(sqrt(rnew_tot)<tol)then
		exit
	endif
	p(2:M+1,2:N+1) = r + (rnew_tot/rold_tot)*p(2:M+1,2:N+1)
	rold_tot=rnew_tot
enddo

	
end subroutine ConjGrad

end module matrix
