module poisson

use grid
use mpi_interface
use matrix

contains

subroutine GradP(Px,Pz,P,hx,hz,M,N)

!No ghost cells around P; we adjust for periodicity at the left boundary

implicit none
real*8, dimension(M,N+1) :: P
real*8, dimension(M,N) :: Px
real*8, dimension(M-1,N) :: Pz
real*8 :: hx,hz
integer :: M,N,j

!Px(:,1)=(1.0/hx)*(P(:,1)-P(:,N))

do j=2,N+1
	Px(:,j-1)=(1.0/hx)*(P(:,j)-P(:,j-1))
enddo


do j=1,M-1
	Pz(j,:)=(1.0/hz)*(P(j+1,2:N+1)-P(j,2:N+1))
enddo


end subroutine GradP

subroutine Div(D,U,W,hx,hz,M,N)

!Layer of ghost cells at top and bottom for W
!Layer of ghost cells at right boundary for U

implicit none
real*8, dimension(M,N) :: D
real*8, dimension(M,N+1), intent(in) :: U
real*8, dimension(M+1,N), intent(in) :: W
real*8 :: hx,hz
integer :: M,N,i
real*8, dimension(M,N) :: DivX,DivZ

do i=1,N
	DivX(:,i)=(1.0/hx)*(U(:,i+1)-U(:,i))
enddo

do i=1,M
	DivZ(i,:)=(1.0/hz)*(W(i+1,:)-W(i,:))
enddo

D=DivX+DivZ

end subroutine Div

subroutine DivOfDiffusion(DoD,U,W,hx,hz,M,N)

implicit none
real*8, dimension(M,N) :: DoD
real*8, dimension(M+2,N+3) :: U
real*8, dimension(M+3,N+2) :: W
real*8, dimension(M,N) :: Uxxx,Uzzx,Wxxz,Wzzz
real*8 :: hx,hz
integer :: M,N,i,j

do i=1,M
	do j=1,N
		Uxxx(i,j)=(1.0/hx**3)*((U(i+1,j+3)-2.0*U(i+1,j+2)+U(i+1,j+1))-(U(i+1,j+2)-2.0*U(i+1,j+1)+U(i+1,j)))
	enddo
enddo

do i=1,M
	do j=1,N
		Uzzx(i,j)=(1.0/hx)*(1.0/hz**2)*((U(i,j+2)-2.0*U(i+1,j+2)+U(i+2,j+2))-(U(i,j+1)-2.0*U(i+1,j+1)+U(i+2,j+1)))
	enddo
enddo

do i=1,M
	do j=1,N
		Wxxz(i,j)=(1.0/hz)*(1.0/hx**2)*((W(i+1,j+2)-2.0*W(i+1,j+1)+W(i+1,j))-(W(i+2,j+2)-2.0*W(i+2,j+1)+W(i+2,j)))
	enddo
enddo

do i=1,M
	do j=1,N
		Wzzz(i,j)=(1.0/hz**3)*((W(i,j+1)-2.0*W(i+1,j+1)+W(i+2,j+1))-(W(i+1,j+1)-2.0*W(i+2,j+1)+W(i+3,j+1)))
	enddo
enddo

DoD=ev*(Uxxx+Uzzx+Wxxz+Wzzz)

end subroutine DivOfDiffusion

subroutine PoissonSolve(U,W,U1,W1,P,N,M,dt,tol,hx,hz)

implicit none
real*8, dimension(M+2,N+2), intent(out) :: P
real*8, dimension(M,N+1) :: U
real*8, dimension(M+1,N) :: W
real*8, dimension(M+2,N+3) :: U1
real*8, dimension(M+3,N+2) :: W1
real*8, dimension(M,N) :: DoD,D,LP
real*8, dimension(M+2,N+2) :: PWB
real*8 :: hx,hz,dt,tol
integer :: M,N,grid_num
integer :: i,j

PWB=P

!if(grid_num==1)then
!	BcT=2
!else
!	BcT=2
!endif

call Div(D,U,W,hx,hz,M,N)

D=(1.0/dt)*D


!call DivOfDiffusion(DoD,U1,W1,hx,hz,M,N)

!call MatMult(LP,P,M,N,hx,hz)

!if(grid_num>1)then
!!LP(M,:)=LP(M,:)-(1.0/hz**2)*P(M+2,:)+(1.0/hz**2)*fu
!D(M,:)=D(M,:)!-(1.0/hz**2)*P(M+2,:)
!endif


call ConjGrad(PWB,D,tol,N,M,hx,hz,dt)

P(2:M+1,2:N+1)=PWB(2:M+1,2:N+1)

end subroutine PoissonSolve


end module poisson

