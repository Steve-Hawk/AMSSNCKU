!$Id: kodiss.f90,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $

subroutine setparameters(a2,r0,phi0,sigma)
implicit none
real*8,intent(out) :: a2,r0,phi0,sigma

a2 =1.d1
r0=120.d0
sigma=8.d0
phi0=1.d0/4.d1

return

end subroutine setparameters
!===================================================================
function phi(X,Y,Z) result(gont)
implicit none

double precision,intent(in)::X
double precision,intent(in)::Y
double precision,intent(in)::Z
real*8 :: gont

real*8 ::r
real*8 :: a2,r0,phi0,sigma

  call setparameters(a2,r0,phi0,sigma)
  r=dsqrt(X*X+Y*Y+Z*Z)
! configuration 1  
  gont = phi0*dtanh((r-r0)/sigma)
! configuration 2  
!  gont = phi0*dexp(-(r-r0)**2/sigma)

return

end function phi

function dphi(X,Y,Z) result(gont)
implicit none

double precision,intent(in)::X
double precision,intent(in)::Y
double precision,intent(in)::Z
real*8 :: gont

real*8 ::r
real*8 :: a2,r0,phi0,sigma

  call setparameters(a2,r0,phi0,sigma)
  r=dsqrt(X*X+Y*Y+Z*Z)
! configuration 1  
  gont = phi0/sigma*(1-(dtanh((r-r0)/sigma))**2)
! configuration 2
!  gont = -2.d0*phi0*(r-r0)/sigma*exp(-(r-r0)**2/sigma)

return

end function dphi
!==================================================================
function potential(X,Y,Z) result(gont)
implicit none

double precision,intent(in)::X
double precision,intent(in)::Y
double precision,intent(in)::Z
real*8 :: gont

real*8 :: phi
real*8 :: PI,v

real*8 :: a2,r0,phi0,sigma

  call setparameters(a2,r0,phi0,sigma)
  PI = dacos(-1.d0)

  v = phi(X,Y,Z)

!  gont = dexp(-8.d0*dsqrt(PI/3)*v)*(1-dexp(4*dsqrt(PI/3)*v))**2/32/PI/a2
  gont = 0.d0
return

end function potential

subroutine set_rho_adm2(ex,rho,X,Y,Z)

implicit none
! argument variables
integer,intent(in)::ex
double precision,intent(in),dimension(ex)::X
double precision,intent(in),dimension(ex)::Y
double precision,intent(in),dimension(ex)::Z
double precision,intent(out),dimension(ex)::rho

integer :: i
real*8 :: dphi
real*8 :: a2,r0,phi0,sigma

  call setparameters(a2,r0,phi0,sigma)
  do i=1,ex
    rho(i) = dphi(X,Y,Z)
    rho(i) = rho(i)*rho(i)
  enddo

  return

end subroutine set_rho_adm2

subroutine set_rho_adm1(ex,rho,X,Y,Z)

implicit none
! argument variables
integer,intent(in)::ex
double precision,intent(in),dimension(ex)::X
double precision,intent(in),dimension(ex)::Y
double precision,intent(in),dimension(ex)::Z
double precision,intent(out),dimension(ex)::rho

real*8 :: potential
integer :: i

  do i=1,ex
    rho(i) = potential(X(i),Y(i),Z(i))
  enddo

  return

end subroutine set_rho_adm1

subroutine set_rho_adm(ex,rho,X,Y,Z)

implicit none
! argument variables
integer,intent(in)::ex
double precision,intent(in),dimension(ex)::X
double precision,intent(in),dimension(ex)::Y
double precision,intent(in),dimension(ex)::Z
! in psivac, out rho_adm
double precision,intent(inout),dimension(ex)::rho

double precision,dimension(ex)::rho1,rho2

  call set_rho_adm1(ex,rho1,X,Y,Z)
  call set_rho_adm2(ex,rho2,X,Y,Z)

  rho = rho**4
  rho = rho**2*rho1+rho*rho2

  return

end subroutine set_rho_adm
