!$Id: getnp4.f90,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
!-----------------------------------------------------------------------------
!
! compute the Newman-Penrose Weyl scalar Psi4
! for BSSN dynamical variables
!
!-----------------------------------------------------------------------------

  subroutine getnp4(ex, X, Y, Z,                            &
               chi, trK, Sphi,&
               dxx,gxy,gxz,dyy,gyz,dzz, &
               Axx,Axy,Axz,Ayy,Ayz,Azz, &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,&
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,&
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,&
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,&
               Rpsi4, Ipsi4, &
               symmetry)

  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3),symmetry
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dxx,gxy,gxz,dyy,gyz,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: chi,trK,Sphi
! physical second kind of connection  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: Gamzyy, Gamzyz, Gamzzz
! physical Ricci tensor  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out):: Rpsi4,Ipsi4

!~~~~~~> Other variables:

  real*8, dimension(ex(1),ex(2),ex(3)) :: f,fx,fy,fz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: chix,chiy,chiz,chipn1
  real*8, dimension(ex(1),ex(2),ex(3)) :: vx,vy,vz,ux,uy,uz,wx,wy,wz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Exx,Exy,Exz,Eyy,Eyz,Ezz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Bxx,Bxy,Bxz,Byy,Byz,Bzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Axxx,Axxy,Axxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Axyx,Axyy,Axyz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Axzx,Axzy,Axzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Ayyx,Ayyy,Ayyz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Ayzx,Ayzy,Ayzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Azzx,Azzy,Azzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupyy,gupyz,gupzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: uuwwxx,uuwwxy,uuwwxz,uuwwyy,uuwwyz,uuwwzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: uwxx,uwxy,uwxz,uwyy,uwyz,uwzz

  real*8, parameter :: ZEO = 0.d0, ONE = 1.d0, TWO = 2.d0
  real*8, parameter :: F1o3 = 1.d0/3.d0, FOUR = 4.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  real*8            :: dX, dY, dZ
  integer::i,j,k
  real*8,parameter::TINYRR=1.d-14
  real*8 :: PI

  PI = dacos(-ONE)

  dX = X(2) - X(1)
  dY = Y(2) - Y(1)
  dZ = Z(2) - Z(1)

  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE
  chipn1= chi + ONE

! invert tilted metric
  gupzz =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
           gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / gupzz
  gupxy = - ( gxy * gzz - gyz * gxz ) / gupzz
  gupxz =   ( gxy * gyz - gyy * gxz ) / gupzz
  gupyy =   ( gxx * gzz - gxz * gxz ) / gupzz
  gupyz = - ( gxx * gyz - gxy * gxz ) / gupzz
  gupzz =   ( gxx * gyy - gxy * gxy ) / gupzz

! initialize U, V, W vetors	
  do i=1,ex(1)
  do j=1,ex(2)
  do k=1,ex(3)
     if(abs(X(i)) < TINYRR .and. abs(Y(j)) < TINYRR .and. abs(Z(k)) < TINYRR)then
        vx(i,j,k) = TINYRR
        vy(i,j,k) = TINYRR
        vz(i,j,k) = TINYRR
     else
        vx(i,j,k) = X(i)
        vy(i,j,k) = Y(j)
        vz(i,j,k) = Z(k)
     endif
     if(abs(X(i)) < TINYRR .and. abs(Y(j)) < TINYRR)then
        ux(i,j,k) = - TINYRR
        uy(i,j,k) = TINYRR
        uz(i,j,k) = ZEO
        wx(i,j,k) = TINYRR*Z(k)
        wy(i,j,k) = TINYRR*Z(k)
        wz(i,j,k) = -2*TINYRR*TINYRR
     else
        ux(i,j,k) = - Y(j)
        uy(i,j,k) = X(i)
        uz(i,j,k) = ZEO
        wx(i,j,k) = X(i)*Z(k)
        wy(i,j,k) = Y(j)*Z(k)
        wz(i,j,k) = -(X(i)*X(i) + Y(j)*Y(j))
     endif
  enddo
  enddo
  enddo

  f = 1.d0/chipn1

  fx = gxx*vx*vx + gyy*vy*vy + gzz*vz*vz &
     +(gxy*vx*vy + gxz*vx*vz + gyz*vy*vz)*TWO
  fx = dsqrt(fx*f)
  vx = vx/fx
  vy = vy/fx
  vz = vz/fx

  fx = gxx*vx*ux + gxy*vx*uy + gxz*vx*uz + &
       gxy*vy*ux + gyy*vy*uy + gyz*vy*uz + &
       gxz*vz*ux + gyz*vz*uy + gzz*vz*uz
  fx = fx*f
  ux = ux - fx*vx
  uy = uy - fx*vy
  uz = uz - fx*vz
  fx = gxx*ux*ux + gyy*uy*uy + gzz*uz*uz &
     +(gxy*ux*uy + gxz*ux*uz + gyz*uy*uz)*TWO
  fx = dsqrt(fx*f) 
  ux = ux/fx
  uy = uy/fx
  uz = uz/fx

  fx = gxx*vx*wx + gxy*vx*wy + gxz*vx*wz + &
       gxy*vy*wx + gyy*vy*wy + gyz*vy*wz + &
       gxz*vz*wx + gyz*vz*wy + gzz*vz*wz
  fx = fx*f       
  wx = wx - fx*vx
  wy = wy - fx*vy
  wz = wz - fx*vz
  fx = gxx*ux*wx + gxy*ux*wy + gxz*ux*wz + &
       gxy*uy*wx + gyy*uy*wy + gyz*uy*wz + &
       gxz*uz*wx + gyz*uz*wy + gzz*uz*wz
  fx = fx*f       
  wx = wx - fx*ux
  wy = wy - fx*uy
  wz = wz - fx*uz
  fx = gxx*wx*wx + gyy*wy*wy + gzz*wz*wz &
     +(gxy*wx*wy + gxz*wx*wz + gyz*wy*wz)*TWO
  fx = dsqrt(fx*f)
  wx = wx/fx
  wy = wy/fx
  wz = wz/fx

  call fderivs(ex,Axx,Axxx,Axxy,Axxz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fderivs(ex,Axy,Axyx,Axyy,Axyz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,0)
  call fderivs(ex,Axz,Axzx,Axzy,Axzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,0)
  call fderivs(ex,Ayy,Ayyx,Ayyy,Ayyz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fderivs(ex,Ayz,Ayzx,Ayzy,Ayzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,0)
  call fderivs(ex,Azz,Azzx,Azzy,Azzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)

  call fderivs(ex,chi,chix,chiy,chiz,X,Y,Z,SYM,SYM,SYM,symmetry,0)
  call fderivs(ex,trK,fx,fy,fz,X,Y,Z,SYM,SYM,SYM,symmetry,0)
! compute D_k K_ij up to chi^-1
  Axxx = Axxx - (Gamxxx*Axx + Gamyxx*Axy + Gamzxx*Axz)*TWO - chix/chipn1*Axx + F1o3*gxx*fx
  Axxy = Axxy - (Gamxxy*Axx + Gamyxy*Axy + Gamzxy*Axz)*TWO - chiy/chipn1*Axx + F1o3*gxx*fy
  Axxz = Axxz - (Gamxxz*Axx + Gamyxz*Axy + Gamzxz*Axz)*TWO - chiz/chipn1*Axx + F1o3*gxx*fz
  Ayyx = Ayyx - (Gamxxy*Axy + Gamyxy*Ayy + Gamzxy*Ayz)*TWO - chix/chipn1*Ayy + F1o3*gyy*fx
  Ayyy = Ayyy - (Gamxyy*Axy + Gamyyy*Ayy + Gamzyy*Ayz)*TWO - chiy/chipn1*Ayy + F1o3*gyy*fy
  Ayyz = Ayyz - (Gamxyz*Axy + Gamyyz*Ayy + Gamzyz*Ayz)*TWO - chiz/chipn1*Ayy + F1o3*gyy*fz
  Azzx = Azzx - (Gamxxz*Axz + Gamyxz*Ayz + Gamzxz*Azz)*TWO - chix/chipn1*Azz + F1o3*gzz*fx
  Azzy = Azzy - (Gamxyz*Axz + Gamyyz*Ayz + Gamzyz*Azz)*TWO - chiy/chipn1*Azz + F1o3*gzz*fy
  Azzz = Azzz - (Gamxzz*Axz + Gamyzz*Ayz + Gamzzz*Azz)*TWO - chiz/chipn1*Azz + F1o3*gzz*fz
  Axyx = Axyx - (Gamxxy*Axx + Gamyxy*Axy + Gamzxy*Axz + &
                 Gamxxx*Axy + Gamyxx*Ayy + Gamzxx*Ayz) - chix/chipn1*Axy + F1o3*gxy*fx
  Axyy = Axyy - (Gamxyy*Axx + Gamyyy*Axy + Gamzyy*Axz + &
                 Gamxxy*Axy + Gamyxy*Ayy + Gamzxy*Ayz) - chiy/chipn1*Axy + F1o3*gxy*fy
  Axyz = Axyz - (Gamxyz*Axx + Gamyyz*Axy + Gamzyz*Axz + &
                 Gamxxz*Axy + Gamyxz*Ayy + Gamzxz*Ayz) - chiz/chipn1*Axy + F1o3*gxy*fz
  Axzx = Axzx - (Gamxxz*Axx + Gamyxz*Axy + Gamzxz*Axz + &
                 Gamxxx*Axz + Gamyxx*Ayz + Gamzxx*Azz) - chix/chipn1*Axz + F1o3*gxz*fx
  Axzy = Axzy - (Gamxyz*Axx + Gamyyz*Axy + Gamzyz*Axz + &
                 Gamxxy*Axz + Gamyxy*Ayz + Gamzxy*Azz) - chiy/chipn1*Axz + F1o3*gxz*fy
  Axzz = Axzz - (Gamxzz*Axx + Gamyzz*Axy + Gamzzz*Axz + &
                 Gamxxz*Axz + Gamyxz*Ayz + Gamzxz*Azz) - chiz/chipn1*Axz + F1o3*gxz*fz
  Ayzx = Ayzx - (Gamxxz*Axy + Gamyxz*Ayy + Gamzxz*Ayz + &
                 Gamxxy*Axz + Gamyxy*Ayz + Gamzxy*Azz) - chix/chipn1*Ayz + F1o3*gyz*fx
  Ayzy = Ayzy - (Gamxyz*Axy + Gamyyz*Ayy + Gamzyz*Ayz + &
                 Gamxyy*Axz + Gamyyy*Ayz + Gamzyy*Azz) - chiy/chipn1*Ayz + F1o3*gyz*fy
  Ayzz = Ayzz - (Gamxzz*Axy + Gamyzz*Ayy + Gamzzz*Ayz + &
                 Gamxyz*Axz + Gamyyz*Ayz + Gamzyz*Azz) - chiz/chipn1*Ayz + F1o3*gyz*fz
! symmetrize B_ij = v^k (D_k K_ij -D_j K_ik)  
  Bxx =(vy*(Axxy - Axyx) + vz*(Axxz - Axzx))*f
  Byy =(vx*(Ayyx - Axyy) + vz*(Ayyz - Ayzy))*f
  Bzz =(vx*(Azzx - Axzz) + vy*(Azzy - Ayzz))*f
  Bxy =(vx*(Axyx - (Axxy+Axyx)/TWO) + vz*(Axyz - (Axzy+Ayzx)/TWO))*f
  Bxz =(vx*(Axzx - (Axxz+Axzx)/TWO) + vy*(Axzy - (Axyz+Ayzx)/TWO))*f
  Byz =(vx*(Ayzx - (Axyz+Axzy)/TWO) + vy*(Ayzy - (Ayyz+Ayzy)/TWO))*f
! E_ij = R_ij - K_ik * K^k_j + K * K_ij

! K_ij up to chi^-1
  Axxx = Axx + F1o3*trK*gxx
  Axyx = Axy + F1o3*trK*gxy
  Axzx = Axz + F1o3*trK*gxz
  Ayyx = Ayy + F1o3*trK*gyy
  Ayzx = Ayz + F1o3*trK*gyz
  Azzx = Azz + F1o3*trK*gzz
! gup and A_ijk cancel a chi^-1
  Exx =       gupxx * Axxx * Axxx + gupyy * Axyx * Axyx + gupzz * Axzx * Axzx + &
       TWO * (gupxy * Axxx * Axyx + gupxz * Axxx * Axzx + gupyz * Axyx * Axzx)
  Eyy =       gupxx * Axyx * Axyx + gupyy * Ayyx * Ayyx + gupzz * Ayzx * Ayzx + &
       TWO * (gupxy * Axyx * Ayyx + gupxz * Axyx * Ayzx + gupyz * Ayyx * Ayzx)
  Ezz =       gupxx * Axzx * Axzx + gupyy * Ayzx * Ayzx + gupzz * Azzx * Azzx + &
       TWO * (gupxy * Axzx * Ayzx + gupxz * Axzx * Azzx + gupyz * Ayzx * Azzx)
  Exy =       gupxx * Axxx * Axyx + gupyy * Axyx * Ayyx + gupzz * Axzx * Ayzx + &
              gupxy *(Axxx * Ayyx + Axyx * Axyx)                            + &
              gupxz *(Axxx * Ayzx + Axzx * Axyx)                            + &
              gupyz *(Axyx * Ayzx + Axzx * Ayyx)
  Exz =       gupxx * Axxx * Axzx + gupyy * Axyx * Ayzx + gupzz * Axzx * Azzx + &
              gupxy *(Axxx * Ayzx + Axyx * Axzx)                            + &
              gupxz *(Axxx * Azzx + Axzx * Axzx)                            + &
              gupyz *(Axyx * Azzx + Axzx * Ayzx)
  Eyz =       gupxx * Axyx * Axzx + gupyy * Ayyx * Ayzx + gupzz * Ayzx * Azzx + &
              gupxy *(Axyx * Ayzx + Ayyx * Axzx)                            + &
              gupxz *(Axyx * Azzx + Ayzx * Axzx)                            + &
              gupyz *(Ayyx * Azzx + Ayzx * Ayzx)

  Exx = Rxx - (Exx - Axxx*trK)*f - Bxx
  Exy = Rxy - (Exy - Axyx*trK)*f - Bxy
  Exz = Rxz - (Exz - Axzx*trK)*f - Bxz
  Eyy = Ryy - (Eyy - Ayyx*trK)*f - Byy
  Eyz = Ryz - (Eyz - Ayzx*trK)*f - Byz
  Ezz = Rzz - (Ezz - Azzx*trK)*f - Bzz
!set m = (u - iw)/sqrt(2) following Frans, PRD 75, 124018(2007)
! compute uuww^ij = u^i * u^j - w^i * w^j
  uuwwxx = ux * ux - wx * wx
  uuwwxy = ux * uy - wx * wy
  uuwwxz = ux * uz - wx * wz
  uuwwyy = uy * uy - wy * wy
  uuwwyz = uy * uz - wy * wz
  uuwwzz = uz * uz - wz * wz

! compute uw^ij = u^i * w^j + w^i * u^j
  uwxx = ux * wx + wx * ux
  uwxy = ux * wy + wx * uy
  uwxz = ux * wz + wx * uz
  uwyy = uy * wy + wy * uy
  uwyz = uy * wz + wy * uz
  uwzz = uz * wz + wz * uz
!the real part of Psi4
  Rpsi4 =   Exx * uuwwxx + Eyy * uuwwyy + Ezz * uuwwzz &
         + (Exy * uuwwxy + Exz * uuwwxz + Eyz * uuwwyz) * TWO

!the imaginary part of Psi4
  Ipsi4 =   Exx * uwxx + Eyy * uwyy + Ezz * uwzz &
         + (Exy * uwxy + Exz * uwxz + Eyz * uwyz) * TWO
         
!multiply with -1/2
  Rpsi4 = - Rpsi4/TWO
  Ipsi4 = - Ipsi4/TWO

  Rpsi4 = dexp(-FOUR*dsqrt(PI/3)*Sphi)*Rpsi4
  Ipsi4 = dexp(-FOUR*dsqrt(PI/3)*Sphi)*Ipsi4

  return

  end subroutine getnp4
!-----------------------------------------------------------------------------
!
! compute the Newman-Penrose Weyl scalar Psi4
! for BSSN dynamical variables for shell
!
!-----------------------------------------------------------------------------

  subroutine getnp4_ss(ex,crho,sigma,R, X, Y, Z,                               &
               drhodx, drhody, drhodz,                                         &
               dsigmadx,dsigmady,dsigmadz,                                     &
               dRdx,dRdy,dRdz,                                                 &
               drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                &
               dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,    &
               dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,                            &
               chi, trK, &
               dxx,gxy,gxz,dyy,gyz,dzz, &
               Axx,Axy,Axz,Ayy,Ayz,Azz, &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,&
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,&
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,&
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,&
               Rpsi4, Ipsi4, &
               symmetry,sst)

  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3),symmetry,sst
  double precision,intent(in),dimension(ex(1))::crho
  double precision,intent(in),dimension(ex(2))::sigma
  double precision,intent(in),dimension(ex(3))::R
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: X,Y,Z
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::drhodx, drhody, drhodz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::dsigmadx,dsigmady,dsigmadz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::dRdx,dRdy,dRdz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dxx,gxy,gxz,dyy,gyz,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: chi,trK
! physical second kind of connection  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: Gamzyy, Gamzyz, Gamzzz
! physical Ricci tensor  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out):: Rpsi4,Ipsi4

!~~~~~~> Other variables:

  real*8, dimension(ex(1),ex(2),ex(3)) :: f,fx,fy,fz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: chix,chiy,chiz,chipn1
  real*8, dimension(ex(1),ex(2),ex(3)) :: vx,vy,vz,ux,uy,uz,wx,wy,wz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Exx,Exy,Exz,Eyy,Eyz,Ezz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Bxx,Bxy,Bxz,Byy,Byz,Bzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Axxx,Axxy,Axxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Axyx,Axyy,Axyz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Axzx,Axzy,Axzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Ayyx,Ayyy,Ayyz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Ayzx,Ayzy,Ayzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Azzx,Azzy,Azzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupyy,gupyz,gupzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: uuwwxx,uuwwxy,uuwwxz,uuwwyy,uuwwyz,uuwwzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: uwxx,uwxy,uwxz,uwyy,uwyz,uwzz

  real*8, parameter :: ZEO = 0.d0, ONE = 1.d0, TWO = 2.d0
  real*8, parameter :: F1o3 = 1.d0/3.d0, FOUR = 4.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  integer::i,j,k
  real*8,parameter::TINYRR=1.d-14

  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE
  chipn1= chi + ONE

! invert tilted metric
  gupzz =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
           gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / gupzz
  gupxy = - ( gxy * gzz - gyz * gxz ) / gupzz
  gupxz =   ( gxy * gyz - gyy * gxz ) / gupzz
  gupyy =   ( gxx * gzz - gxz * gxz ) / gupzz
  gupyz = - ( gxx * gyz - gxy * gxz ) / gupzz
  gupzz =   ( gxx * gyy - gxy * gxy ) / gupzz

! initialize U, V, W vetors	
  do i=1,ex(1)
  do j=1,ex(2)
  do k=1,ex(3)
     if(abs(X(i,j,k)) < TINYRR .and. abs(Y(i,j,k)) < TINYRR .and. abs(Z(i,j,k)) < TINYRR)then
        vx(i,j,k) = TINYRR
        vy(i,j,k) = TINYRR
        vz(i,j,k) = TINYRR
     else
        vx(i,j,k) = X(i,j,k)
        vy(i,j,k) = Y(i,j,k)
        vz(i,j,k) = Z(i,j,k)
     endif
     if(abs(X(i,j,k)) < TINYRR .and. abs(Y(i,j,k)) < TINYRR)then
        ux(i,j,k) = - TINYRR
        uy(i,j,k) = TINYRR
        uz(i,j,k) = ZEO
        wx(i,j,k) = TINYRR*Z(i,j,k)
        wy(i,j,k) = TINYRR*Z(i,j,k)
        wz(i,j,k) = -2*TINYRR*TINYRR
     else
        ux(i,j,k) = - Y(i,j,k)
        uy(i,j,k) = X(i,j,k)
        uz(i,j,k) = ZEO
        wx(i,j,k) = X(i,j,k)*Z(i,j,k)
        wy(i,j,k) = Y(i,j,k)*Z(i,j,k)
        wz(i,j,k) = -(X(i,j,k)*X(i,j,k) + Y(i,j,k)*Y(i,j,k))
     endif
  enddo
  enddo
  enddo

  f = 1.d0/chipn1

  fx = gxx*vx*vx + gyy*vy*vy + gzz*vz*vz &
     +(gxy*vx*vy + gxz*vx*vz + gyz*vy*vz)*TWO
  fx = dsqrt(fx*f)
  vx = vx/fx
  vy = vy/fx
  vz = vz/fx

  fx = gxx*vx*ux + gxy*vx*uy + gxz*vx*uz + &
       gxy*vy*ux + gyy*vy*uy + gyz*vy*uz + &
       gxz*vz*ux + gyz*vz*uy + gzz*vz*uz
  fx = fx*f
  ux = ux - fx*vx
  uy = uy - fx*vy
  uz = uz - fx*vz
  fx = gxx*ux*ux + gyy*uy*uy + gzz*uz*uz &
     +(gxy*ux*uy + gxz*ux*uz + gyz*uy*uz)*TWO
  fx = dsqrt(fx*f) 
  ux = ux/fx
  uy = uy/fx
  uz = uz/fx

  fx = gxx*vx*wx + gxy*vx*wy + gxz*vx*wz + &
       gxy*vy*wx + gyy*vy*wy + gyz*vy*wz + &
       gxz*vz*wx + gyz*vz*wy + gzz*vz*wz
  fx = fx*f       
  wx = wx - fx*vx
  wy = wy - fx*vy
  wz = wz - fx*vz
  fx = gxx*ux*wx + gxy*ux*wy + gxz*ux*wz + &
       gxy*uy*wx + gyy*uy*wy + gyz*uy*wz + &
       gxz*uz*wx + gyz*uz*wy + gzz*uz*wz
  fx = fx*f       
  wx = wx - fx*ux
  wy = wy - fx*uy
  wz = wz - fx*uz
  fx = gxx*wx*wx + gyy*wy*wy + gzz*wz*wz &
     +(gxy*wx*wy + gxz*wx*wz + gyz*wy*wz)*TWO
  fx = dsqrt(fx*f)
  wx = wx/fx
  wy = wy/fx
  wz = wz/fx

  call fderivs_shc(ex,Axx,Axxx,Axxy,Axxz,crho,sigma,R, SYM, SYM,SYM,Symmetry,0,sst,         &
                       drhodx, drhody, drhodz,                                              &
                       dsigmadx,dsigmady,dsigmadz,                                          &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,Axy,Axyx,Axyy,Axyz,crho,sigma,R,ANTI,ANTI,SYM,Symmetry,0,sst,         &
                       drhodx, drhody, drhodz,                                              &
                       dsigmadx,dsigmady,dsigmadz,                                          &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,Axz,Axzx,Axzy,Axzz,crho,sigma,R,ANTI,SYM ,ANTI,Symmetry,0,sst,        &
                       drhodx, drhody, drhodz,                                              &
                       dsigmadx,dsigmady,dsigmadz,                                          &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,Ayy,Ayyx,Ayyy,Ayyz,crho,sigma,R, SYM, SYM,SYM,Symmetry,0,sst,         &
                       drhodx, drhody, drhodz,                                              &
                       dsigmadx,dsigmady,dsigmadz,                                          &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,Ayz,Ayzx,Ayzy,Ayzz,crho,sigma,R,SYM ,ANTI,ANTI,Symmetry,0,sst,        &
                       drhodx, drhody, drhodz,                                              &
                       dsigmadx,dsigmady,dsigmadz,                                          &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,Azz,Azzx,Azzy,Azzz,crho,sigma,R, SYM, SYM,SYM,Symmetry,0,sst,         &
                       drhodx, drhody, drhodz,                                              &
                       dsigmadx,dsigmady,dsigmadz,                                          &
                       dRdx,dRdy,dRdz)

  call fderivs_shc(ex,chi,chix,chiy,chiz,crho,sigma,R, SYM, SYM,SYM,Symmetry,0,sst,         &
                       drhodx, drhody, drhodz,                                              &
                       dsigmadx,dsigmady,dsigmadz,                                          &
                       dRdx,dRdy,dRdz)

  call fderivs_shc(ex,trK,fx,fy,fz,crho,sigma,R, SYM, SYM,SYM,Symmetry,0,sst,         &
                       drhodx, drhody, drhodz,                                        &
                       dsigmadx,dsigmady,dsigmadz,                                    &
                       dRdx,dRdy,dRdz)
! compute D_k K_ij up to chi^-1
  Axxx = Axxx - (Gamxxx*Axx + Gamyxx*Axy + Gamzxx*Axz)*TWO - chix/chipn1*Axx + F1o3*gxx*fx
  Axxy = Axxy - (Gamxxy*Axx + Gamyxy*Axy + Gamzxy*Axz)*TWO - chiy/chipn1*Axx + F1o3*gxx*fy
  Axxz = Axxz - (Gamxxz*Axx + Gamyxz*Axy + Gamzxz*Axz)*TWO - chiz/chipn1*Axx + F1o3*gxx*fz
  Ayyx = Ayyx - (Gamxxy*Axy + Gamyxy*Ayy + Gamzxy*Ayz)*TWO - chix/chipn1*Ayy + F1o3*gyy*fx
  Ayyy = Ayyy - (Gamxyy*Axy + Gamyyy*Ayy + Gamzyy*Ayz)*TWO - chiy/chipn1*Ayy + F1o3*gyy*fy
  Ayyz = Ayyz - (Gamxyz*Axy + Gamyyz*Ayy + Gamzyz*Ayz)*TWO - chiz/chipn1*Ayy + F1o3*gyy*fz
  Azzx = Azzx - (Gamxxz*Axz + Gamyxz*Ayz + Gamzxz*Azz)*TWO - chix/chipn1*Azz + F1o3*gzz*fx
  Azzy = Azzy - (Gamxyz*Axz + Gamyyz*Ayz + Gamzyz*Azz)*TWO - chiy/chipn1*Azz + F1o3*gzz*fy
  Azzz = Azzz - (Gamxzz*Axz + Gamyzz*Ayz + Gamzzz*Azz)*TWO - chiz/chipn1*Azz + F1o3*gzz*fz
  Axyx = Axyx - (Gamxxy*Axx + Gamyxy*Axy + Gamzxy*Axz + &
                 Gamxxx*Axy + Gamyxx*Ayy + Gamzxx*Ayz) - chix/chipn1*Axy + F1o3*gxy*fx
  Axyy = Axyy - (Gamxyy*Axx + Gamyyy*Axy + Gamzyy*Axz + &
                 Gamxxy*Axy + Gamyxy*Ayy + Gamzxy*Ayz) - chiy/chipn1*Axy + F1o3*gxy*fy
  Axyz = Axyz - (Gamxyz*Axx + Gamyyz*Axy + Gamzyz*Axz + &
                 Gamxxz*Axy + Gamyxz*Ayy + Gamzxz*Ayz) - chiz/chipn1*Axy + F1o3*gxy*fz
  Axzx = Axzx - (Gamxxz*Axx + Gamyxz*Axy + Gamzxz*Axz + &
                 Gamxxx*Axz + Gamyxx*Ayz + Gamzxx*Azz) - chix/chipn1*Axz + F1o3*gxz*fx
  Axzy = Axzy - (Gamxyz*Axx + Gamyyz*Axy + Gamzyz*Axz + &
                 Gamxxy*Axz + Gamyxy*Ayz + Gamzxy*Azz) - chiy/chipn1*Axz + F1o3*gxz*fy
  Axzz = Axzz - (Gamxzz*Axx + Gamyzz*Axy + Gamzzz*Axz + &
                 Gamxxz*Axz + Gamyxz*Ayz + Gamzxz*Azz) - chiz/chipn1*Axz + F1o3*gxz*fz
  Ayzx = Ayzx - (Gamxxz*Axy + Gamyxz*Ayy + Gamzxz*Ayz + &
                 Gamxxy*Axz + Gamyxy*Ayz + Gamzxy*Azz) - chix/chipn1*Ayz + F1o3*gyz*fx
  Ayzy = Ayzy - (Gamxyz*Axy + Gamyyz*Ayy + Gamzyz*Ayz + &
                 Gamxyy*Axz + Gamyyy*Ayz + Gamzyy*Azz) - chiy/chipn1*Ayz + F1o3*gyz*fy
  Ayzz = Ayzz - (Gamxzz*Axy + Gamyzz*Ayy + Gamzzz*Ayz + &
                 Gamxyz*Axz + Gamyyz*Ayz + Gamzyz*Azz) - chiz/chipn1*Ayz + F1o3*gyz*fz
! symmetrize B_ij = v^k (D_k K_ij -D_j K_ik)  
  Bxx =(vy*(Axxy - Axyx) + vz*(Axxz - Axzx))*f
  Byy =(vx*(Ayyx - Axyy) + vz*(Ayyz - Ayzy))*f
  Bzz =(vx*(Azzx - Axzz) + vy*(Azzy - Ayzz))*f
  Bxy =(vx*(Axyx - (Axxy+Axyx)/TWO) + vz*(Axyz - (Axzy+Ayzx)/TWO))*f
  Bxz =(vx*(Axzx - (Axxz+Axzx)/TWO) + vy*(Axzy - (Axyz+Ayzx)/TWO))*f
  Byz =(vx*(Ayzx - (Axyz+Axzy)/TWO) + vy*(Ayzy - (Ayyz+Ayzy)/TWO))*f
! E_ij = R_ij - K_ik * K^k_j + K * K_ij

! K_ij up to chi^-1
  Axxx = Axx + F1o3*trK*gxx
  Axyx = Axy + F1o3*trK*gxy
  Axzx = Axz + F1o3*trK*gxz
  Ayyx = Ayy + F1o3*trK*gyy
  Ayzx = Ayz + F1o3*trK*gyz
  Azzx = Azz + F1o3*trK*gzz
! gup and A_ijk cancel a chi^-1
  Exx =       gupxx * Axxx * Axxx + gupyy * Axyx * Axyx + gupzz * Axzx * Axzx + &
       TWO * (gupxy * Axxx * Axyx + gupxz * Axxx * Axzx + gupyz * Axyx * Axzx)
  Eyy =       gupxx * Axyx * Axyx + gupyy * Ayyx * Ayyx + gupzz * Ayzx * Ayzx + &
       TWO * (gupxy * Axyx * Ayyx + gupxz * Axyx * Ayzx + gupyz * Ayyx * Ayzx)
  Ezz =       gupxx * Axzx * Axzx + gupyy * Ayzx * Ayzx + gupzz * Azzx * Azzx + &
       TWO * (gupxy * Axzx * Ayzx + gupxz * Axzx * Azzx + gupyz * Ayzx * Azzx)
  Exy =       gupxx * Axxx * Axyx + gupyy * Axyx * Ayyx + gupzz * Axzx * Ayzx + &
              gupxy *(Axxx * Ayyx + Axyx * Axyx)                            + &
              gupxz *(Axxx * Ayzx + Axzx * Axyx)                            + &
              gupyz *(Axyx * Ayzx + Axzx * Ayyx)
  Exz =       gupxx * Axxx * Axzx + gupyy * Axyx * Ayzx + gupzz * Axzx * Azzx + &
              gupxy *(Axxx * Ayzx + Axyx * Axzx)                            + &
              gupxz *(Axxx * Azzx + Axzx * Axzx)                            + &
              gupyz *(Axyx * Azzx + Axzx * Ayzx)
  Eyz =       gupxx * Axyx * Axzx + gupyy * Ayyx * Ayzx + gupzz * Ayzx * Azzx + &
              gupxy *(Axyx * Ayzx + Ayyx * Axzx)                            + &
              gupxz *(Axyx * Azzx + Ayzx * Axzx)                            + &
              gupyz *(Ayyx * Azzx + Ayzx * Ayzx)

  Exx = Rxx - (Exx - Axxx*trK)*f - Bxx
  Exy = Rxy - (Exy - Axyx*trK)*f - Bxy
  Exz = Rxz - (Exz - Axzx*trK)*f - Bxz
  Eyy = Ryy - (Eyy - Ayyx*trK)*f - Byy
  Eyz = Ryz - (Eyz - Ayzx*trK)*f - Byz
  Ezz = Rzz - (Ezz - Azzx*trK)*f - Bzz
!set m = (u - iw)/sqrt(2) following Frans, PRD 75, 124018(2007)
! compute uuww^ij = u^i * u^j - w^i * w^j
  uuwwxx = ux * ux - wx * wx
  uuwwxy = ux * uy - wx * wy
  uuwwxz = ux * uz - wx * wz
  uuwwyy = uy * uy - wy * wy
  uuwwyz = uy * uz - wy * wz
  uuwwzz = uz * uz - wz * wz

! compute uw^ij = u^i * w^j + w^i * u^j
  uwxx = ux * wx + wx * ux
  uwxy = ux * wy + wx * uy
  uwxz = ux * wz + wx * uz
  uwyy = uy * wy + wy * uy
  uwyz = uy * wz + wy * uz
  uwzz = uz * wz + wz * uz
!the real part of Psi4
  Rpsi4 =   Exx * uuwwxx + Eyy * uuwwyy + Ezz * uuwwzz &
         + (Exy * uuwwxy + Exz * uuwwxz + Eyz * uuwwyz) * TWO

!the imaginary part of Psi4
  Ipsi4 =   Exx * uwxx + Eyy * uwyy + Ezz * uwzz &
         + (Exy * uwxy + Exz * uwxz + Eyz * uwyz) * TWO
         
!multiply with -1/2
  Rpsi4 = - Rpsi4/TWO
  Ipsi4 = - Ipsi4/TWO

  return

  end subroutine getnp4_ss
