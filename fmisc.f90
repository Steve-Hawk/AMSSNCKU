!$Id: fmisc.f90,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $

#include "microdef.fh"

#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
!---------------------------------------------------------------------------------------------------
! copy a point of data into data target for vertext center code
!---------------------------------------------------------------------------------------------------
  subroutine pointcopy(wei,llbout,uubout,ext_out,data_out,xx,yy,zz,dv)
  implicit none
  integer,intent(in) :: wei
  integer,dimension(3),intent(in) ::ext_out
  real*8,dimension(3) :: llbout,uubout
  real*8,dimension(ext_out(1),ext_out(2),ext_out(3)),intent(inout)::data_out
  real*8,intent(in) :: xx,yy,zz,dv

  real*8,dimension(3) :: ho
  integer :: i,j,k

!sanity check
  if(wei.ne.3)then
     write(*,*)"fmisc.f90::pointcopy: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

!!!   
  if(any(ext_out == 1))then
       write(*,*)"fmisc.f90::pointcopy: meets iolated points for out data"
       write(*,*) llbout,uubout
       stop
  else
    ho = (uubout-llbout)/(ext_out-1)
  endif    
  i = idint((xx-llbout(1))/ho(1)+0.4)+1
  j = idint((yy-llbout(2))/ho(2)+0.4)+1
  k = idint((zz-llbout(3))/ho(3)+0.4)+1

  if(i<1 .or. i>ext_out(1) .or. &
     j<1 .or. j>ext_out(2) .or. &
     k<1 .or. k>ext_out(3) )then
     write(*,*)"i,j,k = ",i,j,k
     write(*,*)"ext = ",ext_out
     stop
  endif
  if(dabs(llbout(1)+(i-1)*ho(1)-xx)>ho(1)/2 .or. &
     dabs(llbout(2)+(j-1)*ho(2)-yy)>ho(2)/2 .or. &
     dabs(llbout(3)+(k-1)*ho(3)-zz)>ho(3)/2 )then    
     write(*,*)"fmisc.f90::pointcopy: llbout = ",llbout
     write(*,*)"fmisc.f90::pointcopy: ho = ",ho
     write(*,*)"fmisc.f90::pointcopy: x,y,z = ",llbout(1)+(i-1)*ho(1),llbout(2)+(j-1)*ho(2),llbout(3)+(k-1)*ho(3)
     write(*,*)"fmisc.f90::pointcopy: point = ",xx,yy,zz
     stop
   endif

  data_out(i,j,k)=dv

  return

  end subroutine pointcopy
!---------------------------------------------------------------------------------------------------
! copy a part of data from data source, for vertex center code
!---------------------------------------------------------------------------------------------------
  subroutine copy(wei,llbout,uubout,ext_out,data_out,llbin,uubin,ext_in,data_in,lcopy,ucopy)
  implicit none
  integer,intent(in) :: wei
  integer,dimension(3),intent(in) ::ext_out,ext_in
  real*8,dimension(3),intent(in) :: lcopy,ucopy
  real*8,dimension(3) :: llbout,uubout,llbin,uubin
  real*8,dimension(ext_out(1),ext_out(2),ext_out(3)),intent(inout)::data_out
  real*8,dimension(ext_in(1),ext_in(2),ext_in(3)),intent(in)::data_in

  real*8,dimension(3) :: ho,hi
  integer,dimension(3) :: illo,iuuo,illi,iuui

!sanity check
  if(wei.ne.3)then
     write(*,*)"fmisc.f90::copy: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

!!!   
  if(any(ext_out == 1))then
    if(any(ext_in == 1))then
       write(*,*)"fmisc.f90::copy: meets iolated points for both in and out data"
       write(*,*) llbin,uubin
       write(*,*) llbout,uubout
       stop
     else
       hi = (uubin-llbin)/(ext_in-1)
       ho = hi
    endif
  else
    ho = (uubout-llbout)/(ext_out-1)
    if(any(ext_in == 1))then
      hi = ho
    else
      hi = (uubin-llbin)/(ext_in-1)
      if(any(abs(hi-ho) > min(hi,ho)/2))then
         write(*,*)"fmisc.f90::copy: meets copy reqest for different numerical grid"
         write(*,*)hi,ho
         stop
      endif
    endif
  endif    
  illo = idint((lcopy-llbout)/ho+0.4)+1
  iuuo = ext_out - idint((uubout-ucopy)/ho+0.4)
  illi = idint((lcopy-llbin)/hi+0.4)+1
  iuui = ext_in - idint((uubin-ucopy)/hi+0.4)

  if(any(llbout-lcopy>ho/2) .or. any(ucopy-uubout>ho/2))then
     write(*,*)"fmisc.f90::copy: llbout = ",llbout
     write(*,*)"fmisc.f90::copy: uubout = ",uubout
     write(*,*)"fmisc.f90::copy: llbcopy = ",lcopy
     write(*,*)"fmisc.f90::copy: uubcopy = ",ucopy
     write(*,*)"fmisc.f90::copy: ho = ",ho
     write(*,*)llbout-lcopy,ucopy-uubout
     stop
  elseif(any(llbin -lcopy>hi/2) .or. any(ucopy-uubin >hi/2))then
     write(*,*)"fmisc.f90::copy: llbin = ",llbin
     write(*,*)"fmisc.f90::copy: uubin = ",uubin
     write(*,*)"fmisc.f90::copy: llbcopy = ",lcopy
     write(*,*)"fmisc.f90::copy: uubcopy = ",ucopy
     stop
  elseif(any(illo<1) .or. any(illi<1) .or. any(illo-iuuo>0) .or. any(illi-iuui>0) .or. &
         any(iuui-ext_in>0) .or. any(iuuo-ext_out>0))then
     write(*,*)"fmisc.f90::copy: illi = ",illi
     write(*,*)"fmisc.f90::copy: iuui = ",iuui
     write(*,*)"fmisc.f90::copy: illo = ",illo
     write(*,*)"fmisc.f90::copy: iuuo = ",iuuo     
     write(*,*)"fmisc.f90::copy: llbout = ",llbout
     write(*,*)"fmisc.f90::copy: uubout = ",uubout
     write(*,*)"fmisc.f90::copy: llbin = ",llbin
     write(*,*)"fmisc.f90::copy: uubin = ",uubin
     write(*,*)"fmisc.f90::copy: llbcopy = ",lcopy
     write(*,*)"fmisc.f90::copy: uubcopy = ",ucopy
     stop
   endif

  data_out(illo(1):iuuo(1),illo(2):iuuo(2),illo(3):iuuo(3))=data_in(illi(1):iuui(1),illi(2):iuui(2),illi(3):iuui(3))

  return

  end subroutine copy
!-----------------------------------------------------------------------------------------------------------------
! three dimensional interpolation for vertex center grid structure  
  subroutine global_interp(ex,X,Y,Z,f,f_int,x1,y1,z1,ORDN,SoA,symmetry)
  implicit none

!~~~~~~> Input parameters:

  integer,                             intent(in) :: ex(1:3), symmetry,ORDN
  real*8,intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8,                              intent(out):: f_int
  real*8,                              intent(in) :: x1,y1,z1
  real*8, dimension(3),                intent(in) :: SoA

!~~~~~~> Other parameters:

  integer :: j,m,imin,jmin,kmin
  integer,dimension(3) :: cxB,cxT,cxI,cmin,cmax
  real*8,dimension(3) :: cx
  real*8, dimension(1:ORDN) :: x1a
  real*8, dimension(1:ORDN,1:ORDN,1:ORDN) :: ya
  integer, parameter :: NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2
  real*8 :: dX,dY,dZ,ddy
  real*8, parameter :: ONE=1.d0
  logical::decide3d

  imin = lbound(f,1)
  jmin = lbound(f,2)
  kmin = lbound(f,3)

  dX = X(imin+1)-X(imin)
  dY = Y(jmin+1)-Y(jmin)
  dZ = Z(kmin+1)-Z(kmin)

  forall( j = 1:ordn ) x1a(j) = ( j - 1 )* ONE

  cxI(1) = idint((x1-X(1))/dX+0.4)+1
  cxI(2) = idint((y1-Y(1))/dY+0.4)+1
  cxI(3) = idint((z1-Z(1))/dZ+0.4)+1

  cxB = cxI - ORDN/2+1
  cxT = cxB + ORDN - 1
       
  cmin = 1
  cmax = ex
  if(Symmetry == OCTANT  .and.dabs(X(1))<dX) cmin(1) = -ORDN/2+2
  if(Symmetry == OCTANT  .and.dabs(Y(1))<dY) cmin(2) = -ORDN/2+2
  if(Symmetry /= NO_SYMM .and.dabs(Z(1))<dZ) cmin(3) = -ORDN/2+2
  do m =1,3
   if(cxB(m) < cmin(m))then
      cxB(m) = cmin(m)
      cxT(m) = cxB(m) + ORDN - 1
   endif
   if(cxT(m) > cmax(m))then
      cxT(m) = cmax(m)
      cxB(m) = cxT(m) + 1 - ORDN
   endif
 enddo
 if(cxB(1)>0)then
  cx(1) = (x1 - X(cxB(1)))/dX
 else
  cx(1) = (x1 + X(2-cxB(1)))/dX
 endif
 if(cxB(2)>0)then
  cx(2) = (y1 - Y(cxB(2)))/dY
 else
  cx(2) = (y1 + Y(2-cxB(2)))/dY
 endif
 if(cxB(3)>0)then
  cx(3) = (z1 - Z(cxB(3)))/dZ
 else
  cx(3) = (z1 + Z(2-cxB(3)))/dZ
 endif

  if(decide3d(ex,f,f,cxB,cxT,SoA,ya,ORDN,Symmetry))then
     write(*,*)"global_interp position: ",x1,y1,z1
     write(*,*)"data range: ",X(1),X(ex(1)),Y(1),Y(ex(2)),Z(1),Z(ex(3))
     stop
  endif
  call polin3(x1a,x1a,x1a,ya,cx(1),cx(2),cx(3),f_int,ddy,ORDN)

  return

  end subroutine global_interp
!----------------------------------------------------------------
! decide which 3d data to be used does not surport PI-Symmetry yet 
!----------------------------------------------------------------
  function decide3d(ex,f,fpi,cxB,cxT,SoA,ya,ORDN,Symmetry)  result(gont)
  implicit none

  integer,                                 intent(in) :: ORDN,Symmetry
  integer,dimension(1:3)                 , intent(in) :: ex,cxB,cxT
  real*8, dimension(1:3)                 , intent(in) :: SoA
  real*8, dimension(ex(1),ex(2),ex(3))   , intent(in) :: f,fpi
  real*8, dimension(cxB(1):cxT(1),cxB(2):cxT(2),cxB(3):cxT(3)), intent(out):: ya
  logical::gont

  integer,dimension(1:3) :: fmin1,fmin2,fmax1,fmax2
  integer::i,j,k,m

  gont=.false.
  do m=1,3
! check cxB and cxT are NaN or not  
    if(.not.(iabs(cxB(m)).ge.0)) gont=.true.
    if(.not.(iabs(cxT(m)).ge.0)) gont=.true.
    fmin1(m) = max(1,cxB(m))
    fmax1(m) = cxT(m)
    fmin2(m) = cxB(m)
    fmax2(m) = min(0,cxT(m))
    if((fmin1(m).le.fmax1(m)).and.(  fmin1(m)<1.or.  fmax1(m)>ex(m)))gont=.true.
    if((fmin2(m).le.fmax2(m)).and.(2-fmax2(m)<1.or.2-fmin2(m)>ex(m)))gont=.true.
  enddo
!sanity check
  if(gont)then
          write(*,*)"error in decide3d"
          write(*,*)((fmin1.le.fmax1).and.(  fmin1<1.or.  fmax1>ex))
          write(*,*)((fmin2.le.fmax2).and.(2-fmax2<1.or.2-fmin2>ex))
          write(*,*)"cxB, cxT and data shape:"
          write(*,*)cxB,cxT,ex
          write(*,*)"resulted fmin1, fmax1 and fmin2, fmax2:"
          write(*,*)fmin1,fmax1,fmin2,fmax2
  else

  do k=fmin1(3),fmax1(3)
     do j=fmin1(2),fmax1(2)
        do i=fmin1(1),fmax1(1)
           ya(i,j,k) = f(i,j,k)
        enddo
        do i=fmin2(1),fmax2(1)
           ya(i,j,k) = f(2-i,j,k)*SoA(1)
        enddo
     enddo
     do j=fmin2(2),fmax2(2)
        do i=fmin1(1),fmax1(1)
           ya(i,j,k) = f(i,2-j,k)*SoA(2)
        enddo
        do i=fmin2(1),fmax2(1)
           ya(i,j,k) = f(2-i,2-j,k)*SoA(1)*SoA(2)
        enddo
     enddo
  enddo
           
  do k=fmin2(3),fmax2(3)
     do j=fmin1(2),fmax1(2)
        do i=fmin1(1),fmax1(1)
           ya(i,j,k) = f(i,j,2-k)*SoA(3)
        enddo
        do i=fmin2(1),fmax2(1)
           ya(i,j,k) = f(2-i,j,2-k)*SoA(1)*SoA(3)
        enddo
     enddo
     do j=fmin2(2),fmax2(2)
        do i=fmin1(1),fmax1(1)
           ya(i,j,k) = f(i,2-j,2-k)*SoA(2)*SoA(3)
        enddo
        do i=fmin2(1),fmax2(1)
           ya(i,j,k) = f(2-i,2-j,2-k)*SoA(1)*SoA(2)*SoA(3)
        enddo
     enddo
  enddo

  endif

  end function decide3d

!---------------------------------------------------------------------------------------
subroutine symmetry_bd(ord,extc,func,funcc,SoA)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: ord
  integer,dimension(3),   intent(in) :: extc
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in ):: func
  real*8, dimension(-ord+1:extc(1),-ord+1:extc(2),-ord+1:extc(3)),intent(out):: funcc
  real*8, dimension(1:3), intent(in) :: SoA

  integer::i

  funcc = 0.d0
  funcc(1:extc(1),1:extc(2),1:extc(3)) = func
   do i=0,ord-1
      funcc(-i,1:extc(2),1:extc(3)) = funcc(i+2,1:extc(2),1:extc(3))*SoA(1)
   enddo
   do i=0,ord-1
      funcc(:,-i,1:extc(3)) = funcc(:,i+2,1:extc(3))*SoA(2)
   enddo
   do i=0,ord-1
      funcc(:,:,-i) = funcc(:,:,i+2)*SoA(3)
   enddo

end subroutine symmetry_bd

subroutine symmetry_tbd(ord,extc,func,funcc,SoA)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: ord
  integer,dimension(3),   intent(in) :: extc
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in ):: func
  real*8, dimension(-ord+1:extc(1)+ord,-ord+1:extc(2)+ord,-ord+1:extc(3)+ord),intent(out):: funcc
  real*8, dimension(1:3), intent(in) :: SoA

  integer::i

  funcc = 0.d0
  funcc(1:extc(1),1:extc(2),1:extc(3)) = func
   do i=0,ord-1
      funcc(-i,1:extc(2),1:extc(3)) = funcc(i+2,1:extc(2),1:extc(3))*SoA(1)
      funcc(extc(1)+1+i,1:extc(2),1:extc(3)) = funcc(extc(1)-1-i,1:extc(2),1:extc(3))*SoA(1)
   enddo
   do i=0,ord-1
      funcc(:,-i,1:extc(3)) = funcc(:,i+2,1:extc(3))*SoA(2)
      funcc(:,extc(2)+1+i,1:extc(3)) = funcc(:,extc(2)-1-i,1:extc(3))*SoA(2)
   enddo
   do i=0,ord-1
      funcc(:,:,-i) = funcc(:,:,i+2)*SoA(3)
      funcc(:,:,extc(3)+1+i) = funcc(:,:,extc(3)-1-i)*SoA(3)
   enddo

end subroutine symmetry_tbd

subroutine symmetry_stbd(ord,extc,func,funcc,SoA)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: ord
  integer,dimension(3),   intent(in) :: extc
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in ):: func
  real*8, dimension(-ord+1:extc(1)+ord,-ord+1:extc(2)+ord,extc(3)),intent(out):: funcc
  real*8, dimension(2), intent(in) :: SoA

  integer::i

  funcc = 0.d0
  funcc(1:extc(1),1:extc(2),1:extc(3)) = func
   do i=0,ord-1
      funcc(-i,1:extc(2),1:extc(3)) = funcc(i+2,1:extc(2),1:extc(3))*SoA(1)
      funcc(extc(1)+1+i,1:extc(2),1:extc(3)) = funcc(extc(1)-1-i,1:extc(2),1:extc(3))*SoA(1)
   enddo
   do i=0,ord-1
      funcc(:,-i,1:extc(3)) = funcc(:,i+2,1:extc(3))*SoA(2)
      funcc(:,extc(2)+1+i,1:extc(3)) = funcc(:,extc(2)-1-i,1:extc(3))*SoA(2)
   enddo

end subroutine symmetry_stbd

subroutine symmetry_sntbd(ord,extc,func,funcc,SoA,actd)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: ord,actd
  integer,dimension(3),   intent(in) :: extc
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in ):: func
  real*8, dimension(-ord+1:extc(1)+ord,-ord+1:extc(2)+ord,extc(3)),intent(out):: funcc
  real*8, intent(in) :: SoA

  integer::i

  funcc = 0.d0
  funcc(1:extc(1),1:extc(2),1:extc(3)) = func
  if(actd==0)then
   do i=0,ord-1
      funcc(-i,1:extc(2),1:extc(3)) = funcc(i+2,1:extc(2),1:extc(3))*SoA
      funcc(extc(1)+1+i,1:extc(2),1:extc(3)) = funcc(extc(1)-1-i,1:extc(2),1:extc(3))*SoA
   enddo
  elseif(actd==1)then
   do i=0,ord-1
      funcc(1:extc(1),-i,1:extc(3)) = funcc(1:extc(1),i+2,1:extc(3))*SoA
      funcc(1:extc(1),extc(2)+1+i,1:extc(3)) = funcc(1:extc(1),extc(2)-1-i,1:extc(3))*SoA
   enddo
  else
    write(*,*)"symmetry_sntbd: not recognized actd = ",actd
  endif

end subroutine symmetry_sntbd

#else
#ifdef Cell
!---------------------------------------------------------------------------------------------------
! copy a point of data into data target for vertext center code
!---------------------------------------------------------------------------------------------------
  subroutine pointcopy(wei,llbout,uubout,ext_out,data_out,xx,yy,zz,dv)
  implicit none
  integer,intent(in) :: wei
  integer,dimension(3),intent(in) ::ext_out
  real*8,dimension(3) :: llbout,uubout
  real*8,dimension(ext_out(1),ext_out(2),ext_out(3)),intent(inout)::data_out
  real*8,intent(in) :: xx,yy,zz,dv

  real*8,dimension(3) :: ho
  integer :: i,j,k

!sanity check
  if(wei.ne.3)then
     write(*,*)"fmisc.f90::pointcopy: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

!!!   
  ho = (uubout-llbout)/ext_out
  i = idint((xx-llbout(1))/ho(1)+0.4)+1
  j = idint((yy-llbout(2))/ho(2)+0.4)+1
  k = idint((zz-llbout(3))/ho(3)+0.4)+1

  if(i<1 .or. i>ext_out(1) .or. &
     j<1 .or. j>ext_out(2) .or. &
     k<1 .or. k>ext_out(3) )then
     write(*,*)"i,j,k = ",i,j,k
     write(*,*)"ext = ",ext_out
     stop
  endif
  if(dabs(llbout(1)+(i-0.5)*ho(1)-xx)>ho(1)/2 .or. &
     dabs(llbout(2)+(j-0.5)*ho(2)-yy)>ho(2)/2 .or. &
     dabs(llbout(3)+(k-0.5)*ho(3)-zz)>ho(3)/2 )then    
     write(*,*)"fmisc.f90::pointcopy: llbout = ",llbout
     write(*,*)"fmisc.f90::pointcopy: ho = ",ho
     write(*,*)"fmisc.f90::pointcopy: x,y,z = ",llbout(1)+(i-0.5)*ho(1),llbout(2)+(j-0.5)*ho(2),llbout(3)+(k-0.5)*ho(3)
     write(*,*)"fmisc.f90::pointcopy: point = ",xx,yy,zz
     stop
   endif

  data_out(i,j,k)=dv

  return

  end subroutine pointcopy
!---------------------------------------------------------------------------------------------------
! copy a part of data from data source, for cell center code
!---------------------------------------------------------------------------------------------------
  subroutine copy(wei,llbout,uubout,ext_out,data_out,llbin,uubin,ext_in,data_in,lcopy,ucopy)
  implicit none
  integer,intent(in) :: wei
  integer,dimension(3),intent(in) ::ext_out,ext_in
  real*8,dimension(3),intent(in) :: lcopy,ucopy
  real*8,dimension(3) :: llbout,uubout,llbin,uubin
  real*8,dimension(ext_out(1),ext_out(2),ext_out(3)),intent(inout)::data_out
  real*8,dimension(ext_in(1),ext_in(2),ext_in(3)),intent(in)::data_in

  real*8,dimension(3) :: ho,hi
  integer,dimension(3) :: illo,iuuo,illi,iuui

!sanity check
  if(wei.ne.3)then
     write(*,*)"fmisc.f90::copy: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

!!!   
  ho = (uubout-llbout)/ext_out
  hi = (uubin-llbin)/ext_in
  illo = idint((lcopy-llbout)/ho+0.4)+1
  iuuo = ext_out - idint((uubout-ucopy)/ho+0.4)
  illi = idint((lcopy-llbin)/hi+0.4)+1
  iuui = ext_in - idint((uubin-ucopy)/hi+0.4)

  if(any(llbout-lcopy>ho/2) .or. any(ucopy-uubout>ho/2))then
     write(*,*)"fmisc.f90::copy: llbout = ",llbout
     write(*,*)"fmisc.f90::copy: uubout = ",uubout
     write(*,*)"fmisc.f90::copy: llbcopy = ",lcopy
     write(*,*)"fmisc.f90::copy: uubcopy = ",ucopy
     write(*,*)"fmisc.f90::copy: ho = ",ho
     write(*,*)llbout-lcopy,ucopy-uubout
     stop
  elseif(any(llbin -lcopy>hi/2) .or. any(ucopy-uubin >hi/2))then
     write(*,*)"fmisc.f90::copy: llbin = ",llbin
     write(*,*)"fmisc.f90::copy: uubin = ",uubin
     write(*,*)"fmisc.f90::copy: llbcopy = ",lcopy
     write(*,*)"fmisc.f90::copy: uubcopy = ",ucopy
     stop
  elseif(any(illo<1) .or. any(illi<1) .or. any(illo-iuuo>0) .or. any(illi-iuui>0) .or. &
         any(iuui-ext_in>0) .or. any(iuuo-ext_out>0))then
     write(*,*)"fmisc.f90::copy: illi = ",illi
     write(*,*)"fmisc.f90::copy: iuui = ",iuui
     write(*,*)"fmisc.f90::copy: illo = ",illo
     write(*,*)"fmisc.f90::copy: iuuo = ",iuuo     
     write(*,*)"fmisc.f90::copy: llbout = ",llbout
     write(*,*)"fmisc.f90::copy: uubout = ",uubout
     write(*,*)"fmisc.f90::copy: llbin = ",llbin
     write(*,*)"fmisc.f90::copy: uubin = ",uubin
     write(*,*)"fmisc.f90::copy: llbcopy = ",lcopy
     write(*,*)"fmisc.f90::copy: uubcopy = ",ucopy
     stop
   endif

  data_out(illo(1):iuuo(1),illo(2):iuuo(2),illo(3):iuuo(3))=data_in(illi(1):iuui(1),illi(2):iuui(2),illi(3):iuui(3))

  return

  end subroutine copy
!--------------------------------------------------------------------------
! three dimensional interpolation for cell center grid structure  
  subroutine global_interp(ex,X,Y,Z,f,f_int,x1,y1,z1,ORDN,SoA,symmetry)
  implicit none

!~~~~~~> Input parameters:

  integer,                             intent(in) :: ex(1:3), symmetry,ORDN
  real*8,intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8,                              intent(out):: f_int
  real*8,                              intent(in) :: x1,y1,z1
  real*8, dimension(3),                intent(in) :: SoA

!~~~~~~> Other parameters:

  integer :: j,m,imin,jmin,kmin
  integer,dimension(3) :: cxB,cxT,cxI,cmin,cmax
  real*8,dimension(3) :: cx
  real*8, dimension(1:ORDN) :: x1a
  real*8, dimension(1:ORDN,1:ORDN,1:ORDN) :: ya
  integer, parameter :: NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2
  real*8 :: dX,dY,dZ,ddy
  real*8, parameter :: ONE=1.d0
  logical::decide3d

  imin = lbound(f,1)
  jmin = lbound(f,2)
  kmin = lbound(f,3)

  dX = X(imin+1)-X(imin)
  dY = Y(jmin+1)-Y(jmin)
  dZ = Z(kmin+1)-Z(kmin)

  forall( j = 1:ordn ) x1a(j) = ( j - 1 )* ONE

  cxI(1) = idint((x1-X(1))/dX+0.4)+1
  cxI(2) = idint((y1-Y(1))/dY+0.4)+1
  cxI(3) = idint((z1-Z(1))/dZ+0.4)+1

  cxB = cxI - ORDN/2+1
  cxT = cxB + ORDN - 1
       
  cmin = 1
  cmax = ex
  if(Symmetry == OCTANT  .and.dabs(X(1))<dX) cmin(1) = -ORDN/2+1
  if(Symmetry == OCTANT  .and.dabs(Y(1))<dY) cmin(2) = -ORDN/2+1
  if(Symmetry /= NO_SYMM .and.dabs(Z(1))<dZ) cmin(3) = -ORDN/2+1

  do m =1,3
   if(cxB(m) < cmin(m))then
      cxB(m) = cmin(m)
      cxT(m) = cxB(m) + ORDN - 1
   endif
   if(cxT(m) > cmax(m))then
      cxT(m) = cmax(m)
      cxB(m) = cxT(m) + 1 - ORDN
   endif
 enddo
 if(cxB(1)>0)then
  cx(1) = (x1 - X(cxB(1)))/dX
 else
  cx(1) = (x1 + X(1-cxB(1)))/dX
 endif
 if(cxB(2)>0)then
  cx(2) = (y1 - Y(cxB(2)))/dY
 else
  cx(2) = (y1 + Y(1-cxB(2)))/dY
 endif
 if(cxB(3)>0)then
  cx(3) = (z1 - Z(cxB(3)))/dZ
 else
  cx(3) = (z1 + Z(1-cxB(3)))/dZ
 endif

  if(decide3d(ex,f,f,cxB,cxT,SoA,ya,ORDN,Symmetry))then
     write(*,*)"global_interp position: ",x1,y1,z1
     write(*,*)"data range: ",X(1),X(ex(1)),Y(1),Y(ex(2)),Z(1),Z(ex(3))
     stop
  endif
  call polin3(x1a,x1a,x1a,ya,cx(1),cx(2),cx(3),f_int,ddy,ORDN)

  return

  end subroutine global_interp
!----------------------------------------------------------------
! decide which 3d data to be used does not surport PI-Symmetry yet 
!----------------------------------------------------------------
  function decide3d(ex,f,fpi,cxB,cxT,SoA,ya,ORDN,Symmetry)  result(gont)
  implicit none

  integer,                                 intent(in) :: ORDN,Symmetry
  integer,dimension(1:3)                 , intent(in) :: ex,cxB,cxT
  real*8, dimension(1:3)                 , intent(in) :: SoA
  real*8, dimension(ex(1),ex(2),ex(3))   , intent(in) :: f,fpi
  real*8, dimension(cxB(1):cxT(1),cxB(2):cxT(2),cxB(3):cxT(3)), intent(out):: ya
  logical::gont

  integer,dimension(1:3) :: fmin1,fmin2,fmax1,fmax2
  integer::i,j,k,m

  gont=.false.
  do m=1,3
! check cxB and cxT are NaN or not  
    if(.not.(iabs(cxB(m)).ge.0)) gont=.true.
    if(.not.(iabs(cxT(m)).ge.0)) gont=.true.
    fmin1(m) = max(1,cxB(m))
    fmax1(m) = cxT(m)
    fmin2(m) = cxB(m)
    fmax2(m) = min(0,cxT(m))
    if((fmin1(m).le.fmax1(m)).and.(  fmin1(m)<1.or.  fmax1(m)>ex(m)))gont=.true.
    if((fmin2(m).le.fmax2(m)).and.(1-fmax2(m)<1.or.1-fmin2(m)>ex(m)))gont=.true.
  enddo
!sanity check
  if(gont)then
          write(*,*)"error in decide3d"
          write(*,*)((fmin1.le.fmax1).and.(  fmin1<1.or.  fmax1>ex))
          write(*,*)((fmin2.le.fmax2).and.(1-fmax2<1.or.1-fmin2>ex))
          write(*,*)"cxB, cxT and data shape:"
          write(*,*)cxB,cxT,ex
          write(*,*)"resulted fmin1, fmax1 and fmin2, fmax2:"
          write(*,*)fmin1,fmax1,fmin2,fmax2
  else

  do k=fmin1(3),fmax1(3)
     do j=fmin1(2),fmax1(2)
        do i=fmin1(1),fmax1(1)
           ya(i,j,k) = f(i,j,k)
        enddo
        do i=fmin2(1),fmax2(1)
           ya(i,j,k) = f(1-i,j,k)*SoA(1)
        enddo
     enddo
     do j=fmin2(2),fmax2(2)
        do i=fmin1(1),fmax1(1)
           ya(i,j,k) = f(i,1-j,k)*SoA(2)
        enddo
        do i=fmin2(1),fmax2(1)
           ya(i,j,k) = f(1-i,1-j,k)*SoA(1)*SoA(2)
        enddo
     enddo
  enddo
           
  do k=fmin2(3),fmax2(3)
     do j=fmin1(2),fmax1(2)
        do i=fmin1(1),fmax1(1)
           ya(i,j,k) = f(i,j,1-k)*SoA(3)
        enddo
        do i=fmin2(1),fmax2(1)
           ya(i,j,k) = f(1-i,j,1-k)*SoA(1)*SoA(3)
        enddo
     enddo
     do j=fmin2(2),fmax2(2)
        do i=fmin1(1),fmax1(1)
           ya(i,j,k) = f(i,1-j,1-k)*SoA(2)*SoA(3)
        enddo
        do i=fmin2(1),fmax2(1)
           ya(i,j,k) = f(1-i,1-j,1-k)*SoA(1)*SoA(2)*SoA(3)
        enddo
     enddo
  enddo

  endif

  end function decide3d

!---------------------------------------------------------------------------------------
subroutine symmetry_bd(ord,extc,func,funcc,SoA)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: ord
  integer,dimension(3),   intent(in) :: extc
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in ):: func
  real*8, dimension(-ord+1:extc(1),-ord+1:extc(2),-ord+1:extc(3)),intent(out):: funcc
  real*8, dimension(1:3), intent(in) :: SoA

  integer::i

  funcc = 0.d0
  funcc(1:extc(1),1:extc(2),1:extc(3)) = func
   do i=0,ord-1
      funcc(-i,1:extc(2),1:extc(3)) = funcc(i+1,1:extc(2),1:extc(3))*SoA(1)
   enddo
   do i=0,ord-1
      funcc(:,-i,1:extc(3)) = funcc(:,i+1,1:extc(3))*SoA(2)
   enddo
   do i=0,ord-1
      funcc(:,:,-i) = funcc(:,:,i+1)*SoA(3)
   enddo

end subroutine symmetry_bd

subroutine symmetry_tbd(ord,extc,func,funcc,SoA)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: ord
  integer,dimension(3),   intent(in) :: extc
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in ):: func
  real*8, dimension(-ord+1:extc(1)+ord,-ord+1:extc(2)+ord,-ord+1:extc(3)+ord),intent(out):: funcc
  real*8, dimension(1:3), intent(in) :: SoA

  integer::i

  funcc = 0.d0
  funcc(1:extc(1),1:extc(2),1:extc(3)) = func
   do i=0,ord-1
      funcc(-i,1:extc(2),1:extc(3)) = funcc(i+1,1:extc(2),1:extc(3))*SoA(1)
      funcc(extc(1)+1+i,1:extc(2),1:extc(3)) = funcc(extc(1)-i,1:extc(2),1:extc(3))*SoA(1)
   enddo
   do i=0,ord-1
      funcc(:,-i,1:extc(3)) = funcc(:,i+1,1:extc(3))*SoA(2)
      funcc(:,extc(2)+1+i,1:extc(3)) = funcc(:,extc(2)-i,1:extc(3))*SoA(2)
   enddo
   do i=0,ord-1
      funcc(:,:,-i) = funcc(:,:,i+1)*SoA(3)
      funcc(:,:,extc(3)+1+i) = funcc(:,:,extc(3)-i)*SoA(3)
   enddo

end subroutine symmetry_tbd

subroutine symmetry_stbd(ord,extc,func,funcc,SoA)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: ord
  integer,dimension(3),   intent(in) :: extc
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in ):: func
  real*8, dimension(-ord+1:extc(1)+ord,-ord+1:extc(2)+ord,extc(3)),intent(out):: funcc
  real*8, dimension(2), intent(in) :: SoA

  integer::i

  funcc = 0.d0
  funcc(1:extc(1),1:extc(2),1:extc(3)) = func
   do i=0,ord-1
      funcc(-i,1:extc(2),1:extc(3)) = funcc(i+1,1:extc(2),1:extc(3))*SoA(1)
      funcc(extc(1)+1+i,1:extc(2),1:extc(3)) = funcc(extc(1)-i,1:extc(2),1:extc(3))*SoA(1)
   enddo
   do i=0,ord-1
      funcc(:,-i,1:extc(3)) = funcc(:,i+1,1:extc(3))*SoA(2)
      funcc(:,extc(2)+1+i,1:extc(3)) = funcc(:,extc(2)-i,1:extc(3))*SoA(2)
   enddo

end subroutine symmetry_stbd

subroutine symmetry_sntbd(ord,extc,func,funcc,SoA,actd)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: ord,actd
  integer,dimension(3),   intent(in) :: extc
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in ):: func
  real*8, dimension(-ord+1:extc(1)+ord,-ord+1:extc(2)+ord,extc(3)),intent(out):: funcc
  real*8, intent(in) :: SoA

  integer::i

  funcc = 0.d0
  funcc(1:extc(1),1:extc(2),1:extc(3)) = func
  if(actd==0)then
   do i=0,ord-1
      funcc(-i,1:extc(2),1:extc(3)) = funcc(i+1,1:extc(2),1:extc(3))*SoA
      funcc(extc(1)+1+i,1:extc(2),1:extc(3)) = funcc(extc(1)-i,1:extc(2),1:extc(3))*SoA
   enddo
  elseif(actd==1)then
   do i=0,ord-1
      funcc(1:extc(1),-i,1:extc(3)) = funcc(1:extc(1),i+1,1:extc(3))*SoA
      funcc(1:extc(1),extc(2)+1+i,1:extc(3)) = funcc(1:extc(1),extc(2)-i,1:extc(3))*SoA
   enddo
  else
    write(*,*)"symmetry_sntbd: not recognized actd = ",actd
  endif

end subroutine symmetry_sntbd

#else
#error Not define Vertex nor Cell
#endif  
#endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! common code for cell and vertex
!------------------------------------------------------------------------------
! Lagrangian polynomial interpolation
!------------------------------------------------------------------------------

  subroutine polint(xa,ya,x,y,dy,ordn)

  implicit none

!~~~~~~> Input Parameter:
  integer,intent(in) :: ordn
  real*8, dimension(ordn), intent(in) :: xa,ya
  real*8, intent(in) :: x
  real*8, intent(out) :: y,dy

!~~~~~~> Other parameter:

  integer :: m,n,ns
  real*8, dimension(ordn) :: c,d,den,ho
  real*8 :: dif,dift

!~~~~~~>

  n=ordn
  m=ordn

  c=ya
  d=ya
  ho=xa-x

  ns=1
  dif=abs(x-xa(1))
  do m=1,n
   dift=abs(x-xa(m))
   if(dift < dif) then
    ns=m
    dif=dift
   end if
  end do

  y=ya(ns)
  ns=ns-1
  do m=1,n-1
    den(1:n-m)=ho(1:n-m)-ho(1+m:n)
    if (any(den(1:n-m) == 0.0))then
      write(*,*) 'failure in polint for point',x
      stop
    endif
    den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
    d(1:n-m)=ho(1+m:n)*den(1:n-m)
    c(1:n-m)=ho(1:n-m)*den(1:n-m)
    if (2*ns < n-m) then
      dy=c(ns+1)
    else
      dy=d(ns)
      ns=ns-1
    end if
    y=y+dy
  end do

  return

  end subroutine polint
!------------------------------------------------------------------------------
!
! interpolation in 3 dimensions, follow zyx order
!
!------------------------------------------------------------------------------
  subroutine polin3(x1a,x2a,x3a,ya,x1,x2,x3,y,dy,ordn)

  implicit none

!~~~~~~> Input parameters:
  integer,intent(in) :: ordn
  real*8, dimension(1:ordn), intent(in) :: x1a,x2a,x3a
  real*8, dimension(1:ordn,1:ordn,1:ordn), intent(in) :: ya
  real*8, intent(in) :: x1,x2,x3
  real*8, intent(out) :: y,dy

!~~~~~~> Other parameters:

  integer  :: i,j,m,n
  real*8, dimension(ordn,ordn) :: yatmp
  real*8, dimension(ordn) :: ymtmp
  real*8, dimension(ordn) :: yntmp
  real*8, dimension(ordn) :: yqtmp

  m=size(x1a)
  n=size(x2a)
  
  do i=1,m
   do j=1,n

    yqtmp=ya(i,j,:)
    call polint(x3a,yqtmp,x3,yatmp(i,j),dy,ordn)

   end do

    yntmp=yatmp(i,:)
    call polint(x2a,yntmp,x2,ymtmp(i),dy,ordn)

  end do

  call polint(x1a,ymtmp,x1,y,dy,ordn)

  return

  end subroutine polin3
!--------------------------------------------------------------------------------------
! calculate L2norm  
  subroutine l2normhelper(ex, X, Y, Z,xmin,ymin,zmin,xmax,ymax,zmax,&
                          f,f_out,gw)

  implicit none
!~~~~~~> Input parameters:
  integer,intent(in ):: ex(1:3)
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3)),xmin,ymin,zmin,xmax,ymax,zmax
  integer,intent(in)::gw
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8, intent(out) :: f_out
!~~~~~~> Other variables:

  real*8, parameter :: ZEO = 0.D0
  real*8            :: dX, dY, dZ
  integer::imin,jmin,kmin
  integer::imax,jmax,kmax
  integer::i,j,k

  dX = X(2) - X(1)
  dY = Y(2) - Y(1)
  dZ = Z(2) - Z(1)

! for ghost zone
   imin = gw+1
   jmin = gw+1
   kmin = gw+1

   imax = ex(1) - gw
   jmax = ex(2) - gw
   kmax = ex(3) - gw

!for patch boundary (i.e., not ghost boundary)

if(dabs(X(ex(1))-xmax) < dX) imax = ex(1)
if(dabs(Y(ex(2))-ymax) < dY) jmax = ex(2)
if(dabs(Z(ex(3))-zmax) < dZ) kmax = ex(3)
if(dabs(X(1)-xmin) < dX) imin = 1
if(dabs(Y(1)-ymin) < dY) jmin = 1
if(dabs(Z(1)-zmin) < dZ) kmin = 1

f_out = sum(f(imin:imax,jmin:jmax,kmin:kmax)*f(imin:imax,jmin:jmax,kmin:kmax))

f_out = f_out*dX*dY*dZ

  return

  end subroutine l2normhelper
! locating the position of NaN
  subroutine ScaneNaN(ext,X,Y,Z,f)
  implicit none
  integer,dimension(3),intent(in) :: ext
  real*8,dimension(ext(1)) :: X
  real*8,dimension(ext(2)) :: Y
  real*8,dimension(ext(3)) :: Z
  real*8,dimension(ext(1),ext(2),ext(3)) :: f

  integer :: i,j,k

  do k=1,ext(3)
  do j=1,ext(2)
  do i=1,ext(1)
     if(abs(f(i,j,k)) .ne. abs(f(i,j,k))) write(*,*)X(i),Y(j),Z(k),f(i,j,k)
  enddo
  enddo
  enddo

  end subroutine ScaneNaN
! fortran version writefile
  subroutine fwritefile(time,nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,NN,filename,data_out)
  implicit none

  real*8,intent(in) :: time,xmin,xmax,ymin,ymax,zmin,zmax
  integer,intent(in) :: nx,ny,nz,NN
  real*8,dimension(nx,ny,nz),intent(in) :: data_out
  Character(Len=NN) :: filename

 
!  Open( 12 , File = filename,form='BINARY', access="SEQUENTIAL",status="replace",action='Write')
  Open( 12 , File = filename,form='UNFORMATTED', access="DIRECT",status="replace",action='Write')
  Write( 12 ) time,nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,data_out

  Close( 12 )

  end subroutine fwritefile
!--------------------------------------------------------------------------
!
! average for interface 
!
!--------------------------------------------------------------------------
  subroutine average(ext,f1,f2,fout)
  implicit none
  integer,dimension(3),   intent(in) :: ext
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in):: f1,f2
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out):: fout

  real*8,parameter::HLF=0.5d0

  fout = HLF*(f1+f2)

  return

  end subroutine average
!-----------------------------------------------------------------------------  
  subroutine average3(ext,f1,f2,fout)
  implicit none
  integer,dimension(3),   intent(in) :: ext
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in):: f1,f2
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out):: fout
! f1 ----------                ^
!                fout ------p  | t
! f2 ----------                |
! 2 points, 1st order interpolation
! 1   2
! f2  f1
! *---*--> t
!    ^
! f=3/4*f_1 + 1/4*f_2

  real*8,parameter::C1=0.75d0,C2=0.25d0

  fout = C1*f1+C2*f2

  return

  end subroutine average3
!-----------------------------------------------------------------------------  
  subroutine average2(ext,f1,f2,f3,fout)
  implicit none
  integer,dimension(3),   intent(in) :: ext
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in):: f1,f2,f3
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out):: fout
! f1 ----------                ^
!                fout ------   |
! f2 ----------                | t
!                              |
! f3 ----------                |
! 3 points, 2nd order interpolation
! 1   2   3
! f3  f2  f1
! *---*---*--> t
!       ^
! f=3/8*f_1 + 3/4*f_2 - 1/8*f_3

  real*8,parameter::C1=3.d0/8.d0,C2=3.d0/4.d0,C3=-1.d0/8.d0

  fout = C1*f1+C2*f2+C3*f3

  return

  end subroutine average2
!-----------------------------------------------------------------------------  
  subroutine average2p(ext,f1,f2,f3,fout)
  implicit none
  integer,dimension(3),   intent(in) :: ext
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in):: f1,f2,f3
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out):: fout
! f1 ----------                ^
!                fout ------p  |
! f2 ----------                | t
!                              |
! f3 ----------                |
! 3 points, 2nd order interpolation
! 1   2   3
! f3  f2  f1
! *---*---*--> t
!        ^
! f=21/32*f_1 + 7/16*f_2 - 3/32*f_3

  real*8,parameter::C1=5.d0/3.2d1,C2=1.5d1/1.6d1,C3=-3.d0/3.2d1

  fout = C1*f1+C2*f2+C3*f3

  return

  end subroutine average2p
!-----------------------------------------------------------------------------  
  subroutine average2m(ext,f1,f2,f3,fout)
  implicit none
  integer,dimension(3),   intent(in) :: ext
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in):: f1,f2,f3
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out):: fout
! f1 ----------                ^
!                fout ------m  |
! f2 ----------                | t
!                              |
! f3 ----------                |
! 3 points, 2nd order interpolation
! 1   2   3
! f3  f2  f1
! *---*---*--> t
!      ^
! f=5/32*f_1 + 15/16*f_2 - 3/32*f_3

  real*8,parameter::C1=5.d0/3.2d1,C2=1.5d1/1.6d1,C3=-3.d0/3.2d1

  fout = C1*f1+C2*f2+C3*f3

  return

  end subroutine average2m
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine lowerboundset(ex,chi0,TINNY)
  implicit none

!~~~~~~% Input parameters:

  integer ,intent(in):: ex(1:3)
  real*8  ,intent(in):: TINNY
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) ::chi0

  where(chi0 < TINNY) chi0 = TINNY

  return   

  end subroutine lowerboundset
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!global interpolation with given index and coeffients
  subroutine global_interpind(ex,X,Y,Z,f,f_int,x1,y1,z1,ORDN,SoA,symmetry,inds,coef,sst)
  implicit none

!~~~~~~> Input parameters:

  integer,                             intent(in) :: ex(1:3), symmetry,ORDN,sst
  real*8,intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8,                              intent(out):: f_int
  real*8,                              intent(in) :: x1,y1,z1
  real*8, dimension(3),                intent(in) :: SoA
  integer,dimension(3),                intent(in) :: inds
  real*8, dimension(3*ORDN),           intent(in) :: coef

!~~~~~~> Other parameters:

  real*8, dimension(-ORDN+1:ex(1)+ORDN,-ORDN+1:ex(2)+ORDN,-ORDN+1:ex(3)+ORDN) :: fh
  integer :: m
  integer,dimension(3) :: cxB,cxT
  real*8, dimension(ORDN,ORDN,ORDN) :: ya
  real*8, dimension(ORDN,ORDN) :: tmp2
  real*8, dimension(ORDN) :: tmp1
  real*8, dimension(3) :: SoAh

! +1 because c++ gives 0 for first point
  cxB = inds+1  
  cxT = cxB + ORDN - 1
      
  if(all(cxB>0).and.all(cxT<ex+1))then
     ya=f(cxB(1):cxT(1),cxB(2):cxT(2),cxB(3):cxT(3))
  elseif(any(cxB<-ORDN+1).or.any(cxT>ex+ORDN))then
     write(*,*)"error in global_interpind, cxB = ",cxB
     write(*,*)"                           cxT = ",cxT 
     write(*,*)"                           ext = ",ex
     stop
  else
     if(sst==-1)then
        SoAh = SoA
        if(any(cxT>ex)) write(*,*)"error global_interpind sst =",sst
     elseif(sst==0.or.sst==1)then
        SoAh = SoA
        SoAh(3) = 0
        if(cxB(3)<1.or.cxT(3)>ex(3)) write(*,*)"error global_interpind sst =",sst
     elseif(sst==2.or.sst==3)then
        SoAh(1) = SoA(2)
        SoAh(2) = SoA(3)
        SoAh(3) = 0
        if(cxB(3)<1.or.cxT(3)>ex(3)) write(*,*)"error global_interpind sst =",sst
     elseif(sst==4.or.sst==5)then
        SoAh(1) = SoA(1)
        SoAh(2) = SoA(3)
        SoAh(3) = 0
        if(cxB(3)<1.or.cxT(3)>ex(3)) write(*,*)"error global_interpind sst =",sst,cxB(3),cxT(3)
     endif
     call symmetry_tbd(ORDN,ex,f,fh,SoAh)
     ya=fh(cxB(1):cxT(1),cxB(2):cxT(2),cxB(3):cxT(3))
  endif 

  tmp2=0
  do m=1,ORDN
    tmp2 = tmp2 + coef(2*ORDN+m)*ya(:,:,m)
  enddo

  tmp1=0
  do m=1,ORDN
    tmp1 = tmp1 + coef(ORDN+m)*tmp2(:,m)
  enddo

  f_int=0
  do m=1,ORDN
    f_int = f_int + coef(m)*tmp1(m)
  enddo

  return

  end subroutine global_interpind
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!global interpolation with given index and coeffients
! special for shell to shell
  subroutine global_interpind2d(ex,X,Y,Z,f,f_int,x1,y1,z1,ORDN,SoA,symmetry,inds,coef,sst)
  implicit none

!~~~~~~> Input parameters:

  integer,                             intent(in) :: ex(1:3), symmetry,ORDN,sst
  real*8,intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8,                              intent(out):: f_int
  real*8,                              intent(in) :: x1,y1,z1
  real*8, dimension(3),                intent(in) :: SoA
  integer,dimension(3),                intent(in) :: inds
  real*8, dimension(2*ORDN),           intent(in) :: coef

!~~~~~~> Other parameters:

  real*8, dimension(-ORDN+1:ex(1)+ORDN,-ORDN+1:ex(2)+ORDN,ex(3)) :: fh
  integer :: m
  integer,dimension(2) :: cxB,cxT
  real*8, dimension(ORDN,ORDN) :: ya
  real*8, dimension(ORDN) :: tmp1
  real*8, dimension(2) :: SoAh

! +1 because c++ gives 0 for first point
  cxB = inds(1:2)+1  
  cxT = cxB + ORDN - 1
      
  if(all(cxB>0).and.all(cxT<ex(1:2)+1))then
     ya=f(cxB(1):cxT(1),cxB(2):cxT(2),inds(3))
  elseif(any(cxB<-ORDN+1).or.any(cxT>ex(1:2)+ORDN))then
     write(*,*)"error in global_interpind2d, cxB = ",cxB
     write(*,*)"                             cxT = ",cxT 
     write(*,*)"                             ext = ",ex(1:2)
     stop
  else
     if(sst==-1)then
       write(*,*)"error in global_interpind2d, sst = ",sst
       stop
     elseif(sst==0.or.sst==1)then
        SoAh = SoA(1:2)
     elseif(sst==2.or.sst==3)then
        SoAh(1) = SoA(2)
        SoAh(2) = SoA(3)
     elseif(sst==4.or.sst==5)then
        SoAh(1) = SoA(1)
        SoAh(2) = SoA(3)
     endif
     call symmetry_stbd(ORDN,ex,f,fh,SoAh)
     ya=fh(cxB(1):cxT(1),cxB(2):cxT(2),inds(3))
  endif 

  tmp1=0
  do m=1,ORDN
    tmp1 = tmp1 + coef(ORDN+m)*ya(:,m)
  enddo

  f_int=0
  do m=1,ORDN
    f_int = f_int + coef(m)*tmp1(m)
  enddo

  return

  end subroutine global_interpind2d
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!global interpolation with given index and coeffients
! special for shell to shell
  subroutine global_interpind1d(ex,X,Y,Z,f,f_int,x1,y1,z1,ORDN,SoA,symmetry,indsi,coef,sst,dumyd)
  implicit none

!~~~~~~> Input parameters:

  integer,                             intent(in) :: ex(1:3),symmetry,ORDN,sst,dumyd
  real*8,intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8,                              intent(out):: f_int
  real*8,                              intent(in) :: x1,y1,z1
  real*8, dimension(3),                intent(in) :: SoA
  integer,dimension(3),                intent(in) :: indsi
  real*8, dimension(ORDN),             intent(in) :: coef

!~~~~~~> Other parameters:

  real*8, dimension(-ORDN+1:ex(1)+ORDN,-ORDN+1:ex(2)+ORDN,ex(3)) :: fh
  integer :: m
  integer :: cxB,cxT
  real*8, dimension(ORDN) :: ya
  real*8 :: SoAh
  integer,dimension(3) :: inds

! +1 because c++ gives 0 for first point
  inds = indsi + 1
  cxB = inds(1)  
  cxT = cxB + ORDN - 1
      
  if(dumyd==1)then

  if(cxB>0.and.cxT<ex(1)+1)then
     ya=f(cxB:cxT,inds(2),inds(3))
  elseif(cxB<-ORDN+1.or.cxT>ex(1)+ORDN)then
     write(*,*)"error in global_interpind1d, cxB = ",cxB
     write(*,*)"                             cxT = ",cxT 
     write(*,*)"                             ext = ",ex(1)
     stop
  else
     if(sst==-1)then
       write(*,*)"error in global_interpind1d, sst = ",sst
       stop
     elseif(sst==0.or.sst==1)then
        SoAh = SoA(1)
     elseif(sst==2.or.sst==3)then
        SoAh = SoA(2)
     elseif(sst==4.or.sst==5)then
        SoAh = SoA(1)
     endif
     call symmetry_sntbd(ORDN,ex,f,fh,SoAh,1-dumyd)
     ya=fh(cxB:cxT,inds(2),inds(3))
  endif 

  elseif(dumyd==0)then

  if(cxB>0.and.cxT<ex(2)+1)then
     ya=f(inds(2),cxB:cxT,inds(3))
  elseif(cxB<-ORDN+1.or.cxT>ex(2)+ORDN)then
     write(*,*)"error in global_interpind1d, cxB = ",cxB
     write(*,*)"                             cxT = ",cxT 
     write(*,*)"                             ext = ",ex(2)
     stop
  else
     if(sst==-1)then
       write(*,*)"error in global_interpind1d, sst = ",sst
       stop
     elseif(sst==0.or.sst==1)then
        SoAh = SoA(2)
     elseif(sst==2.or.sst==3)then
        SoAh = SoA(1)
     elseif(sst==4.or.sst==5)then
        SoAh = SoA(2)
     endif
     call symmetry_sntbd(ORDN,ex,f,fh,SoAh,1-dumyd)
     ya=fh(inds(2),cxB:cxT,inds(3))
  endif 

  else
          write(*,*)"error in global_interpind1d, not recognized dumyd = ",dumyd
  endif

  f_int=0
  do m=1,ORDN
    f_int = f_int + coef(m)*ya(m)
  enddo

  return

  end subroutine global_interpind1d
!-----------------------------------------------------------------------------------------------------------------
! three dimensional interpolation for both vertex and cell center grid structure  
! for distinguishing shell and Cartesian
  subroutine global_interp_ss(ex,X,Y,Z,f,f_int,x1,y1,z1,ORDN,SoA,symmetry,sst)
  implicit none

!~~~~~~> Input parameters:

  integer,                             intent(in) :: ex(1:3), symmetry,ORDN,sst
  real*8,intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8,                              intent(out):: f_int
  real*8,                              intent(in) :: x1,y1,z1
  real*8, dimension(3),                intent(in) :: SoA

!~~~~~~> Other parameters:

  real*8, dimension(-ORDN+1:ex(1)+ORDN,-ORDN+1:ex(2)+ORDN,-ORDN+1:ex(3)+ORDN) :: fh
  real*8, dimension(3) :: SoAh
  integer :: j,m,imin,jmin,kmin
  integer,dimension(3) :: cxB,cxT,cxI,cmin,cmax
  real*8,dimension(3) :: cx
  real*8, dimension(1:ORDN) :: x1a
  real*8, dimension(1:ORDN,1:ORDN,1:ORDN) :: ya
  integer, parameter :: NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2
  real*8 :: dX,dY,dZ,ddy
  real*8, parameter :: ONE=1.d0

  imin = lbound(f,1)
  jmin = lbound(f,2)
  kmin = lbound(f,3)

  dX = X(imin+1)-X(imin)
  dY = Y(jmin+1)-Y(jmin)
  dZ = Z(kmin+1)-Z(kmin)

  forall( j = 1:ordn ) x1a(j) = ( j - 1 )* ONE

  cxI(1) = idint((x1-X(1))/dX+0.4)+1
  cxI(2) = idint((y1-Y(1))/dY+0.4)+1
  cxI(3) = idint((z1-Z(1))/dZ+0.4)+1

  cxB = cxI - ORDN/2+1
  cxT = cxB + ORDN - 1

  cmin = 1
  cmax = ex

  if(sst==-1)then
     SoAh = SoA
     cmin = -ORDN+1
  elseif(sst==0.or.sst==1)then
     SoAh = SoA
     SoAh(3) = 0
     cmin(1:2) = -ORDN+1
     cmax(1:2) = ex(1:2)+ORDN
  elseif(sst==2.or.sst==3)then
     SoAh(1) = SoA(2)
     SoAh(2) = SoA(3)
     SoAh(3) = 0
     cmin(1:2) = -ORDN+1
     cmax(1:2) = ex(1:2)+ORDN
  elseif(sst==4.or.sst==5)then
     SoAh(1) = SoA(1)
     SoAh(2) = SoA(3)
     SoAh(3) = 0
     cmin(1:2) = -ORDN+1
     cmax(1:2) = ex(1:2)+ORDN
  endif
  do m =1,3
   if(cxB(m) < cmin(m))then
      cxB(m) = cmin(m)
      cxT(m) = cxB(m) + ORDN - 1
   endif
   if(cxT(m) > cmax(m))then
      cxT(m) = cmax(m)
      cxB(m) = cxT(m) + 1 - ORDN
   endif
 enddo
 cx(1) = (x1 - X(1))/dX-cxB(1)+1
 cx(2) = (y1 - Y(1))/dY-cxB(2)+1
 cx(3) = (z1 - Z(1))/dZ-cxB(3)+1

  call symmetry_tbd(ORDN,ex,f,fh,SoAh)
  ya=fh(cxB(1):cxT(1),cxB(2):cxT(2),cxB(3):cxT(3))

  call polin3(x1a,x1a,x1a,ya,cx(1),cx(2),cx(3),f_int,ddy,ORDN)

  return

  end subroutine global_interp_ss
