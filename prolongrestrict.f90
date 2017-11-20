!$Id: prolongrestrict.f90,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $

! Because of overlap determination, source region is always larger than target
! region

#include "microdef.fh"

#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
!--------------------------------------------------------------------------------------
!
! Restrict from finner grids to coarser grids ignore the boundary point
! this routine is valid for all orders finite difference
!
! 1   2   3   4
! *---*---*---*
!     ^
! COPY directly!
!--------------------------------------------------------------------------------------

  subroutine restrict3(wei,llbc,uubc,extc,func,&
                       llbf,uubf,extf,funf,&
                       llbr,uubr,SoA,Symmetry)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: wei
!                                       coarse    fine       coarse
  real*8,dimension(3),   intent(in) :: llbc,uubc,llbf,uubf,llbr,uubr
  integer,dimension(3),  intent(in) :: extc,extf
  real*8, dimension(extc(1),extc(2),extc(3)),intent(inout):: func
  real*8, dimension(extf(1),extf(2),extf(3)),intent(in):: funf
  real*8, dimension(1:3), intent(in) :: SoA
  integer,intent(in)::Symmetry

!~~~~~~> local variables

  real*8, dimension(1:3) :: base
  integer,dimension(3) :: lbc,ubc,lbf,ubf,lbr,ubr,lbrf,ubrf
  integer,dimension(3) :: cxB,cxT,cxI
  integer :: i,j,k
  real*8,  dimension(4,4,4) :: ya
  real*8, dimension(4,4) :: tmp2
  real*8, dimension(4) :: tmp1
  real*8, parameter :: C1=-1.d0/1.6d1,C2=9.d0/1.6d1

  integer::imini,imaxi,jmini,jmaxi,kmini,kmaxi
  integer::imino,imaxo,jmino,jmaxo,kmino,kmaxo
  real*8,dimension(3) :: CD,FD

  if(wei.ne.3)then
     write(*,*)"prolongrestrict.f90::restrict3: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif
! it's possible a iolated point for target but not for source
  FD = (uubf-llbf)/(extf-1)
  CD = 2*FD

!take care the mismatch of the two segments of grid
  do i=1,3
     if(llbc(i) <= llbf(i))then
        base(i) = llbc(i)
     else
        j=idint((llbc(i)-llbf(i))/FD(i)+0.4)
        if(j/2*2 == j)then
           base(i) = llbf(i)
        else
           base(i) = llbf(i) - CD(i)/2
        endif
     endif
  enddo
!!! function idint:
!If A is of type REAL and |A| < 1, INT(A) equals 0. If |A| \geq 1, 
!then INT(A) equals the largest integer that does not exceed the range of A 
!and whose sign is the same as the sign of A.

    lbf = idint((llbf-base)/FD+0.4)+1
    ubf = idint((uubf-base)/FD+0.4)+1
    lbc = idint((llbc-base)/CD+0.4)+1
    ubc = idint((uubc-base)/CD+0.4)+1
    lbr = idint((llbr-base)/CD+0.4)+1
    lbrf = idint((llbr-base)/FD+0.4)+1
    ubr = idint((uubr-base)/CD+0.4)+1
    ubrf = idint((uubr-base)/FD+0.4)+1

!sanity check
  imino=lbr(1)-lbc(1) + 1
  imaxo=ubr(1)-lbc(1) + 1
  jmino=lbr(2)-lbc(2) + 1
  jmaxo=ubr(2)-lbc(2) + 1
  kmino=lbr(3)-lbc(3) + 1
  kmaxo=ubr(3)-lbc(3) + 1

  imini=lbrf(1)-lbf(1) + 1
  imaxi=ubrf(1)-lbf(1) + 1
  jmini=lbrf(2)-lbf(2) + 1
  jmaxi=ubrf(2)-lbf(2) + 1
  kmini=lbrf(3)-lbf(3) + 1
  kmaxi=ubrf(3)-lbf(3) + 1

  if(imino.lt.1.or.jmino.lt.1.or.kmino.lt.1.or.&
     imini.lt.1.or.jmini.lt.1.or.kmini.lt.1.or.&
     imaxo.gt.extc(1).or.jmaxo.gt.extc(2).or.kmaxo.gt.extc(3).or.&
     imaxi.gt.extf(1).or.jmaxi.gt.extf(2).or.kmaxi.gt.extf(3))then
          write(*,*)"error in restrict for"
          write(*,*)"mino = ",imino,jmino,kmino
          write(*,*)"maxo = ",imaxo,jmaxo,kmaxo
          write(*,*)"extc = ",extc
          write(*,*)"CD = ",CD
          write(*,*)"mini = ",imini,jmini,kmini
          write(*,*)"maxi = ",imaxi,jmaxi,kmaxi
          write(*,*)"extf = ",extf
          write(*,*)"FD = ",FD
          write(*,*)"from"
          write(*,*)lbf,ubf,extf
          write(*,*)"to"
          write(*,*)lbc,ubc,extc
          write(*,*)"want"
          write(*,*)lbr,ubr,lbrf,ubrf
          write(*,*)"llbf = ",llbf
          write(*,*)"uubf = ",uubf
          write(*,*)"llbc = ",llbc
          write(*,*)"uubc = ",uubc
          write(*,*)"llbr = ",llbr
          write(*,*)"uubr = ",uubr
          stop
  endif

!~~~~~~> restriction start...
  do k = kmino,kmaxo
   do j = jmino,jmaxo
    do i = imino,imaxo

       cxI(1) = i
       cxI(2) = j
       cxI(3) = k
! change to fine level reference
!|*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*| 
!|x===============x===============x===============x========|
       cxI = 2*(cxI+lbc-1) - 1
! change to array index      
       cxI = cxI - lbf + 1

       func(i,j,k)= funf(cxI(1),cxI(2),cxI(3))
    enddo
   enddo
  enddo
  
  return

  end subroutine restrict3
#endif  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! for different finite difference order usage
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if (ghost_width == 2)
! second order code
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
!--------------------------------------------------------------------------
!
! Prolongation from coarser grids to finer grids
! 4 points, 3rd order interpolation
! 1   2   3   4
! *---*---*---*
!       ^
! f=-1/16*f_1 + 9/16*f_2
!   -1/16*f_4 + 9/16*f_3
!--------------------------------------------------------------------------

  subroutine prolong3(wei,llbc,uubc,extc,func,&
                      llbf,uubf,extf,funf,&
                      llbp,uubp,SoA,Symmetry)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: wei
!                                       coarse    fine       coarse
  real*8,dimension(3),   intent(in) :: llbc,uubc,llbf,uubf,llbp,uubp
  integer,dimension(3),  intent(in) :: extc,extf
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in)   :: func
  real*8, dimension(extf(1),extf(2),extf(3)),intent(inout):: funf
  real*8, dimension(1:3), intent(in) :: SoA
  integer,intent(in)::Symmetry

!~~~~~~> local variables

  real*8, dimension(1:3) :: base
  integer,dimension(3) :: lbc,ubc,lbf,ubf,lbp,ubp,lbpc,ubpc
  integer,dimension(3) :: cxB,cxT,cxI
  integer :: i,j,k,ii,jj,kk
  real*8,  dimension(4,4,4) :: ya
  real*8, dimension(4,4) :: tmp2
  real*8, dimension(4) :: tmp1

  real*8, parameter :: C1=-1.d0/16,C2=9.d0/16
  real*8, parameter :: C4=C1,C3=C2

  integer::imini,imaxi,jmini,jmaxi,kmini,kmaxi
  integer::imino,imaxo,jmino,jmaxo,kmino,kmaxo
  logical::decide3d

  real*8,dimension(3) :: CD,FD

  if(wei.ne.3)then
     write(*,*)"prolongrestrict.f90::prolong3: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

! it's possible a iolated point for target but not for source
  CD = (uubc-llbc)/(extc-1)
  FD = CD/2

!take care the mismatch of the two segments of grid
  do i=1,3
     if(llbc(i) <= llbf(i))then
        base(i) = llbc(i)
     else
        j=idint((llbc(i)-llbf(i))/FD(i)+0.4)
        if(j/2*2 == j)then
           base(i) = llbf(i)
        else
           base(i) = llbf(i) - CD(i)/2
        endif
     endif
  enddo

!!! function idint:
!If A is of type REAL and |A| < 1, INT(A) equals 0. If |A| \geq 1, 
!then INT(A) equals the largest integer that does not exceed the range of A 
!and whose sign is the same as the sign of A.

    lbf = idint((llbf-base)/FD+0.4)+1
    ubf = idint((uubf-base)/FD+0.4)+1
    lbc = idint((llbc-base)/CD+0.4)+1
    ubc = idint((uubc-base)/CD+0.4)+1
    lbp = idint((llbp-base)/FD+0.4)+1
    lbpc = idint((llbp-base)/CD+0.4)+1
    ubp = idint((uubp-base)/FD+0.4)+1
    ubpc = idint((uubp-base)/CD+0.4)+1
!sanity check
  imino=lbp(1)-lbf(1) + 1
  imaxo=ubp(1)-lbf(1) + 1
  jmino=lbp(2)-lbf(2) + 1
  jmaxo=ubp(2)-lbf(2) + 1
  kmino=lbp(3)-lbf(3) + 1
  kmaxo=ubp(3)-lbf(3) + 1

  imini=lbpc(1)-lbc(1) + 1
  imaxi=ubpc(1)-lbc(1) + 1
  jmini=lbpc(2)-lbc(2) + 1
  jmaxi=ubpc(2)-lbc(2) + 1
  kmini=lbpc(3)-lbc(3) + 1
  kmaxi=ubpc(3)-lbc(3) + 1

  if(imino.lt.1.or.jmino.lt.1.or.kmino.lt.1.or.&
     imini.lt.1.or.jmini.lt.1.or.kmini.lt.1.or.&
     imaxo.gt.extf(1).or.jmaxo.gt.extf(2).or.kmaxo.gt.extf(3).or.&
     imaxi.gt.extc(1)-1.or.jmaxi.gt.extc(2)-1.or.kmaxi.gt.extc(3)-1)then
          write(*,*)"error in prolongation for"
          write(*,*)"from"
          write(*,*)llbc,uubc
          write(*,*)lbc,ubc
          write(*,*)"to"
          write(*,*)llbf,uubf
          write(*,*)lbf,ubf
          write(*,*)"want"
          write(*,*)llbp,uubp
          write(*,*)lbp,ubp,lbpc,ubpc
          if(imini.lt.1) write(*,*)"imini = ",imini
          if(jmini.lt.1) write(*,*)"jmini = ",jmini       
          if(kmini.lt.1) write(*,*)"kmini = ",kmini      
          if(imino.lt.1) write(*,*)"imino = ",imino       
          if(jmino.lt.1) write(*,*)"jmino = ",jmino       
          if(kmino.lt.1) write(*,*)"kmino = ",kmino   
          if(imaxi.gt.extc(1)) write(*,*)"imaxi = ",imaxi,"extc(1) = ",extc(1) 
          if(jmaxi.gt.extc(2)) write(*,*)"jmaxi = ",jmaxi,"extc(2) = ",extc(2) 
          if(kmaxi.gt.extc(3)) write(*,*)"kmaxi = ",kmaxi,"extc(3) = ",extc(3) 
          if(imaxo.gt.extf(1)) write(*,*)"imaxo = ",imaxo,"extf(1) = ",extf(1) 
          if(jmaxo.gt.extf(2)) write(*,*)"jmaxo = ",jmaxo,"extf(2) = ",extf(2) 
          if(kmaxo.gt.extf(3)) write(*,*)"kmaxo = ",kmaxo,"extf(3) = ",extf(3) 
          return
  endif
!~~~~~~> prolongation start...
  do k = kmino,kmaxo
   do j = jmino,jmaxo
    do i = imino,imaxo
! change to coarse level reference
!|*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*| 
!|x===============x===============x===============x========|
!       if(i/2*2 == i)then
!         cxI(1) = (i+lbf(1)-1)/2
!       else
!         cxI(1) = (i+lbf(1)-1)/2+1
!       endif      
!       if(j/2*2 == j)then
!         cxI(2) = (j+lbf(2)-1)/2
!       else
!         cxI(2) = (j+lbf(2)-1)/2+1
!       endif
!       if(k/2*2 == k)then
!         cxI(3) = (k+lbf(3)-1)/2
!       else
!         cxI(3) = (k+lbf(3)-1)/2+1
!       endif
! above code segment is equivalent to 
       cxI(1) = i
       cxI(2) = j
       cxI(3) = k
       cxI = (cxI+lbf)/2
! change to array index      
       cxI = cxI - lbc + 1

       ii=i+lbf(1)-1
       jj=j+lbf(2)-1
       kk=k+lbf(3)-1

       if(any(cxI+2 > extc)) write(*,*)"error in prolong"
       if(ii/2*2==ii)then
         if(jj/2*2==jj)then
           if(kk/2*2==kk)then
!  due to ghost zone, we can deal with symmetry boundary like this                   
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= C1*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)-1)+&
                      C2*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)  )+&
                      C3*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+1)+&
                      C4*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+2)
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= C1*ya(:,:,1)+C2*ya(:,:,2)+C3*ya(:,:,3)+C4*ya(:,:,4)
             endif
             tmp1= C1*tmp2(:,1)+C2*tmp2(:,2)+C3*tmp2(:,3)+C4*tmp2(:,4)
             funf(i,j,k)= C1*tmp1(1)+C2*tmp1(2)+C3*tmp1(3)+C4*tmp1(4)
           else
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3))
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= ya(:,:,2)
             endif
             tmp1= C1*tmp2(:,1)+C2*tmp2(:,2)+C3*tmp2(:,3)+C4*tmp2(:,4)
             funf(i,j,k)=  C1*tmp1(1)+C2*tmp1(2)+C3*tmp1(3)+C4*tmp1(4)
           endif
         else
           if(kk/2*2==kk)then
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= C1*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)-1)+&
                      C2*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)  )+&
                      C3*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+1)+&
                      C4*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+2)
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= C1*ya(:,:,1)+C2*ya(:,:,2)+C3*ya(:,:,3)+C4*ya(:,:,4)
             endif
             tmp1= tmp2(:,2)
             funf(i,j,k)= C1*tmp1(1)+C2*tmp1(2)+C3*tmp1(3)+C4*tmp1(4)
           else
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3))
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= ya(:,:,2)
             endif
             tmp1= tmp2(:,2)
             funf(i,j,k)=  C1*tmp1(1)+C2*tmp1(2)+C3*tmp1(3)+C4*tmp1(4)
           endif
         endif
       else
         if(jj/2*2==jj)then
           if(kk/2*2==kk)then               
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= C1*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)-1)+&
                      C2*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)  )+&
                      C3*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+1)+&
                      C4*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+2)
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= C1*ya(:,:,1)+C2*ya(:,:,2)+C3*ya(:,:,3)+C4*ya(:,:,4)
             endif
             tmp1= C1*tmp2(:,1)+C2*tmp2(:,2)+C3*tmp2(:,3)+C4*tmp2(:,4)
             funf(i,j,k)= tmp1(2)
           else
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3))
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= ya(:,:,2)
             endif
             tmp1= C1*tmp2(:,1)+C2*tmp2(:,2)+C3*tmp2(:,3)+C4*tmp2(:,4)
             funf(i,j,k)= tmp1(2)
           endif
         else
           if(kk/2*2==kk)then
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= C1*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)-1)+&
                      C2*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)  )+&
                      C3*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+1)+&
                      C4*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+2)
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= C1*ya(:,:,1)+C2*ya(:,:,2)+C3*ya(:,:,3)+C4*ya(:,:,4)
             endif
             tmp1= tmp2(:,2)
             funf(i,j,k)= tmp1(2)
           else
             funf(i,j,k)= func(cxI(1),cxI(2),cxI(3))
           endif
         endif
       endif
    enddo
   enddo
  enddo

  return

  end subroutine prolong3

#else
#ifdef Cell

!--------------------------------------------------------------------------------------
!
! Restrict from finner grids to coarser grids ignore the boundary point
!
! 4 points, 3rd order interpolation
! 1   2   3   4
! *---*---*---*
!       ^
! f=-1/16*(f_1+f_4) + 9/16*(f_2+f_3)
!--------------------------------------------------------------------------------------

  subroutine restrict3(wei,llbc,uubc,extc,func,&
                       llbf,uubf,extf,funf,&
                       llbr,uubr,SoA,Symmetry)
  implicit none

!~~~~~~> input arguments
  integer,intent(in)::wei
!                                       coarse    fine       coarse
  real*8,dimension(3),   intent(in) :: llbc,uubc,llbf,uubf,llbr,uubr
  integer,dimension(3),  intent(in) :: extc,extf
  real*8, dimension(extc(1),extc(2),extc(3)),intent(inout):: func
  real*8, dimension(extf(1),extf(2),extf(3)),intent(in):: funf
  real*8, dimension(1:3), intent(in) :: SoA
  integer,intent(in)::Symmetry

!~~~~~~> local variables

  real*8, dimension(1:3) :: base
  integer,dimension(3) :: lbc,ubc,lbf,ubf,lbr,ubr,lbrf,ubrf
  integer,dimension(3) :: cxB,cxT,cxI
  integer :: i,j,k
  real*8,  dimension(4,4,4) :: ya
  real*8, dimension(4,4) :: tmp2
  real*8, dimension(4) :: tmp1
  real*8, parameter :: C1=-1.d0/1.6d1,C2=9.d0/1.6d1

  integer::imini,imaxi,jmini,jmaxi,kmini,kmaxi
  integer::imino,imaxo,jmino,jmaxo,kmino,kmaxo
  logical::decide3d
  real*8,dimension(3) :: CD,FD
  
  if(wei.ne.3)then
     write(*,*)"prolongrestrict.f90::restrict3: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

  CD = (uubc-llbc)/extc
  FD = (uubf-llbf)/extf

!take care the mismatch of the two segments of grid
  do i=1,3
     if(llbc(i) <= llbf(i))then
        base(i) = llbc(i)
     else
        j=idint((llbc(i)-llbf(i))/FD(i)+0.4)
        if(j/2*2 == j)then
           base(i) = llbf(i)
        else
           base(i) = llbf(i) - CD(i)/2
        endif
     endif
  enddo
!!! function idint:
!If A is of type REAL and |A| < 1, INT(A) equals 0. If |A| \geq 1, 
!then INT(A) equals the largest integer that does not exceed the range of A 
!and whose sign is the same as the sign of A.

! note say base = 0, llbf = 0, uubf = 2
! llbf->1 and uubf->2
    lbf = idint((llbf-base)/FD+0.4)+1
    ubf = idint((uubf-base)/FD+0.4)
    lbc = idint((llbc-base)/CD+0.4)+1
    ubc = idint((uubc-base)/CD+0.4)
    lbr = idint((llbr-base)/CD+0.4)+1
    lbrf = idint((llbr-base)/FD+0.4)+1
    ubr = idint((uubr-base)/CD+0.4)
    ubrf = idint((uubr-base)/FD+0.4)

!sanity check
  imino=lbr(1)-lbc(1) + 1
  imaxo=ubr(1)-lbc(1) + 1
  jmino=lbr(2)-lbc(2) + 1
  jmaxo=ubr(2)-lbc(2) + 1
  kmino=lbr(3)-lbc(3) + 1
  kmaxo=ubr(3)-lbc(3) + 1

  imini=lbrf(1)-lbf(1) + 1
  imaxi=ubrf(1)-lbf(1) + 1
  jmini=lbrf(2)-lbf(2) + 1
  jmaxi=ubrf(2)-lbf(2) + 1
  kmini=lbrf(3)-lbf(3) + 1
  kmaxi=ubrf(3)-lbf(3) + 1

  if(imino.lt.1.or.jmino.lt.1.or.kmino.lt.1.or.&
     imini.lt.1.or.jmini.lt.1.or.kmini.lt.1.or.&
     imaxo.gt.extc(1).or.jmaxo.gt.extc(2).or.kmaxo.gt.extc(3).or.&
     imaxi.gt.extf(1)-1.or.jmaxi.gt.extf(2)-1.or.kmaxi.gt.extf(3)-1)then
          write(*,*)"error in restrict for"
          write(*,*)"from"
          write(*,*)lbf,ubf
          write(*,*)"to"
          write(*,*)lbc,ubc
          write(*,*)"want"
          write(*,*)lbr,ubr,lbrf,ubrf
          write(*,*)"llbf = ",llbf
          write(*,*)"uubf = ",uubf
          write(*,*)"llbc = ",llbc
          write(*,*)"uubc = ",uubc
          write(*,*)"llbr = ",llbr
          write(*,*)"uubr = ",uubr
          write(*,*)"base = ",base
          stop
  endif

!~~~~~~> restriction start...
  do k = kmino,kmaxo
   do j = jmino,jmaxo
    do i = imino,imaxo

       cxI(1) = i
       cxI(2) = j
       cxI(3) = k
! change to fine level reference
!|---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*---| 
!|=======x===============x===============x===============x=======|
       cxI = 2*(cxI+lbc-1) - 1
! change to array index      
       cxI = cxI - lbf + 1

       if(any(cxI+2 > extf)) write(*,*)"error in restrict"
!  due to ghost zone, we can deal with symmetry boundary like this                   
       if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
          tmp2= C1*(funf(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)-1)+funf(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+2))&
               +C2*(funf(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)  )+funf(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+1))
       else
          cxB=cxI-1
          cxT=cxI+2
          if(decide3d(extf,funf,funf,cxB,cxT,SoA,ya,4,Symmetry))then
             write(*,*)"restrict3 position index: ",i+lbc(1)-1,j+lbc(2)-1,k+lbc(3)-1
             stop
          endif
          tmp2=C1*(ya(:,:,1)+ya(:,:,4))+C2*(ya(:,:,2)+ya(:,:,3))
       endif
       tmp1= C1*(tmp2(:,1)+tmp2(:,4))+C2*(tmp2(:,2)+tmp2(:,3))
       func(i,j,k)= C1*(tmp1(1)+tmp1(4))+C2*(tmp1(2)+tmp1(3))
    enddo
   enddo
  enddo
  
  return

  end subroutine restrict3
!--------------------------------------------------------------------------
!
! Prolongation from coarser grids to finer grids
! 4 points, 3rd order interpolation
! 1   2   3   4
! *---*---*---*
!      ^
! f=-7/128*f_1 + 105/128*f_2
!   -5/128*f_4 +  35/128*f_3
!--------------------------------------------------------------------------

  subroutine prolong3(wei,llbc,uubc,extc,func,&
                      llbf,uubf,extf,funf,&
                      llbp,uubp,SoA,Symmetry)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: wei
!                                       coarse    fine       coarse
  real*8,dimension(3),   intent(in) :: llbc,uubc,llbf,uubf,llbp,uubp
  integer,dimension(3),  intent(in) :: extc,extf
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in)   :: func
  real*8, dimension(extf(1),extf(2),extf(3)),intent(inout):: funf
  real*8, dimension(1:3), intent(in) :: SoA
  integer,intent(in)::Symmetry

!~~~~~~> local variables

  real*8, dimension(1:3) :: base
  integer,dimension(3) :: lbc,ubc,lbf,ubf,lbp,ubp,lbpc,ubpc
  integer,dimension(3) :: cxB,cxT,cxI
  integer :: i,j,k,ii,jj,kk
  real*8,  dimension(4,4,4) :: ya
  real*8, dimension(4,4) :: tmp2
  real*8, dimension(4) :: tmp1

  real*8, parameter :: C1=-7.d0/1.28d2,C2=1.05d2/1.28d2
  real*8, parameter :: C4=-5.d0/1.28d2,C3= 3.5d1/1.28d2

  integer::imini,imaxi,jmini,jmaxi,kmini,kmaxi
  integer::imino,imaxo,jmino,jmaxo,kmino,kmaxo
  logical::decide3d

  real*8,dimension(3) :: CD,FD
  
  if(wei.ne.3)then
     write(*,*)"prolongrestrict.f90::prolong3: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

  CD = (uubc-llbc)/extc
  FD = (uubf-llbf)/extf

!take care the mismatch of the two segments of grid
  do i=1,3
     if(llbc(i) <= llbf(i))then
        base(i) = llbc(i)
     else
        j=idint((llbc(i)-llbf(i))/FD(i)+0.4)
        if(j/2*2 == j)then
           base(i) = llbf(i)
        else
           base(i) = llbf(i) - CD(i)/2
        endif
     endif
  enddo

!!! function idint:
!If A is of type REAL and |A| < 1, INT(A) equals 0. If |A| \geq 1, 
!then INT(A) equals the largest integer that does not exceed the range of A 
!and whose sign is the same as the sign of A.

    lbf = idint((llbf-base)/FD+0.4)+1
    ubf = idint((uubf-base)/FD+0.4)
    lbc = idint((llbc-base)/CD+0.4)+1
    ubc = idint((uubc-base)/CD+0.4)
    lbp = idint((llbp-base)/FD+0.4)+1
    lbpc = idint((llbp-base)/CD+0.4)+1
    ubp = idint((uubp-base)/FD+0.4)
    ubpc = idint((uubp-base)/CD+0.4)

!sanity check
  imino=lbp(1)-lbf(1) + 1
  imaxo=ubp(1)-lbf(1) + 1
  jmino=lbp(2)-lbf(2) + 1
  jmaxo=ubp(2)-lbf(2) + 1
  kmino=lbp(3)-lbf(3) + 1
  kmaxo=ubp(3)-lbf(3) + 1

  imini=lbpc(1)-lbc(1) + 1
  imaxi=ubpc(1)-lbc(1) + 1
  jmini=lbpc(2)-lbc(2) + 1
  jmaxi=ubpc(2)-lbc(2) + 1
  kmini=lbpc(3)-lbc(3) + 1
  kmaxi=ubpc(3)-lbc(3) + 1

  if(imino.lt.1.or.jmino.lt.1.or.kmino.lt.1.or.&
     imini.lt.1.or.jmini.lt.1.or.kmini.lt.1.or.&
     imaxo.gt.extf(1).or.jmaxo.gt.extf(2).or.kmaxo.gt.extf(3).or.&
     imaxi.gt.extc(1)-1.or.jmaxi.gt.extc(2)-1.or.kmaxi.gt.extc(3)-1)then
          write(*,*)"error in prolongation for"
          write(*,*)"from"
          write(*,*)llbc,uubc
          write(*,*)lbc,ubc
          write(*,*)"to"
          write(*,*)llbf,uubf
          write(*,*)lbf,ubf
          write(*,*)"want"
          write(*,*)llbp,uubp
          write(*,*)lbp,ubp,lbpc,ubpc
          return
  endif
!~~~~~~> prolongation start...
  do k = kmino,kmaxo
   do j = jmino,jmaxo
    do i = imino,imaxo

       cxI(1) = i
       cxI(2) = j
       cxI(3) = k
! change to coarse level reference
!|---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*---| 
!|=======x===============x===============x===============x=======|
       cxI = (cxI+lbf-1)/2
! change to array index      
       cxI = cxI - lbc + 1

       ii=i+lbf(1)-1
       jj=j+lbf(2)-1
       kk=k+lbf(3)-1

       if(any(cxI+2 > extc)) write(*,*)"error in prolong"
       if(ii/2*2==ii)then
         if(jj/2*2==jj)then
           if(kk/2*2==kk)then
!  due to ghost zone, we can deal with symmetry boundary like this                   
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= C1*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)-1)+&
                      C2*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)  )+&
                      C3*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+1)+&
                      C4*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+2)
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= C1*ya(:,:,1)+C2*ya(:,:,2)+C3*ya(:,:,3)+C4*ya(:,:,4)
             endif
             tmp1= C1*tmp2(:,1)+C2*tmp2(:,2)+C3*tmp2(:,3)+C4*tmp2(:,4)
             funf(i,j,k)= C1*tmp1(1)+C2*tmp1(2)+C3*tmp1(3)+C4*tmp1(4)
           else
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= C4*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)-1)+&
                      C3*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)  )+&
                      C2*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+1)+&
                      C1*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+2)
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= C4*ya(:,:,1)+C3*ya(:,:,2)+C2*ya(:,:,3)+C1*ya(:,:,4)
             endif
             tmp1= C1*tmp2(:,1)+C2*tmp2(:,2)+C3*tmp2(:,3)+C4*tmp2(:,4)
             funf(i,j,k)=  C1*tmp1(1)+C2*tmp1(2)+C3*tmp1(3)+C4*tmp1(4)
           endif
         else
           if(kk/2*2==kk)then
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= C1*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)-1)+&
                      C2*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)  )+&
                      C3*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+1)+&
                      C4*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+2)
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= C1*ya(:,:,1)+C2*ya(:,:,2)+C3*ya(:,:,3)+C4*ya(:,:,4)
             endif
             tmp1= C4*tmp2(:,1)+C3*tmp2(:,2)+C2*tmp2(:,3)+C1*tmp2(:,4)
             funf(i,j,k)= C1*tmp1(1)+C2*tmp1(2)+C3*tmp1(3)+C4*tmp1(4)
           else
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= C4*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)-1)+&
                      C3*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)  )+&
                      C2*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+1)+&
                      C1*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+2)
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= C4*ya(:,:,1)+C3*ya(:,:,2)+C2*ya(:,:,3)+C1*ya(:,:,4)
             endif
             tmp1= C4*tmp2(:,1)+C3*tmp2(:,2)+C2*tmp2(:,3)+C1*tmp2(:,4)
             funf(i,j,k)=  C1*tmp1(1)+C2*tmp1(2)+C3*tmp1(3)+C4*tmp1(4)
           endif
         endif
       else
         if(jj/2*2==jj)then
           if(kk/2*2==kk)then               
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= C1*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)-1)+&
                      C2*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)  )+&
                      C3*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+1)+&
                      C4*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+2)
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= C1*ya(:,:,1)+C2*ya(:,:,2)+C3*ya(:,:,3)+C4*ya(:,:,4)
             endif
             tmp1= C1*tmp2(:,1)+C2*tmp2(:,2)+C3*tmp2(:,3)+C4*tmp2(:,4)
             funf(i,j,k)= C4*tmp1(1)+C3*tmp1(2)+C2*tmp1(3)+C1*tmp1(4)
           else
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= C4*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)-1)+&
                      C3*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)  )+&
                      C2*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+1)+&
                      C1*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+2)
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= C4*ya(:,:,1)+C3*ya(:,:,2)+C2*ya(:,:,3)+C1*ya(:,:,4)
             endif
             tmp1= C1*tmp2(:,1)+C2*tmp2(:,2)+C3*tmp2(:,3)+C4*tmp2(:,4)
             funf(i,j,k)=  C4*tmp1(1)+C3*tmp1(2)+C2*tmp1(3)+C1*tmp1(4)
           endif
         else
           if(kk/2*2==kk)then
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= C1*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)-1)+&
                      C2*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)  )+&
                      C3*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+1)+&
                      C4*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+2)
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= C1*ya(:,:,1)+C2*ya(:,:,2)+C3*ya(:,:,3)+C4*ya(:,:,4)
             endif
             tmp1= C4*tmp2(:,1)+C3*tmp2(:,2)+C2*tmp2(:,3)+C1*tmp2(:,4)
             funf(i,j,k)= C4*tmp1(1)+C3*tmp1(2)+C2*tmp1(3)+C1*tmp1(4)
           else
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= C4*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)-1)+&
                      C3*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)  )+&
                      C2*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+1)+&
                      C1*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+2)
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= C4*ya(:,:,1)+C3*ya(:,:,2)+C2*ya(:,:,3)+C1*ya(:,:,4)
             endif
             tmp1= C4*tmp2(:,1)+C3*tmp2(:,2)+C2*tmp2(:,3)+C1*tmp2(:,4)
             funf(i,j,k)=  C4*tmp1(1)+C3*tmp1(2)+C2*tmp1(3)+C1*tmp1(4)
           endif
         endif
       endif
    enddo
   enddo
  enddo

  return

  end subroutine prolong3

  subroutine prolong3new(wei,llbc,uubc,extc,func,&
                      llbf,uubf,extf,funf,&
                      llbp,uubp,SoA,Symmetry)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: wei
!                                       coarse    fine       coarse
  real*8,dimension(3),   intent(in) :: llbc,uubc,llbf,uubf,llbp,uubp
  integer,dimension(3),  intent(in) :: extc,extf
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in)   :: func
  real*8, dimension(extf(1),extf(2),extf(3)),intent(inout):: funf
  real*8, dimension(1:3), intent(in) :: SoA
  integer,intent(in)::Symmetry

!~~~~~~> local variables

  real*8, dimension(1:3) :: base
  integer,dimension(3) :: lbc,ubc,lbf,ubf,lbp,ubp,lbpc,ubpc
  integer,dimension(3) :: cxB,cxT,cxI
  integer :: i,j,k,ii,jj,kk
  real*8,  dimension(4,4,4) :: ya
  real*8, dimension(4,4) :: tmp2
  real*8, dimension(4) :: tmp1

  real*8, parameter :: C1=-7.d0/1.28d2,C2=1.05d2/1.28d2
  real*8, parameter :: C4=-5.d0/1.28d2,C3= 3.5d1/1.28d2

  integer::imini,imaxi,jmini,jmaxi,kmini,kmaxi
  integer::imino,imaxo,jmino,jmaxo,kmino,kmaxo
  logical::decide3d

  real*8,dimension(3) :: CD,FD
  real*8,dimension(3,4) :: CC
  
  if(wei.ne.3)then
     write(*,*)"prolongrestrict.f90::prolong3: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

  CD = (uubc-llbc)/extc
  FD = (uubf-llbf)/extf

!take care the mismatch of the two segments of grid
  do i=1,3
     if(llbc(i) <= llbf(i))then
        base(i) = llbc(i)
     else
        j=idint((llbc(i)-llbf(i))/FD(i)+0.4)
        if(j/2*2 == j)then
           base(i) = llbf(i)
        else
           base(i) = llbf(i) - CD(i)/2
        endif
     endif
  enddo

!!! function idint:
!If A is of type REAL and |A| < 1, INT(A) equals 0. If |A| \geq 1, 
!then INT(A) equals the largest integer that does not exceed the range of A 
!and whose sign is the same as the sign of A.

    lbf = idint((llbf-base)/FD+0.4)+1
    ubf = idint((uubf-base)/FD+0.4)
    lbc = idint((llbc-base)/CD+0.4)+1
    ubc = idint((uubc-base)/CD+0.4)
    lbp = idint((llbp-base)/FD+0.4)+1
    lbpc = idint((llbp-base)/CD+0.4)+1
    ubp = idint((uubp-base)/FD+0.4)
    ubpc = idint((uubp-base)/CD+0.4)

!sanity check
  imino=lbp(1)-lbf(1) + 1
  imaxo=ubp(1)-lbf(1) + 1
  jmino=lbp(2)-lbf(2) + 1
  jmaxo=ubp(2)-lbf(2) + 1
  kmino=lbp(3)-lbf(3) + 1
  kmaxo=ubp(3)-lbf(3) + 1

  imini=lbpc(1)-lbc(1) + 1
  imaxi=ubpc(1)-lbc(1) + 1
  jmini=lbpc(2)-lbc(2) + 1
  jmaxi=ubpc(2)-lbc(2) + 1
  kmini=lbpc(3)-lbc(3) + 1
  kmaxi=ubpc(3)-lbc(3) + 1

  if(imino.lt.1.or.jmino.lt.1.or.kmino.lt.1.or.&
     imini.lt.1.or.jmini.lt.1.or.kmini.lt.1.or.&
     imaxo.gt.extf(1).or.jmaxo.gt.extf(2).or.kmaxo.gt.extf(3).or.&
     imaxi.gt.extc(1)-2.or.jmaxi.gt.extc(2)-2.or.kmaxi.gt.extc(3)-2)then
          write(*,*)"error in prolongation for"
          write(*,*)"from"
          write(*,*)llbc,uubc
          write(*,*)lbc,ubc
          write(*,*)"to"
          write(*,*)llbf,uubf
          write(*,*)lbf,ubf
          write(*,*)"want"
          write(*,*)llbp,uubp
          write(*,*)lbp,ubp,lbpc,ubpc
          return
  endif
!~~~~~~> prolongation start...
  do i=1,3
     if(lbp(i)/2*2 == lbp(i))then
       CC(i,1) = C1
       CC(i,2) = C2
       CC(i,3) = C3
       CC(i,4) = C4
     else
       CC(i,1) = C4
       CC(i,2) = C3
       CC(i,3) = C2
       CC(i,4) = C1
     endif
  enddo

  do k = kmino,kmaxo,2
   do j = jmino,jmaxo,2
    do i = imino,imaxo,2

       cxI(1) = i
       cxI(2) = j
       cxI(3) = k
! change to coarse level reference
!|---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*---| 
!|=======x===============x===============x===============x=======|
       cxI = (cxI+lbf-1)/2
! change to array index      
       cxI = cxI - lbc + 1
!  due to ghost zone, we can deal with symmetry boundary like this                   
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= CC(3,1)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)-1)+&
                      CC(3,2)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)  )+&
                      CC(3,3)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+1)+&
                      CC(3,4)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+2)
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= CC(3,1)*ya(:,:,1)+CC(3,2)*ya(:,:,2)+CC(3,3)*ya(:,:,3)+CC(3,4)*ya(:,:,4)
             endif
             tmp1= CC(2,1)*tmp2(:,1)+CC(2,2)*tmp2(:,2)+CC(2,3)*tmp2(:,3)+CC(2,4)*tmp2(:,4)
             funf(i,j,k)= CC(1,1)*tmp1(1)+CC(1,2)*tmp1(2)+CC(1,3)*tmp1(3)+CC(1,4)*tmp1(4)
    enddo
   enddo
  enddo

  do k = kmino+1,kmaxo,2
   do j = jmino,jmaxo,2
    do i = imino,imaxo,2

       cxI(1) = i
       cxI(2) = j
       cxI(3) = k
! change to coarse level reference
!|---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*---| 
!|=======x===============x===============x===============x=======|
       cxI = (cxI+lbf-1)/2
! change to array index      
       cxI = cxI - lbc + 1
!  due to ghost zone, we can deal with symmetry boundary like this                   
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= CC(3,4)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)-1)+&
                      CC(3,3)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)  )+&
                      CC(3,2)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+1)+&
                      CC(3,1)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+2)
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= CC(3,4)*ya(:,:,1)+CC(3,3)*ya(:,:,2)+CC(3,2)*ya(:,:,3)+CC(3,1)*ya(:,:,4)
             endif
             tmp1= CC(2,1)*tmp2(:,1)+CC(2,2)*tmp2(:,2)+CC(2,3)*tmp2(:,3)+CC(2,4)*tmp2(:,4)
             funf(i,j,k)= CC(1,1)*tmp1(1)+CC(1,2)*tmp1(2)+CC(1,3)*tmp1(3)+CC(1,4)*tmp1(4)
    enddo
   enddo
  enddo

  do k = kmino,kmaxo,2
   do j = jmino+1,jmaxo,2
    do i = imino,imaxo,2

       cxI(1) = i
       cxI(2) = j
       cxI(3) = k
! change to coarse level reference
!|---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*---| 
!|=======x===============x===============x===============x=======|
       cxI = (cxI+lbf-1)/2
! change to array index      
       cxI = cxI - lbc + 1
!  due to ghost zone, we can deal with symmetry boundary like this                   
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= CC(3,1)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)-1)+&
                      CC(3,2)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)  )+&
                      CC(3,3)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+1)+&
                      CC(3,4)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+2)
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= CC(3,1)*ya(:,:,1)+CC(3,2)*ya(:,:,2)+CC(3,3)*ya(:,:,3)+CC(3,4)*ya(:,:,4)
             endif
             tmp1= CC(2,4)*tmp2(:,1)+CC(2,3)*tmp2(:,2)+CC(2,2)*tmp2(:,3)+CC(2,1)*tmp2(:,4)
             funf(i,j,k)= CC(1,1)*tmp1(1)+CC(1,2)*tmp1(2)+CC(1,3)*tmp1(3)+CC(1,4)*tmp1(4)
    enddo
   enddo
  enddo

  do k = kmino+1,kmaxo,2
   do j = jmino+1,jmaxo,2
    do i = imino,imaxo,2

       cxI(1) = i
       cxI(2) = j
       cxI(3) = k
! change to coarse level reference
!|---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*---| 
!|=======x===============x===============x===============x=======|
       cxI = (cxI+lbf-1)/2
! change to array index      
       cxI = cxI - lbc + 1
!  due to ghost zone, we can deal with symmetry boundary like this                   
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= CC(3,4)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)-1)+&
                      CC(3,3)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)  )+&
                      CC(3,2)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+1)+&
                      CC(3,1)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+2)
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= CC(3,4)*ya(:,:,1)+CC(3,3)*ya(:,:,2)+CC(3,2)*ya(:,:,3)+CC(3,1)*ya(:,:,4)
             endif
             tmp1= CC(2,4)*tmp2(:,1)+CC(2,3)*tmp2(:,2)+CC(2,2)*tmp2(:,3)+CC(2,1)*tmp2(:,4)
             funf(i,j,k)= CC(1,1)*tmp1(1)+CC(1,2)*tmp1(2)+CC(1,3)*tmp1(3)+CC(1,4)*tmp1(4)
    enddo
   enddo
  enddo

  do k = kmino,kmaxo,2
   do j = jmino,jmaxo,2
    do i = imino+1,imaxo,2

       cxI(1) = i
       cxI(2) = j
       cxI(3) = k
! change to coarse level reference
!|---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*---| 
!|=======x===============x===============x===============x=======|
       cxI = (cxI+lbf-1)/2
! change to array index      
       cxI = cxI - lbc + 1
!  due to ghost zone, we can deal with symmetry boundary like this                   
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= CC(3,1)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)-1)+&
                      CC(3,2)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)  )+&
                      CC(3,3)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+1)+&
                      CC(3,4)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+2)
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= CC(3,1)*ya(:,:,1)+CC(3,2)*ya(:,:,2)+CC(3,3)*ya(:,:,3)+CC(3,4)*ya(:,:,4)
             endif
             tmp1= CC(2,1)*tmp2(:,1)+CC(2,2)*tmp2(:,2)+CC(2,3)*tmp2(:,3)+CC(2,4)*tmp2(:,4)
             funf(i,j,k)= CC(1,4)*tmp1(1)+CC(1,3)*tmp1(2)+CC(1,2)*tmp1(3)+CC(1,1)*tmp1(4)
    enddo
   enddo
  enddo

  do k = kmino+1,kmaxo,2
   do j = jmino,jmaxo,2
    do i = imino+1,imaxo,2

       cxI(1) = i
       cxI(2) = j
       cxI(3) = k
! change to coarse level reference
!|---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*---| 
!|=======x===============x===============x===============x=======|
       cxI = (cxI+lbf-1)/2
! change to array index      
       cxI = cxI - lbc + 1
!  due to ghost zone, we can deal with symmetry boundary like this                   
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= CC(3,4)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)-1)+&
                      CC(3,3)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)  )+&
                      CC(3,2)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+1)+&
                      CC(3,1)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+2)
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= CC(3,4)*ya(:,:,1)+CC(3,3)*ya(:,:,2)+CC(3,2)*ya(:,:,3)+CC(3,1)*ya(:,:,4)
             endif
             tmp1= CC(2,1)*tmp2(:,1)+CC(2,2)*tmp2(:,2)+CC(2,3)*tmp2(:,3)+CC(2,4)*tmp2(:,4)
             funf(i,j,k)= CC(1,4)*tmp1(1)+CC(1,3)*tmp1(2)+CC(1,2)*tmp1(3)+CC(1,1)*tmp1(4)
    enddo
   enddo
  enddo


  do k = kmino,kmaxo,2
   do j = jmino+1,jmaxo,2
    do i = imino+1,imaxo,2

       cxI(1) = i
       cxI(2) = j
       cxI(3) = k
! change to coarse level reference
!|---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*---| 
!|=======x===============x===============x===============x=======|
       cxI = (cxI+lbf-1)/2
! change to array index      
       cxI = cxI - lbc + 1
!  due to ghost zone, we can deal with symmetry boundary like this                   
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= CC(3,1)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)-1)+&
                      CC(3,2)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)  )+&
                      CC(3,3)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+1)+&
                      CC(3,4)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+2)
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= CC(3,1)*ya(:,:,1)+CC(3,2)*ya(:,:,2)+CC(3,3)*ya(:,:,3)+CC(3,4)*ya(:,:,4)
             endif
             tmp1= CC(2,4)*tmp2(:,1)+CC(2,3)*tmp2(:,2)+CC(2,2)*tmp2(:,3)+CC(2,1)*tmp2(:,4)
             funf(i,j,k)= CC(1,4)*tmp1(1)+CC(1,3)*tmp1(2)+CC(1,2)*tmp1(3)+CC(1,1)*tmp1(4)
    enddo
   enddo
  enddo

  do k = kmino+1,kmaxo,2
   do j = jmino+1,jmaxo,2
    do i = imino+1,imaxo,2

       cxI(1) = i
       cxI(2) = j
       cxI(3) = k
! change to coarse level reference
!|---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*---| 
!|=======x===============x===============x===============x=======|
       cxI = (cxI+lbf-1)/2
! change to array index      
       cxI = cxI - lbc + 1
!  due to ghost zone, we can deal with symmetry boundary like this                   
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= CC(3,4)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)-1)+&
                      CC(3,3)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)  )+&
                      CC(3,2)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+1)+&
                      CC(3,1)*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+2)
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= CC(3,4)*ya(:,:,1)+CC(3,3)*ya(:,:,2)+CC(3,2)*ya(:,:,3)+CC(3,1)*ya(:,:,4)
             endif
             tmp1= CC(2,4)*tmp2(:,1)+CC(2,3)*tmp2(:,2)+CC(2,2)*tmp2(:,3)+CC(2,1)*tmp2(:,4)
             funf(i,j,k)= CC(1,4)*tmp1(1)+CC(1,3)*tmp1(2)+CC(1,2)*tmp1(3)+CC(1,1)*tmp1(4)
    enddo
   enddo
  enddo

  return

  end subroutine prolong3new
#else
#error Not define Vertex nor Cell
#endif  
#endif

#elif (ghost_width == 3)
! fourth order code
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
!--------------------------------------------------------------------------
!
! Prolongation from coarser grids to finer grids
! 6 points, 5th order interpolation
! 1   2   3   4   5   6
! *---*---*---*---*---*
!           ^
! f=3/256*(f_1+f_6) - 25/256*(f_2+f_5) + 75/128*(f_3+f_4)
!--------------------------------------------------------------------------

  subroutine prolong3(wei,llbc,uubc,extc,func,&
                      llbf,uubf,extf,funf,&
                      llbp,uubp,SoA,Symmetry)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: wei
!                                       coarse    fine       coarse
  real*8,dimension(3),   intent(in) :: llbc,uubc,llbf,uubf,llbp,uubp
  integer,dimension(3),  intent(in) :: extc,extf
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in)   :: func
  real*8, dimension(extf(1),extf(2),extf(3)),intent(inout):: funf
  real*8, dimension(1:3), intent(in) :: SoA
  integer,intent(in)::Symmetry

!~~~~~~> local variables

  real*8, dimension(1:3) :: base
  integer,dimension(3) :: lbc,ubc,lbf,ubf,lbp,ubp,lbpc,ubpc
  integer,dimension(3) :: cxB,cxT,cxI
  integer :: i,j,k,ii,jj,kk
  real*8, dimension(6,6) :: tmp2
  real*8, dimension(6) :: tmp1

  real*8, parameter :: C1=3.d0/2.56d2,C2=-2.5d1/2.56d2,C3=7.5d1/1.28d2

  integer::imini,imaxi,jmini,jmaxi,kmini,kmaxi
  integer::imino,imaxo,jmino,jmaxo,kmino,kmaxo
  real*8, dimension(-1:extc(1),-1:extc(2),-1:extc(3))   :: funcc

  real*8,dimension(3) :: CD,FD

  if(wei.ne.3)then
     write(*,*)"prolongrestrict.f90::prolong3: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

! it's possible a iolated point for target but not for source
  CD = (uubc-llbc)/(extc-1)
  FD = CD/2

!take care the mismatch of the two segments of grid
  do i=1,3
     if(llbc(i) <= llbf(i))then
        base(i) = llbc(i)
     else
        j=idint((llbc(i)-llbf(i))/FD(i)+0.4)
        if(j/2*2 == j)then
           base(i) = llbf(i)
        else
           base(i) = llbf(i) - CD(i)/2
        endif
     endif
  enddo

!!! function idint:
!If A is of type REAL and |A| < 1, INT(A) equals 0. If |A| \geq 1, 
!then INT(A) equals the largest integer that does not exceed the range of A 
!and whose sign is the same as the sign of A.

    lbf = idint((llbf-base)/FD+0.4)+1
    ubf = idint((uubf-base)/FD+0.4)+1
    lbc = idint((llbc-base)/CD+0.4)+1
    ubc = idint((uubc-base)/CD+0.4)+1
    lbp = idint((llbp-base)/FD+0.4)+1
    lbpc = idint((llbp-base)/CD+0.4)+1
    ubp = idint((uubp-base)/FD+0.4)+1
    ubpc = idint((uubp-base)/CD+0.4)+1
!sanity check
  imino=lbp(1)-lbf(1) + 1
  imaxo=ubp(1)-lbf(1) + 1
  jmino=lbp(2)-lbf(2) + 1
  jmaxo=ubp(2)-lbf(2) + 1
  kmino=lbp(3)-lbf(3) + 1
  kmaxo=ubp(3)-lbf(3) + 1

  imini=lbpc(1)-lbc(1) + 1
  imaxi=ubpc(1)-lbc(1) + 1
  jmini=lbpc(2)-lbc(2) + 1
  jmaxi=ubpc(2)-lbc(2) + 1
  kmini=lbpc(3)-lbc(3) + 1
  kmaxi=ubpc(3)-lbc(3) + 1

  if(imino.lt.1.or.jmino.lt.1.or.kmino.lt.1.or.&
     imini.lt.1.or.jmini.lt.1.or.kmini.lt.1.or.&
     imaxo.gt.extf(1).or.jmaxo.gt.extf(2).or.kmaxo.gt.extf(3).or.&
     imaxi.gt.extc(1)-2.or.jmaxi.gt.extc(2)-2.or.kmaxi.gt.extc(3)-2)then
          write(*,*)"error in prolongation for"
          write(*,*)"from"
          write(*,*)llbc,uubc
          write(*,*)lbc,ubc
          write(*,*)"to"
          write(*,*)llbf,uubf
          write(*,*)lbf,ubf
          write(*,*)"want"
          write(*,*)llbp,uubp
          write(*,*)lbp,ubp,lbpc,ubpc
          if(imini.lt.1) write(*,*)"imini = ",imini
          if(jmini.lt.1) write(*,*)"jmini = ",jmini       
          if(kmini.lt.1) write(*,*)"kmini = ",kmini      
          if(imino.lt.1) write(*,*)"imino = ",imino       
          if(jmino.lt.1) write(*,*)"jmino = ",jmino       
          if(kmino.lt.1) write(*,*)"kmino = ",kmino   
          if(imaxi.gt.extc(1)) write(*,*)"imaxi = ",imaxi,"extc(1) = ",extc(1) 
          if(jmaxi.gt.extc(2)) write(*,*)"jmaxi = ",jmaxi,"extc(2) = ",extc(2) 
          if(kmaxi.gt.extc(3)) write(*,*)"kmaxi = ",kmaxi,"extc(3) = ",extc(3) 
          if(imaxo.gt.extf(1)) write(*,*)"imaxo = ",imaxo,"extf(1) = ",extf(1) 
          if(jmaxo.gt.extf(2)) write(*,*)"jmaxo = ",jmaxo,"extf(2) = ",extf(2) 
          if(kmaxo.gt.extf(3)) write(*,*)"kmaxo = ",kmaxo,"extf(3) = ",extf(3) 
          return
  endif
  
  call symmetry_bd(2,extc,func,funcc,SoA)

!~~~~~~> prolongation start...
  do k = kmino,kmaxo
   do j = jmino,jmaxo
    do i = imino,imaxo
! change to coarse level reference        v
!|*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*| 
!|x===============x===============x===============x===============x===============x|
       cxI(1) = i
       cxI(2) = j
       cxI(3) = k
       cxI = (cxI+lbf)/2
! change to array index      
       cxI = cxI - lbc + 1

       ii=i+lbf(1)-1
       jj=j+lbf(2)-1
       kk=k+lbf(3)-1

       if(any(cxI+3 > extc)) write(*,*)"error in prolong"
       if(ii/2*2==ii)then
         if(jj/2*2==jj)then
           if(kk/2*2==kk)then
!  due to ghost zone, we can deal with symmetry boundary like this       
             tmp2= C1*(funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-2)+funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+3))&
                  +C2*(funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-1)+funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+2))&
                  +C3*(funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)  )+funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+1))
             tmp1= C1*(tmp2(:,1)+tmp2(:,6))+C2*(tmp2(:,2)+tmp2(:,5))+C3*(tmp2(:,3)+tmp2(:,4))
             funf(i,j,k)= C1*(tmp1(1)+tmp1(6))+C2*(tmp1(2)+tmp1(5))+C3*(tmp1(3)+tmp1(4))
           else
             tmp2= funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3))            
             tmp1= C1*(tmp2(:,1)+tmp2(:,6))+C2*(tmp2(:,2)+tmp2(:,5))+C3*(tmp2(:,3)+tmp2(:,4))
             funf(i,j,k)= C1*(tmp1(1)+tmp1(6))+C2*(tmp1(2)+tmp1(5))+C3*(tmp1(3)+tmp1(4))
           endif
         else
           if(kk/2*2==kk)then
             tmp2= C1*(funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-2)+funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+3))&
                  +C2*(funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-1)+funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+2))&
                  +C3*(funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)  )+funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+1))
             tmp1= tmp2(:,3)
             funf(i,j,k)= C1*(tmp1(1)+tmp1(6))+C2*(tmp1(2)+tmp1(5))+C3*(tmp1(3)+tmp1(4))
           else
             tmp2= funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3))  
             tmp1= tmp2(:,3)
             funf(i,j,k)= C1*(tmp1(1)+tmp1(6))+C2*(tmp1(2)+tmp1(5))+C3*(tmp1(3)+tmp1(4))
           endif
         endif
       else
         if(jj/2*2==jj)then
           if(kk/2*2==kk)then    
             tmp2= C1*(funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-2)+funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+3))&
                  +C2*(funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-1)+funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+2))&
                  +C3*(funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)  )+funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+1))
             tmp1= C1*(tmp2(:,1)+tmp2(:,6))+C2*(tmp2(:,2)+tmp2(:,5))+C3*(tmp2(:,3)+tmp2(:,4))
             funf(i,j,k)= tmp1(3)
           else
             tmp2= funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3))  
             tmp1= C1*(tmp2(:,1)+tmp2(:,6))+C2*(tmp2(:,2)+tmp2(:,5))+C3*(tmp2(:,3)+tmp2(:,4))
             funf(i,j,k)= tmp1(3)
           endif
         else
           if(kk/2*2==kk)then
             tmp2= C1*(funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-2)+funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+3))&
                  +C2*(funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-1)+funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+2))&
                  +C3*(funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)  )+funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+1))
             tmp1= tmp2(:,3)
             funf(i,j,k)= tmp1(3)
           else
             funf(i,j,k)= funcc(cxI(1),cxI(2),cxI(3))
           endif
         endif
       endif
    enddo
   enddo
  enddo

  return

  end subroutine prolong3
#else
#ifdef Cell
!--------------------------------------------------------------------------
!
! Prolongation from coarser grids to finer grids
! 6 points, 5th order interpolation
! 1   2   3   4   5   6
! *---*---*---*---*---*
!          ^
! f=77/8192*f_1 - 693/8192*f_2 + 3465/4096*f_3 +
!   63/8192*f_6 - 495/8192*f_5 + 1155/4096*f_4
!--------------------------------------------------------------------------

  subroutine prolong3(wei,llbc,uubc,extc,func,&
                      llbf,uubf,extf,funf,&
                      llbp,uubp,SoA,Symmetry)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: wei
!                                       coarse    fine       coarse
  real*8,dimension(3),   intent(in) :: llbc,uubc,llbf,uubf,llbp,uubp
  integer,dimension(3),  intent(in) :: extc,extf
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in)   :: func
  real*8, dimension(extf(1),extf(2),extf(3)),intent(inout):: funf
  real*8, dimension(1:3), intent(in) :: SoA
  integer,intent(in)::Symmetry

!~~~~~~> local variables

  real*8, dimension(1:3) :: base
  integer,dimension(3) :: lbc,ubc,lbf,ubf,lbp,ubp,lbpc,ubpc
! when if=1 -> ic=0, this is different to vertex center grid 
  real*8, dimension(-2:extc(1),-2:extc(2),-2:extc(3))   :: funcc
  integer,dimension(3) :: cxI
  integer :: i,j,k,ii,jj,kk
  real*8, dimension(6,6) :: tmp2
  real*8, dimension(6) :: tmp1

  real*8, parameter :: C1=7.7d1/8.192d3,C2=-6.93d2/8.192d3,C3=3.465d3/4.096d3
  real*8, parameter :: C6=6.3d1/8.192d3,C5=-4.95d2/8.192d3,C4=1.155d3/4.096d3

  integer::imini,imaxi,jmini,jmaxi,kmini,kmaxi
  integer::imino,imaxo,jmino,jmaxo,kmino,kmaxo

  real*8,dimension(3) :: CD,FD
  
  if(wei.ne.3)then
     write(*,*)"prolongrestrict.f90::prolong3: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

  CD = (uubc-llbc)/extc
  FD = (uubf-llbf)/extf

!take care the mismatch of the two segments of grid
  do i=1,3
     if(llbc(i) <= llbf(i))then
        base(i) = llbc(i)
     else
        j=idint((llbc(i)-llbf(i))/FD(i)+0.4)
        if(j/2*2 == j)then
           base(i) = llbf(i)
        else
           base(i) = llbf(i) - CD(i)/2
        endif
     endif
  enddo

!!! function idint:
!If A is of type REAL and |A| < 1, INT(A) equals 0. If |A| \geq 1, 
!then INT(A) equals the largest integer that does not exceed the range of A 
!and whose sign is the same as the sign of A.

    lbf = idint((llbf-base)/FD+0.4)+1
    ubf = idint((uubf-base)/FD+0.4)
    lbc = idint((llbc-base)/CD+0.4)+1
    ubc = idint((uubc-base)/CD+0.4)
    lbp = idint((llbp-base)/FD+0.4)+1
    lbpc = idint((llbp-base)/CD+0.4)+1
    ubp = idint((uubp-base)/FD+0.4)
    ubpc = idint((uubp-base)/CD+0.4)

!sanity check
  imino=lbp(1)-lbf(1) + 1
  imaxo=ubp(1)-lbf(1) + 1
  jmino=lbp(2)-lbf(2) + 1
  jmaxo=ubp(2)-lbf(2) + 1
  kmino=lbp(3)-lbf(3) + 1
  kmaxo=ubp(3)-lbf(3) + 1

  imini=lbpc(1)-lbc(1) + 1
  imaxi=ubpc(1)-lbc(1) + 1
  jmini=lbpc(2)-lbc(2) + 1
  jmaxi=ubpc(2)-lbc(2) + 1
  kmini=lbpc(3)-lbc(3) + 1
  kmaxi=ubpc(3)-lbc(3) + 1

  if(imino.lt.1.or.jmino.lt.1.or.kmino.lt.1.or.&
     imini.lt.1.or.jmini.lt.1.or.kmini.lt.1.or.&
     imaxo.gt.extf(1).or.jmaxo.gt.extf(2).or.kmaxo.gt.extf(3).or.&
     imaxi.gt.extc(1)-2.or.jmaxi.gt.extc(2)-2.or.kmaxi.gt.extc(3)-2)then
          write(*,*)"error in prolongation for"
          write(*,*)"from"
          write(*,*)llbc,uubc
          write(*,*)lbc,ubc
          write(*,*)"to"
          write(*,*)llbf,uubf
          write(*,*)lbf,ubf
          write(*,*)"want"
          write(*,*)llbp,uubp
          write(*,*)lbp,ubp,lbpc,ubpc
          return
  endif

  call symmetry_bd(3,extc,func,funcc,SoA)
     
!~~~~~~> prolongation start...
  do k = kmino,kmaxo
   do j = jmino,jmaxo
    do i = imino,imaxo
       cxI(1) = i
       cxI(2) = j
       cxI(3) = k
! change to coarse level reference
!|---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*---| 
!|=======x===============x===============x===============x=======|
       cxI = (cxI+lbf-1)/2
! change to array index      
       cxI = cxI - lbc + 1

       if(any(cxI+3 > extc)) write(*,*)"error in prolong"
       ii=i+lbf(1)-1
       jj=j+lbf(2)-1
       kk=k+lbf(3)-1
       if(ii/2*2==ii)then
         if(jj/2*2==jj)then
           if(kk/2*2==kk)then
             tmp2= C1*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-2)+&
                   C2*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-1)+&
                   C3*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)  )+&
                   C4*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+1)+&
                   C5*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+2)+&
                   C6*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+3)
             tmp1= C1*tmp2(:,1)+C2*tmp2(:,2)+C3*tmp2(:,3)+C4*tmp2(:,4)+C5*tmp2(:,5)+C6*tmp2(:,6)
             funf(i,j,k)= C1*tmp1(1)+C2*tmp1(2)+C3*tmp1(3)+C4*tmp1(4)+C5*tmp1(5)+C6*tmp1(6)
           else
             tmp2= C6*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-2)+&
                   C5*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-1)+&
                   C4*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)  )+&
                   C3*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+1)+&
                   C2*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+2)+&
                   C1*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+3)
             tmp1= C1*tmp2(:,1)+C2*tmp2(:,2)+C3*tmp2(:,3)+C4*tmp2(:,4)+C5*tmp2(:,5)+C6*tmp2(:,6)
             funf(i,j,k)=  C1*tmp1(1)+C2*tmp1(2)+C3*tmp1(3)+C4*tmp1(4)+C5*tmp1(5)+C6*tmp1(6)
           endif
         else
           if(kk/2*2==kk)then
             tmp2= C1*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-2)+&
                   C2*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-1)+&
                   C3*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)  )+&
                   C4*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+1)+&
                   C5*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+2)+&
                   C6*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+3)
             tmp1= C6*tmp2(:,1)+C5*tmp2(:,2)+C4*tmp2(:,3)+C3*tmp2(:,4)+C2*tmp2(:,5)+C1*tmp2(:,6)
             funf(i,j,k)= C1*tmp1(1)+C2*tmp1(2)+C3*tmp1(3)+C4*tmp1(4)+C5*tmp1(5)+C6*tmp1(6)
           else
             tmp2= C6*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-2)+&
                   C5*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-1)+&
                   C4*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)  )+&
                   C3*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+1)+&
                   C2*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+2)+&
                   C1*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+3)
             tmp1= C6*tmp2(:,1)+C5*tmp2(:,2)+C4*tmp2(:,3)+C3*tmp2(:,4)+C2*tmp2(:,5)+C1*tmp2(:,6)
             funf(i,j,k)=  C1*tmp1(1)+C2*tmp1(2)+C3*tmp1(3)+C4*tmp1(4)+C5*tmp1(5)+C6*tmp1(6)
           endif
         endif
       else
         if(jj/2*2==jj)then
           if(kk/2*2==kk)then               
             tmp2= C1*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-2)+&
                   C2*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-1)+&
                   C3*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)  )+&
                   C4*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+1)+&
                   C5*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+2)+&
                   C6*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+3)
             tmp1= C1*tmp2(:,1)+C2*tmp2(:,2)+C3*tmp2(:,3)+C4*tmp2(:,4)+C5*tmp2(:,5)+C6*tmp2(:,6)
             funf(i,j,k)= C6*tmp1(1)+C5*tmp1(2)+C4*tmp1(3)+C3*tmp1(4)+C2*tmp1(5)+C1*tmp1(6)
           else
             tmp2= C6*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-2)+&
                   C5*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-1)+&
                   C4*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)  )+&
                   C3*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+1)+&
                   C2*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+2)+&
                   C1*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+3)
             tmp1= C1*tmp2(:,1)+C2*tmp2(:,2)+C3*tmp2(:,3)+C4*tmp2(:,4)+C5*tmp2(:,5)+C6*tmp2(:,6)
             funf(i,j,k)=  C6*tmp1(1)+C5*tmp1(2)+C4*tmp1(3)+C3*tmp1(4)+C2*tmp1(5)+C1*tmp1(6)
           endif
         else
           if(kk/2*2==kk)then
             tmp2= C1*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-2)+&
                   C2*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-1)+&
                   C3*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)  )+&
                   C4*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+1)+&
                   C5*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+2)+&
                   C6*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+3)
             tmp1= C6*tmp2(:,1)+C5*tmp2(:,2)+C4*tmp2(:,3)+C3*tmp2(:,4)+C2*tmp2(:,5)+C1*tmp2(:,6)
             funf(i,j,k)= C6*tmp1(1)+C5*tmp1(2)+C4*tmp1(3)+C3*tmp1(4)+C2*tmp1(5)+C1*tmp1(6)
           else
             tmp2= C6*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-2)+&
                   C5*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-1)+&
                   C4*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)  )+&
                   C3*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+1)+&
                   C2*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+2)+&
                   C1*funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+3)
             tmp1= C6*tmp2(:,1)+C5*tmp2(:,2)+C4*tmp2(:,3)+C3*tmp2(:,4)+C2*tmp2(:,5)+C1*tmp2(:,6)
             funf(i,j,k)=  C6*tmp1(1)+C5*tmp1(2)+C4*tmp1(3)+C3*tmp1(4)+C2*tmp1(5)+C1*tmp1(6)
           endif
         endif
       endif
    enddo
   enddo
  enddo

  return

  end subroutine prolong3
!--------------------------------------------------------------------------
!
! Restrict from finner grids to coarser grids ignore the boundary point
!
! 6 points, 5th order interpolation
! 1   2   3   4   5   6
! *---*---*---*---*---*
!           ^
! f=3/256*(f_1+f_6) - 25/256*(f_2+f_5) + 75/128*(f_3+f_4)
!--------------------------------------------------------------------------

  subroutine restrict3(wei,llbc,uubc,extc,func,&
                       llbf,uubf,extf,funf,&
                       llbr,uubr,SoA,Symmetry)
  implicit none

!~~~~~~> input arguments
  integer,intent(in)::wei
!                                       coarse    fine       coarse
  real*8,dimension(3),   intent(in) :: llbc,uubc,llbf,uubf,llbr,uubr
  integer,dimension(3),  intent(in) :: extc,extf
  real*8, dimension(extc(1),extc(2),extc(3)),intent(inout):: func
  real*8, dimension(extf(1),extf(2),extf(3)),intent(in):: funf
  real*8, dimension(1:3), intent(in) :: SoA
  integer,intent(in)::Symmetry

!~~~~~~> local variables

  real*8, dimension(1:3) :: base
  integer,dimension(3) :: lbc,ubc,lbf,ubf,lbr,ubr,lbrf,ubrf
  real*8, dimension(-1:extf(1),-1:extf(2),-1:extf(3)):: funff
  integer,dimension(3) :: cxI
  integer :: i,j,k
  real*8, dimension(6,6) :: tmp2
  real*8, dimension(6) :: tmp1
  real*8, parameter :: C1=3.d0/2.56d2,C2=-2.5d1/2.56d2,C3=7.5d1/1.28d2

  integer::imini,imaxi,jmini,jmaxi,kmini,kmaxi
  integer::imino,imaxo,jmino,jmaxo,kmino,kmaxo

  real*8,dimension(3) :: CD,FD
  
  if(wei.ne.3)then
     write(*,*)"prolongrestrict.f90::restrict3: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

  CD = (uubc-llbc)/extc
  FD = (uubf-llbf)/extf

!take care the mismatch of the two segments of grid
  do i=1,3
     if(llbc(i) <= llbf(i))then
        base(i) = llbc(i)
     else
        j=idint((llbc(i)-llbf(i))/FD(i)+0.4)
        if(j/2*2 == j)then
           base(i) = llbf(i)
        else
           base(i) = llbf(i) - CD(i)/2
        endif
     endif
  enddo
!!! function idint:
!If A is of type REAL and |A| < 1, INT(A) equals 0. If |A| \geq 1, 
!then INT(A) equals the largest integer that does not exceed the range of A 
!and whose sign is the same as the sign of A.

! note say base = 0, llbf = 0, uubf = 2
! llbf->1 and uubf->2
    lbf = idint((llbf-base)/FD+0.4)+1
    ubf = idint((uubf-base)/FD+0.4)
    lbc = idint((llbc-base)/CD+0.4)+1
    ubc = idint((uubc-base)/CD+0.4)
    lbr = idint((llbr-base)/CD+0.4)+1
    lbrf = idint((llbr-base)/FD+0.4)+1
    ubr = idint((uubr-base)/CD+0.4)
    ubrf = idint((uubr-base)/FD+0.4)

!sanity check
  imino=lbr(1)-lbc(1) + 1
  imaxo=ubr(1)-lbc(1) + 1
  jmino=lbr(2)-lbc(2) + 1
  jmaxo=ubr(2)-lbc(2) + 1
  kmino=lbr(3)-lbc(3) + 1
  kmaxo=ubr(3)-lbc(3) + 1

  imini=lbrf(1)-lbf(1) + 1
  imaxi=ubrf(1)-lbf(1) + 1
  jmini=lbrf(2)-lbf(2) + 1
  jmaxi=ubrf(2)-lbf(2) + 1
  kmini=lbrf(3)-lbf(3) + 1
  kmaxi=ubrf(3)-lbf(3) + 1

  if(imino.lt.1.or.jmino.lt.1.or.kmino.lt.1.or.&
     imini.lt.1.or.jmini.lt.1.or.kmini.lt.1.or.&
     imaxo.gt.extc(1).or.jmaxo.gt.extc(2).or.kmaxo.gt.extc(3).or.&
     imaxi.gt.extf(1)-2.or.jmaxi.gt.extf(2)-2.or.kmaxi.gt.extf(3)-2)then
          write(*,*)"error in restrict for"
          write(*,*)"from"
          write(*,*)lbf,ubf
          write(*,*)"to"
          write(*,*)lbc,ubc
          write(*,*)"want"
          write(*,*)lbr,ubr,lbrf,ubrf
          write(*,*)"llbf = ",llbf
          write(*,*)"uubf = ",uubf
          write(*,*)"llbc = ",llbc
          write(*,*)"uubc = ",uubc
          write(*,*)"llbr = ",llbr
          write(*,*)"uubr = ",uubr
          write(*,*)"base = ",base
          stop
  endif

  call symmetry_bd(2,extf,funf,funff,SoA)

!~~~~~~> restriction start...
  do k = kmino,kmaxo
   do j = jmino,jmaxo
    do i = imino,imaxo

       cxI(1) = i
       cxI(2) = j
       cxI(3) = k
! change to fine level reference
!|---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*---| 
!|=======x===============x===============x===============x=======|
       cxI = 2*(cxI+lbc-1) - 1
! change to array index      
       cxI = cxI - lbf + 1

       if(any(cxI+3 > extf)) write(*,*)"error in restrict"
       tmp2= C1*(funff(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-2)+funff(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+3))&
            +C2*(funff(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-1)+funff(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+2))&
            +C3*(funff(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)  )+funff(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+1))
       tmp1= C1*(tmp2(:,1)+tmp2(:,6))+C2*(tmp2(:,2)+tmp2(:,5))+C3*(tmp2(:,3)+tmp2(:,4))
       func(i,j,k)= C1*(tmp1(1)+tmp1(6))+C2*(tmp1(2)+tmp1(5))+C3*(tmp1(3)+tmp1(4))
    enddo
   enddo
  enddo
  
  return

  end subroutine restrict3
#else
#error Not define Vertex nor Cell
#endif  
#endif

#elif (ghost_width == 4)
! sixth order code
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
!--------------------------------------------------------------------------
!
! Prolongation from coarser grids to finer grids
! 8 points, 7th order interpolation
! 1   2   3   4   5   6   7   8
! *---*---*---*---*---*---*---*
!               ^
! f=-5/2048*(f_1+f_8) + 49/2048*(f_2+f_7) - 245/2048*(f_3+f_6) + 1225/2048*(f_4+f_5)
!--------------------------------------------------------------------------

  subroutine prolong3(wei,llbc,uubc,extc,func,&
                      llbf,uubf,extf,funf,&
                      llbp,uubp,SoA,Symmetry)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: wei
!                                       coarse    fine       coarse
  real*8,dimension(3),   intent(in) :: llbc,uubc,llbf,uubf,llbp,uubp
  integer,dimension(3),  intent(in) :: extc,extf
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in)   :: func
  real*8, dimension(extf(1),extf(2),extf(3)),intent(inout):: funf
  real*8, dimension(1:3), intent(in) :: SoA
  integer,intent(in)::Symmetry

!~~~~~~> local variables

  real*8, dimension(1:3) :: base
  integer,dimension(3) :: lbc,ubc,lbf,ubf,lbp,ubp,lbpc,ubpc
  integer,dimension(3) :: cxB,cxT,cxI
  integer :: i,j,k,ii,jj,kk
  real*8,  dimension(8,8,8) :: ya
  real*8, dimension(8,8) :: tmp2
  real*8, dimension(8) :: tmp1

  real*8, parameter :: C1=-5.d0/2.048d3,C2=4.9d1/2.048d3,C3=-2.45d2/2.048d3,C4=1.225d3/2.048d3

  integer::imini,imaxi,jmini,jmaxi,kmini,kmaxi
  integer::imino,imaxo,jmino,jmaxo,kmino,kmaxo
  logical::decide3d

  real*8,dimension(3) :: CD,FD

  if(wei.ne.3)then
     write(*,*)"prolongrestrict.f90::prolong3: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

! it's possible a iolated point for target but not for source
  CD = (uubc-llbc)/(extc-1)
  FD = CD/2

!take care the mismatch of the two segments of grid
  do i=1,3
     if(llbc(i) <= llbf(i))then
        base(i) = llbc(i)
     else
        j=idint((llbc(i)-llbf(i))/FD(i)+0.4)
        if(j/2*2 == j)then
           base(i) = llbf(i)
        else
           base(i) = llbf(i) - CD(i)/2
        endif
     endif
  enddo

!!! function idint:
!If A is of type REAL and |A| < 1, INT(A) equals 0. If |A| \geq 1, 
!then INT(A) equals the largest integer that does not exceed the range of A 
!and whose sign is the same as the sign of A.

    lbf = idint((llbf-base)/FD+0.4)+1
    ubf = idint((uubf-base)/FD+0.4)+1
    lbc = idint((llbc-base)/CD+0.4)+1
    ubc = idint((uubc-base)/CD+0.4)+1
    lbp = idint((llbp-base)/FD+0.4)+1
    lbpc = idint((llbp-base)/CD+0.4)+1
    ubp = idint((uubp-base)/FD+0.4)+1
    ubpc = idint((uubp-base)/CD+0.4)+1
!sanity check
  imino=lbp(1)-lbf(1) + 1
  imaxo=ubp(1)-lbf(1) + 1
  jmino=lbp(2)-lbf(2) + 1
  jmaxo=ubp(2)-lbf(2) + 1
  kmino=lbp(3)-lbf(3) + 1
  kmaxo=ubp(3)-lbf(3) + 1

  imini=lbpc(1)-lbc(1) + 1
  imaxi=ubpc(1)-lbc(1) + 1
  jmini=lbpc(2)-lbc(2) + 1
  jmaxi=ubpc(2)-lbc(2) + 1
  kmini=lbpc(3)-lbc(3) + 1
  kmaxi=ubpc(3)-lbc(3) + 1

  if(imino.lt.1.or.jmino.lt.1.or.kmino.lt.1.or.&
     imini.lt.1.or.jmini.lt.1.or.kmini.lt.1.or.&
     imaxo.gt.extf(1).or.jmaxo.gt.extf(2).or.kmaxo.gt.extf(3).or.&
     imaxi.gt.extc(1)-3.or.jmaxi.gt.extc(2)-3.or.kmaxi.gt.extc(3)-3)then
          write(*,*)"error in prolongation for"
          write(*,*)"from"
          write(*,*)llbc,uubc
          write(*,*)lbc,ubc
          write(*,*)"to"
          write(*,*)llbf,uubf
          write(*,*)lbf,ubf
          write(*,*)"want"
          write(*,*)llbp,uubp
          write(*,*)lbp,ubp,lbpc,ubpc
          if(imini.lt.1) write(*,*)"imini = ",imini
          if(jmini.lt.1) write(*,*)"jmini = ",jmini       
          if(kmini.lt.1) write(*,*)"kmini = ",kmini      
          if(imino.lt.1) write(*,*)"imino = ",imino       
          if(jmino.lt.1) write(*,*)"jmino = ",jmino       
          if(kmino.lt.1) write(*,*)"kmino = ",kmino   
          if(imaxi.gt.extc(1)) write(*,*)"imaxi = ",imaxi,"extc(1) = ",extc(1) 
          if(jmaxi.gt.extc(2)) write(*,*)"jmaxi = ",jmaxi,"extc(2) = ",extc(2) 
          if(kmaxi.gt.extc(3)) write(*,*)"kmaxi = ",kmaxi,"extc(3) = ",extc(3) 
          if(imaxo.gt.extf(1)) write(*,*)"imaxo = ",imaxo,"extf(1) = ",extf(1) 
          if(jmaxo.gt.extf(2)) write(*,*)"jmaxo = ",jmaxo,"extf(2) = ",extf(2) 
          if(kmaxo.gt.extf(3)) write(*,*)"kmaxo = ",kmaxo,"extf(3) = ",extf(3) 
          return
  endif
!~~~~~~> prolongation start...
  do k = kmino,kmaxo
   do j = jmino,jmaxo
    do i = imino,imaxo
! change to coarse level reference                        v
!|*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*| 
!|x===============x===============x===============x===============x===============x===============x===============x|
       cxI(1) = i
       cxI(2) = j
       cxI(3) = k
       cxI = (cxI+lbf)/2
! change to array index      
       cxI = cxI - lbc + 1

       ii=i+lbf(1)-1
       jj=j+lbf(2)-1
       kk=k+lbf(3)-1

       if(any(cxI+4 > extc)) write(*,*)"error in prolong"
       if(ii/2*2==ii)then
         if(jj/2*2==jj)then
           if(kk/2*2==kk)then
!  due to ghost zone, we can deal with symmetry boundary like this                  
             if(cxI(1)>3.and.cxI(2)>3.and.cxI(3)>3)then
               tmp2= C1*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-3)+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+4))&
                    +C2*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-2)+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+3))&
                    +C3*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-1)+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+2))&
                    +C4*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)  )+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+1))
             else
               cxB=cxI-3
               cxT=cxI+4
               if(decide3d(extc,func,func,cxB,cxT,SoA,ya,8,Symmetry))then
                 write(*,*)"prolong3 position: "
                 write(*,*)llbf(1)+(i-1)*FD(1),llbf(2)+(j-1)*FD(2),llbf(3)+(k-1)*FD(3)
                 write(*,*)"llbf = ",llbf
                 stop
               endif
               tmp2=C1*(ya(:,:,1)+ya(:,:,8))+C2*(ya(:,:,2)+ya(:,:,7))+C3*(ya(:,:,3)+ya(:,:,6))+C4*(ya(:,:,4)+ya(:,:,5))
             endif
             tmp1= C1*(tmp2(:,1)+tmp2(:,8))+C2*(tmp2(:,2)+tmp2(:,7))+C3*(tmp2(:,3)+tmp2(:,6))+C4*(tmp2(:,4)+tmp2(:,5))
             funf(i,j,k)= C1*(tmp1(1)+tmp1(8))+C2*(tmp1(2)+tmp1(7))+C3*(tmp1(3)+tmp1(6))+C4*(tmp1(4)+tmp1(5))
           else
             if(cxI(1)>3.and.cxI(2)>3.and.cxI(3)>3)then
                tmp2= func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3))
             else
                cxB=cxI-3
                cxT=cxI+4
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,8,Symmetry))then
                 write(*,*)"prolong3 position: "
                 write(*,*)llbf(1)+(i-1)*FD(1),llbf(2)+(j-1)*FD(2),llbf(3)+(k-1)*FD(3)
                 write(*,*)"llbf = ",llbf
                 stop
                endif
                tmp2= ya(:,:,4)
             endif
             tmp1= C1*(tmp2(:,1)+tmp2(:,8))+C2*(tmp2(:,2)+tmp2(:,7))+C3*(tmp2(:,3)+tmp2(:,6))+C4*(tmp2(:,4)+tmp2(:,5))
             funf(i,j,k)= C1*(tmp1(1)+tmp1(8))+C2*(tmp1(2)+tmp1(7))+C3*(tmp1(3)+tmp1(6))+C4*(tmp1(4)+tmp1(5))
           endif
         else
           if(kk/2*2==kk)then
             if(cxI(1)>3.and.cxI(2)>3.and.cxI(3)>3)then
               tmp2= C1*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-3)+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+4))&
                    +C2*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-2)+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+3))&
                    +C3*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-1)+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+2))&
                    +C4*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)  )+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+1))
             else
                cxB=cxI-3
                cxT=cxI+4
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,8,Symmetry))then
                 write(*,*)"prolong3 position: "
                 write(*,*)llbf(1)+(i-1)*FD(1),llbf(2)+(j-1)*FD(2),llbf(3)+(k-1)*FD(3)
                 write(*,*)"llbf = ",llbf
                 stop
                endif
                tmp2=C1*(ya(:,:,1)+ya(:,:,8))+C2*(ya(:,:,2)+ya(:,:,7))+C3*(ya(:,:,3)+ya(:,:,6))+C4*(ya(:,:,4)+ya(:,:,5))
             endif
             tmp1= tmp2(:,4)
             funf(i,j,k)= C1*(tmp1(1)+tmp1(8))+C2*(tmp1(2)+tmp1(7))+C3*(tmp1(3)+tmp1(6))+C4*(tmp1(4)+tmp1(5))
           else
             if(cxI(1)>3.and.cxI(2)>3.and.cxI(3)>3)then
                tmp2= func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3))
             else
                cxB=cxI-3
                cxT=cxI+4
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,8,Symmetry))then
                 write(*,*)"prolong3 position: "
                 write(*,*)llbf(1)+(i-1)*FD(1),llbf(2)+(j-1)*FD(2),llbf(3)+(k-1)*FD(3)
                 write(*,*)"llbf = ",llbf
                 stop
                endif
                tmp2= ya(:,:,4)
             endif
             tmp1= tmp2(:,4)
             funf(i,j,k)= C1*(tmp1(1)+tmp1(8))+C2*(tmp1(2)+tmp1(7))+C3*(tmp1(3)+tmp1(6))+C4*(tmp1(4)+tmp1(5))
           endif
         endif
       else
         if(jj/2*2==jj)then
           if(kk/2*2==kk)then               
             if(cxI(1)>3.and.cxI(2)>3.and.cxI(3)>3)then
               tmp2= C1*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-3)+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+4))&
                    +C2*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-2)+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+3))&
                    +C3*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-1)+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+2))&
                    +C4*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)  )+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+1))
             else
                cxB=cxI-3
                cxT=cxI+4
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,8,Symmetry))then
                 write(*,*)"prolong3 position: "
                 write(*,*)llbf(1)+(i-1)*FD(1),llbf(2)+(j-1)*FD(2),llbf(3)+(k-1)*FD(3)
                 write(*,*)"llbf = ",llbf
                 stop
                endif
                tmp2=C1*(ya(:,:,1)+ya(:,:,8))+C2*(ya(:,:,2)+ya(:,:,7))+C3*(ya(:,:,3)+ya(:,:,6))+C4*(ya(:,:,4)+ya(:,:,5))
             endif
             tmp1= C1*(tmp2(:,1)+tmp2(:,8))+C2*(tmp2(:,2)+tmp2(:,7))+C3*(tmp2(:,3)+tmp2(:,6))+C4*(tmp2(:,4)+tmp2(:,5))
             funf(i,j,k)= tmp1(4)
           else
             if(cxI(1)>3.and.cxI(2)>3.and.cxI(3)>3)then
                tmp2= func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3))
             else
                cxB=cxI-3
                cxT=cxI+4
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,8,Symmetry))then
                 write(*,*)"prolong3 position: "
                 write(*,*)llbf(1)+(i-1)*FD(1),llbf(2)+(j-1)*FD(2),llbf(3)+(k-1)*FD(3)
                 write(*,*)"llbf = ",llbf
                 stop
                endif
                tmp2= ya(:,:,4)
             endif
             tmp1= C1*(tmp2(:,1)+tmp2(:,8))+C2*(tmp2(:,2)+tmp2(:,7))+C3*(tmp2(:,3)+tmp2(:,6))+C4*(tmp2(:,4)+tmp2(:,5))
             funf(i,j,k)= tmp1(4)
           endif
         else
           if(kk/2*2==kk)then
             if(cxI(1)>3.and.cxI(2)>3.and.cxI(3)>3)then
               tmp2= C1*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-3)+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+4))&
                    +C2*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-2)+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+3))&
                    +C3*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-1)+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+2))&
                    +C4*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)  )+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+1))
             else
                cxB=cxI-3
                cxT=cxI+4
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,8,Symmetry))then
                 write(*,*)"prolong3 position: "
                 write(*,*)llbf(1)+(i-1)*FD(1),llbf(2)+(j-1)*FD(2),llbf(3)+(k-1)*FD(3)
                 write(*,*)"llbf = ",llbf
                 stop
                endif
                tmp2=C1*(ya(:,:,1)+ya(:,:,8))+C2*(ya(:,:,2)+ya(:,:,7))+C3*(ya(:,:,3)+ya(:,:,6))+C4*(ya(:,:,4)+ya(:,:,5))
             endif
             tmp1= tmp2(:,4)
             funf(i,j,k)= tmp1(4)
           else
             funf(i,j,k)= func(cxI(1),cxI(2),cxI(3))
           endif
         endif
       endif
    enddo
   enddo
  enddo

  return

  end subroutine prolong3

#else
#ifdef Cell
!--------------------------------------------------------------------------------------
!
! Restrict from finner grids to coarser grids ignore the boundary point
!
! 8 points, 7th order interpolation
! 1   2   3   4   5   6   7   8
! *---*---*---*---*---*---*---*
!               ^
! f=-5/2048*(f_1+f_8) + 49/2048*(f_2+f_7) - 245/2048*(f_3+f_6) + 1225/2048*(f_4+f_5)
!--------------------------------------------------------------------------------------

  subroutine restrict3(wei,llbc,uubc,extc,func,&
                       llbf,uubf,extf,funf,&
                       llbr,uubr,SoA,Symmetry)
  implicit none

!~~~~~~> input arguments
  integer,intent(in)::wei
!                                       coarse    fine       coarse
  real*8,dimension(3),   intent(in) :: llbc,uubc,llbf,uubf,llbr,uubr
  integer,dimension(3),  intent(in) :: extc,extf
  real*8, dimension(extc(1),extc(2),extc(3)),intent(inout):: func
  real*8, dimension(extf(1),extf(2),extf(3)),intent(in):: funf
  real*8, dimension(1:3), intent(in) :: SoA
  integer,intent(in)::Symmetry

!~~~~~~> local variables

  real*8, dimension(1:3) :: base
  integer,dimension(3) :: lbc,ubc,lbf,ubf,lbr,ubr,lbrf,ubrf
  integer,dimension(3) :: cxB,cxT,cxI
  integer :: i,j,k
  real*8,  dimension(8,8,8) :: ya
  real*8, dimension(8,8) :: tmp2
  real*8, dimension(8) :: tmp1
  real*8, parameter :: C1=-5.d0/2.048d3,C2=4.9d1/2.048d3,C3=-2.45d2/2.048d3,C4=1.225d3/2.048d3

  integer::imini,imaxi,jmini,jmaxi,kmini,kmaxi
  integer::imino,imaxo,jmino,jmaxo,kmino,kmaxo
  logical::decide3d


  real*8,dimension(3) :: CD,FD
  
  if(wei.ne.3)then
     write(*,*)"prolongrestrict.f90::restrict3: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

  CD = (uubc-llbc)/extc
  FD = (uubf-llbf)/extf

!take care the mismatch of the two segments of grid
  do i=1,3
     if(llbc(i) <= llbf(i))then
        base(i) = llbc(i)
     else
        j=idint((llbc(i)-llbf(i))/FD(i)+0.4)
        if(j/2*2 == j)then
           base(i) = llbf(i)
        else
           base(i) = llbf(i) - CD(i)/2
        endif
     endif
  enddo
!!! function idint:
!If A is of type REAL and |A| < 1, INT(A) equals 0. If |A| \geq 1, 
!then INT(A) equals the largest integer that does not exceed the range of A 
!and whose sign is the same as the sign of A.

! note say base = 0, llbf = 0, uubf = 2
! llbf->1 and uubf->2
    lbf = idint((llbf-base)/FD+0.4)+1
    ubf = idint((uubf-base)/FD+0.4)
    lbc = idint((llbc-base)/CD+0.4)+1
    ubc = idint((uubc-base)/CD+0.4)
    lbr = idint((llbr-base)/CD+0.4)+1
    lbrf = idint((llbr-base)/FD+0.4)+1
    ubr = idint((uubr-base)/CD+0.4)
    ubrf = idint((uubr-base)/FD+0.4)

!sanity check
  imino=lbr(1)-lbc(1) + 1
  imaxo=ubr(1)-lbc(1) + 1
  jmino=lbr(2)-lbc(2) + 1
  jmaxo=ubr(2)-lbc(2) + 1
  kmino=lbr(3)-lbc(3) + 1
  kmaxo=ubr(3)-lbc(3) + 1

  imini=lbrf(1)-lbf(1) + 1
  imaxi=ubrf(1)-lbf(1) + 1
  jmini=lbrf(2)-lbf(2) + 1
  jmaxi=ubrf(2)-lbf(2) + 1
  kmini=lbrf(3)-lbf(3) + 1
  kmaxi=ubrf(3)-lbf(3) + 1

  if(imino.lt.1.or.jmino.lt.1.or.kmino.lt.1.or.&
     imini.lt.1.or.jmini.lt.1.or.kmini.lt.1.or.&
     imaxo.gt.extc(1).or.jmaxo.gt.extc(2).or.kmaxo.gt.extc(3).or.&
     imaxi.gt.extf(1)-3.or.jmaxi.gt.extf(2)-3.or.kmaxi.gt.extf(3)-3)then
!-3 is because
!|-x---x-|-x---x-|-x---
!|- -*- -|
          write(*,*)"error in restrict for"
          write(*,*)"from"
          write(*,*)lbf,ubf
          write(*,*)"to"
          write(*,*)lbc,ubc
          write(*,*)"want"
          write(*,*)lbr,ubr,lbrf,ubrf
          write(*,*)"llbf = ",llbf
          write(*,*)"uubf = ",uubf
          write(*,*)"llbc = ",llbc
          write(*,*)"uubc = ",uubc
          write(*,*)"llbr = ",llbr
          write(*,*)"uubr = ",uubr
          write(*,*)"base = ",base
          stop
  endif

!~~~~~~> restriction start...
  do k = kmino,kmaxo
   do j = jmino,jmaxo
    do i = imino,imaxo

       cxI(1) = i
       cxI(2) = j
       cxI(3) = k
! change to fine level reference
!|---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*---| 
!|=======x===============x===============x===============x=======|
       cxI = 2*(cxI+lbc-1) - 1
! change to array index      
       cxI = cxI - lbf + 1
       if(any(cxI+4 > extf)) write(*,*)"error in restrict"
!  due to ghost zone, we can deal with symmetry boundary like this                   
       if(cxI(1)>3.and.cxI(2)>3.and.cxI(3)>3)then
          tmp2= C1*(funf(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-3)+funf(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+4))&
               +C2*(funf(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-2)+funf(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+3))&
               +C3*(funf(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-1)+funf(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+2))&
               +C4*(funf(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)  )+funf(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+1))
       else
          cxB=cxI-3
          cxT=cxI+4
          if(decide3d(extf,funf,funf,cxB,cxT,SoA,ya,8,Symmetry))then
             write(*,*)"restrict3 position index: ",i+lbc(1)-1,j+lbc(2)-1,k+lbc(3)-1
             stop
          endif
          tmp2=C1*(ya(:,:,1)+ya(:,:,8))+C2*(ya(:,:,2)+ya(:,:,7))+C3*(ya(:,:,3)+ya(:,:,6))+C4*(ya(:,:,4)+ya(:,:,5))
       endif
       tmp1= C1*(tmp2(:,1)+tmp2(:,8))+C2*(tmp2(:,2)+tmp2(:,7))+C3*(tmp2(:,3)+tmp2(:,6))+C4*(tmp2(:,4)+tmp2(:,5))
       func(i,j,k)= C1*(tmp1(1)+tmp1(8))+C2*(tmp1(2)+tmp1(7))+C3*(tmp1(3)+tmp1(6))+C4*(tmp1(4)+tmp1(5))
    enddo
   enddo
  enddo
  
  return

  end subroutine restrict3
!--------------------------------------------------------------------------
!
! Prolongation from coarser grids to finer grids
! 8 points, 7th order interpolation
! 1   2   3   4   5   6   7   8
! *---*---*---*---*---*---*---*
!              ^
! f=-495/262144*f_1 + 5005/262144*f_2 - 27027/262144*f_3 + 225225/262144*f_4
!   -429/262144*f_8 + 4095/262144*f_7 - 19305/262144*f_6 + 75075/262144*f_5
!--------------------------------------------------------------------------
  subroutine prolong3(wei,llbc,uubc,extc,func,&
                      llbf,uubf,extf,funf,&
                      llbp,uubp,SoA,Symmetry)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: wei
!                                       coarse    fine       coarse
  real*8,dimension(3),   intent(in) :: llbc,uubc,llbf,uubf,llbp,uubp
  integer,dimension(3),  intent(in) :: extc,extf
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in)   :: func
  real*8, dimension(extf(1),extf(2),extf(3)),intent(inout):: funf
  real*8, dimension(1:3), intent(in) :: SoA
  integer,intent(in)::Symmetry

!~~~~~~> local variables

  real*8, dimension(1:3) :: base
  integer,dimension(3) :: lbc,ubc,lbf,ubf,lbp,ubp,lbpc,ubpc
  integer,dimension(3) :: cxB,cxT,cxI
  integer :: i,j,k,ii,jj,kk
  real*8,  dimension(8,8,8) :: ya
  real*8, dimension(8,8) :: tmp2
  real*8, dimension(8) :: tmp1

  real*8, parameter :: C1=-4.95d2/2.62144d5,C2=5.005d3/2.62144d5,C3=-2.7027d4/2.62144d5,C4=2.25225d5/2.62144d5
  real*8, parameter :: C8=-4.29d2/2.62144d5,C7=4.095d3/2.62144d5,C6=-1.9305d4/2.62144d5,C5=7.5075d4/2.62144d5

  integer::imini,imaxi,jmini,jmaxi,kmini,kmaxi
  integer::imino,imaxo,jmino,jmaxo,kmino,kmaxo
  logical::decide3d

  real*8,dimension(3) :: CD,FD
  
  if(wei.ne.3)then
     write(*,*)"prolongrestrict.f90::prolong3: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

  CD = (uubc-llbc)/extc
  FD = (uubf-llbf)/extf

!take care the mismatch of the two segments of grid
  do i=1,3
     if(llbc(i) <= llbf(i))then
        base(i) = llbc(i)
     else
        j=idint((llbc(i)-llbf(i))/FD(i)+0.4)
        if(j/2*2 == j)then
           base(i) = llbf(i)
        else
           base(i) = llbf(i) - CD(i)/2
        endif
     endif
  enddo

!!! function idint:
!If A is of type REAL and |A| < 1, INT(A) equals 0. If |A| \geq 1, 
!then INT(A) equals the largest integer that does not exceed the range of A 
!and whose sign is the same as the sign of A.

    lbf = idint((llbf-base)/FD+0.4)+1
    ubf = idint((uubf-base)/FD+0.4)
    lbc = idint((llbc-base)/CD+0.4)+1
    ubc = idint((uubc-base)/CD+0.4)
    lbp = idint((llbp-base)/FD+0.4)+1
    lbpc = idint((llbp-base)/CD+0.4)+1
    ubp = idint((uubp-base)/FD+0.4)
    ubpc = idint((uubp-base)/CD+0.4)

!sanity check
  imino=lbp(1)-lbf(1) + 1
  imaxo=ubp(1)-lbf(1) + 1
  jmino=lbp(2)-lbf(2) + 1
  jmaxo=ubp(2)-lbf(2) + 1
  kmino=lbp(3)-lbf(3) + 1
  kmaxo=ubp(3)-lbf(3) + 1

  imini=lbpc(1)-lbc(1) + 1
  imaxi=ubpc(1)-lbc(1) + 1
  jmini=lbpc(2)-lbc(2) + 1
  jmaxi=ubpc(2)-lbc(2) + 1
  kmini=lbpc(3)-lbc(3) + 1
  kmaxi=ubpc(3)-lbc(3) + 1

  if(imino.lt.1.or.jmino.lt.1.or.kmino.lt.1.or.&
     imini.lt.1.or.jmini.lt.1.or.kmini.lt.1.or.&
     imaxo.gt.extf(1).or.jmaxo.gt.extf(2).or.kmaxo.gt.extf(3).or.&
     imaxi.gt.extc(1)-3.or.jmaxi.gt.extc(2)-3.or.kmaxi.gt.extc(3)-3)then
          write(*,*)"error in prolongation for"
          write(*,*)"from"
          write(*,*)llbc,uubc
          write(*,*)lbc,ubc
          write(*,*)"to"
          write(*,*)llbf,uubf
          write(*,*)lbf,ubf
          write(*,*)"want"
          write(*,*)llbp,uubp
          write(*,*)lbp,ubp,lbpc,ubpc
          return
  endif

!~~~~~~> prolongation start...
  do k = kmino,kmaxo
   do j = jmino,jmaxo
    do i = imino,imaxo
       cxI(1) = i
       cxI(2) = j
       cxI(3) = k
! change to coarse level reference
!|---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*---| 
!|=======x===============x===============x===============x=======|
       cxI = (cxI+lbf-1)/2
! change to array index      
       cxI = cxI - lbc + 1

       if(any(cxI+4 > extc)) write(*,*)"error in prolong"
       ii=i+lbf(1)-1
       jj=j+lbf(2)-1
       kk=k+lbf(3)-1
       if(ii/2*2==ii)then
         if(jj/2*2==jj)then
           if(kk/2*2==kk)then
!  due to ghost zone, we can deal with symmetry boundary like this                   
             if(cxI(1)>3.and.cxI(2)>3.and.cxI(3)>3)then
                tmp2= C1*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-3)+&
                      C2*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-2)+&
                      C3*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-1)+&
                      C4*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)  )+&
                      C5*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+1)+&
                      C6*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+2)+&
                      C7*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+3)+&
                      C8*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+4)
             else
                cxB=cxI-3
                cxT=cxI+4
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,8,Symmetry))then
                 write(*,*)"prolong3 position: "
                 write(*,*)llbf(1)+(i-0.5)*FD(1),llbf(2)+(j-0.5)*FD(2),llbf(3)+(k-0.5)*FD(3)
                 write(*,*)"llbf = ",llbf
                 stop
                endif
                tmp2= C1*ya(:,:,1)+C2*ya(:,:,2)+C3*ya(:,:,3)+C4*ya(:,:,4)+& 
                      C5*ya(:,:,5)+C6*ya(:,:,6)+C7*ya(:,:,7)+C8*ya(:,:,8)
             endif
             tmp1= C1*tmp2(:,1)+C2*tmp2(:,2)+C3*tmp2(:,3)+C4*tmp2(:,4)+&
                   C5*tmp2(:,5)+C6*tmp2(:,6)+C7*tmp2(:,7)+C8*tmp2(:,8)
             funf(i,j,k)= C1*tmp1(1)+C2*tmp1(2)+C3*tmp1(3)+C4*tmp1(4)+&
                          C5*tmp1(5)+C6*tmp1(6)+C7*tmp1(7)+C8*tmp1(8)
           else
             if(cxI(1)>3.and.cxI(2)>3.and.cxI(3)>3)then
                tmp2= C8*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-3)+&
                      C7*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-2)+&
                      C6*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-1)+&
                      C5*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)  )+&
                      C4*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+1)+&
                      C3*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+2)+&
                      C2*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+3)+&
                      C1*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+4)
             else
                cxB=cxI-3
                cxT=cxI+4
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,8,Symmetry))then
                 write(*,*)"prolong3 position: "
                 write(*,*)llbf(1)+(i-0.5)*FD(1),llbf(2)+(j-0.5)*FD(2),llbf(3)+(k-0.5)*FD(3)
                 write(*,*)"llbf = ",llbf
                 stop
                endif
                tmp2= C8*ya(:,:,1)+C7*ya(:,:,2)+C6*ya(:,:,3)+C5*ya(:,:,4)+&
                      C4*ya(:,:,5)+C3*ya(:,:,6)+C2*ya(:,:,7)+C1*ya(:,:,8)
             endif
             tmp1= C1*tmp2(:,1)+C2*tmp2(:,2)+C3*tmp2(:,3)+C4*tmp2(:,4)+&
                   C5*tmp2(:,5)+C6*tmp2(:,6)+C7*tmp2(:,7)+C8*tmp2(:,8)
             funf(i,j,k)=  C1*tmp1(1)+C2*tmp1(2)+C3*tmp1(3)+C4*tmp1(4)+&
                           C5*tmp1(5)+C6*tmp1(6)+C7*tmp1(7)+C8*tmp1(8)
           endif
         else
           if(kk/2*2==kk)then
             if(cxI(1)>3.and.cxI(2)>3.and.cxI(3)>3)then
                tmp2= C1*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-3)+&
                      C2*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-2)+&
                      C3*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-1)+&
                      C4*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)  )+&
                      C5*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+1)+&
                      C6*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+2)+&
                      C7*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+3)+&
                      C8*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+4)
             else
                cxB=cxI-3
                cxT=cxI+4
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,8,Symmetry))then
                 write(*,*)"prolong3 position: "
                 write(*,*)llbf(1)+(i-0.5)*FD(1),llbf(2)+(j-0.5)*FD(2),llbf(3)+(k-0.5)*FD(3)
                 write(*,*)"llbf = ",llbf
                 stop
                endif
                tmp2= C1*ya(:,:,1)+C2*ya(:,:,2)+C3*ya(:,:,3)+C4*ya(:,:,4)+&
                      C5*ya(:,:,5)+C6*ya(:,:,6)+C7*ya(:,:,7)+C8*ya(:,:,8)
             endif
             tmp1= C8*tmp2(:,1)+C7*tmp2(:,2)+C6*tmp2(:,3)+C5*tmp2(:,4)+&
                   C4*tmp2(:,5)+C3*tmp2(:,6)+C2*tmp2(:,7)+C1*tmp2(:,8)
             funf(i,j,k)= C1*tmp1(1)+C2*tmp1(2)+C3*tmp1(3)+C4*tmp1(4)+&
                          C5*tmp1(5)+C6*tmp1(6)+C7*tmp1(7)+C8*tmp1(8)
           else
             if(cxI(1)>3.and.cxI(2)>3.and.cxI(3)>3)then
                tmp2= C8*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-3)+&
                      C7*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-2)+&
                      C6*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-1)+&
                      C5*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)  )+&
                      C4*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+1)+&
                      C3*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+2)+&
                      C2*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+3)+&
                      C1*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+4)
             else
                cxB=cxI-3
                cxT=cxI+4
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,8,Symmetry))then
                 write(*,*)"prolong3 position: "
                 write(*,*)llbf(1)+(i-0.5)*FD(1),llbf(2)+(j-0.5)*FD(2),llbf(3)+(k-0.5)*FD(3)
                 write(*,*)"llbf = ",llbf
                 stop
                endif
                tmp2= C8*ya(:,:,1)+C7*ya(:,:,2)+C6*ya(:,:,3)+C5*ya(:,:,4)+&
                      C4*ya(:,:,5)+C3*ya(:,:,6)+C2*ya(:,:,7)+C1*ya(:,:,8)
             endif
             tmp1= C8*tmp2(:,1)+C7*tmp2(:,2)+C6*tmp2(:,3)+C5*tmp2(:,4)+&
                   C4*tmp2(:,5)+C3*tmp2(:,6)+C2*tmp2(:,7)+C1*tmp2(:,8)
             funf(i,j,k)=  C1*tmp1(1)+C2*tmp1(2)+C3*tmp1(3)+C4*tmp1(4)+&
                           C5*tmp1(5)+C6*tmp1(6)+C7*tmp1(7)+C8*tmp1(8)
           endif
         endif
       else
         if(jj/2*2==jj)then
           if(kk/2*2==kk)then               
             if(cxI(1)>3.and.cxI(2)>3.and.cxI(3)>3)then
                tmp2= C1*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-3)+&
                      C2*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-2)+&
                      C3*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-1)+&
                      C4*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)  )+&
                      C5*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+1)+&
                      C6*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+2)+&
                      C7*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+3)+&
                      C8*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+4)
             else
                cxB=cxI-3
                cxT=cxI+4
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,8,Symmetry))then
                 write(*,*)"prolong3 position: "
                 write(*,*)llbf(1)+(i-0.5)*FD(1),llbf(2)+(j-0.5)*FD(2),llbf(3)+(k-0.5)*FD(3)
                 write(*,*)"llbf = ",llbf
                 stop
                endif
                tmp2= C1*ya(:,:,1)+C2*ya(:,:,2)+C3*ya(:,:,3)+C4*ya(:,:,4)+&
                      C5*ya(:,:,5)+C6*ya(:,:,6)+C7*ya(:,:,7)+C8*ya(:,:,8)
             endif
             tmp1= C1*tmp2(:,1)+C2*tmp2(:,2)+C3*tmp2(:,3)+C4*tmp2(:,4)+&
                   C5*tmp2(:,5)+C6*tmp2(:,6)+C7*tmp2(:,7)+C8*tmp2(:,8)
             funf(i,j,k)= C8*tmp1(1)+C7*tmp1(2)+C6*tmp1(3)+C5*tmp1(4)+&
                          C4*tmp1(5)+C3*tmp1(6)+C2*tmp1(7)+C1*tmp1(8)
           else
             if(cxI(1)>3.and.cxI(2)>3.and.cxI(3)>3)then
                tmp2= C8*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-3)+&
                      C7*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-2)+&
                      C6*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-1)+&
                      C5*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)  )+&
                      C4*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+1)+&
                      C3*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+2)+&
                      C2*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+3)+&
                      C1*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+4)
             else
                cxB=cxI-3
                cxT=cxI+4
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,8,Symmetry))then
                 write(*,*)"prolong3 position: "
                 write(*,*)llbf(1)+(i-0.5)*FD(1),llbf(2)+(j-0.5)*FD(2),llbf(3)+(k-0.5)*FD(3)
                 write(*,*)"llbf = ",llbf
                 stop
                endif
                tmp2= C8*ya(:,:,1)+C7*ya(:,:,2)+C6*ya(:,:,3)+C5*ya(:,:,4)+&
                      C4*ya(:,:,5)+C3*ya(:,:,6)+C2*ya(:,:,7)+C1*ya(:,:,8)
             endif
             tmp1= C1*tmp2(:,1)+C2*tmp2(:,2)+C3*tmp2(:,3)+C4*tmp2(:,4)+&
                   C5*tmp2(:,5)+C6*tmp2(:,6)+C7*tmp2(:,7)+C8*tmp2(:,8)
             funf(i,j,k)=  C8*tmp1(1)+C7*tmp1(2)+C6*tmp1(3)+C5*tmp1(4)+&
                           C4*tmp1(5)+C3*tmp1(6)+C2*tmp1(7)+C1*tmp1(8)
           endif
         else
           if(kk/2*2==kk)then
             if(cxI(1)>3.and.cxI(2)>3.and.cxI(3)>3)then
                tmp2= C1*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-3)+&
                      C2*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-2)+&
                      C3*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-1)+&
                      C4*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)  )+&
                      C5*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+1)+&
                      C6*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+2)+&
                      C7*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+3)+&
                      C8*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+4)
             else
                cxB=cxI-3
                cxT=cxI+4
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,8,Symmetry))then
                 write(*,*)"prolong3 position: "
                 write(*,*)llbf(1)+(i-0.5)*FD(1),llbf(2)+(j-0.5)*FD(2),llbf(3)+(k-0.5)*FD(3)
                 write(*,*)"llbf = ",llbf
                 stop
                endif
                tmp2= C1*ya(:,:,1)+C2*ya(:,:,2)+C3*ya(:,:,3)+C4*ya(:,:,4)+&
                      C5*ya(:,:,5)+C6*ya(:,:,6)+C7*ya(:,:,7)+C8*ya(:,:,8)
             endif
             tmp1= C8*tmp2(:,1)+C7*tmp2(:,2)+C6*tmp2(:,3)+C5*tmp2(:,4)+&
                   C4*tmp2(:,5)+C3*tmp2(:,6)+C2*tmp2(:,7)+C1*tmp2(:,8)
             funf(i,j,k)= C8*tmp1(1)+C7*tmp1(2)+C6*tmp1(3)+C5*tmp1(4)+&
                          C4*tmp1(5)+C3*tmp1(6)+C2*tmp1(7)+C1*tmp1(8)
           else
             if(cxI(1)>3.and.cxI(2)>3.and.cxI(3)>3)then
                tmp2= C8*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-3)+&
                      C7*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-2)+&
                      C6*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-1)+&
                      C5*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)  )+&
                      C4*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+1)+&
                      C3*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+2)+&
                      C2*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+3)+&
                      C1*func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+4)
             else
                cxB=cxI-3
                cxT=cxI+4
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,8,Symmetry))then
                 write(*,*)"prolong3 position: "
                 write(*,*)llbf(1)+(i-0.5)*FD(1),llbf(2)+(j-0.5)*FD(2),llbf(3)+(k-0.5)*FD(3)
                 write(*,*)"llbf = ",llbf
                 stop
                endif
                tmp2= C8*ya(:,:,1)+C7*ya(:,:,2)+C6*ya(:,:,3)+C5*ya(:,:,4)+&
                      C4*ya(:,:,5)+C3*ya(:,:,6)+C2*ya(:,:,7)+C1*ya(:,:,8)
             endif
             tmp1= C8*tmp2(:,1)+C7*tmp2(:,2)+C6*tmp2(:,3)+C5*tmp2(:,4)+&
                   C4*tmp2(:,5)+C3*tmp2(:,6)+C2*tmp2(:,7)+C1*tmp2(:,8)
             funf(i,j,k)=  C8*tmp1(1)+C7*tmp1(2)+C6*tmp1(3)+C5*tmp1(4)+&
                           C4*tmp1(5)+C3*tmp1(6)+C2*tmp1(7)+C1*tmp1(8)
           endif
         endif
       endif
    enddo
   enddo
  enddo

  return

  end subroutine prolong3
#else
#error Not define Vertex nor Cell
#endif  
#endif

#elif (ghost_width == 5)
! eighth order code
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
!--------------------------------------------------------------------------
!
! Prolongation from coarser grids to finer grids
! 10 points, 9th order interpolation
! 1   2   3   4   5   6   7   8   9   10
! *---*---*---*---*---*---*---*---*---*
!                   ^
! f=35/65536(f_1+f_10)-405/65536*(f_2+f_9) + 567/16384*(f_3+f_8) - 2205/16384*(f_4+f_7) + 19845/32768*(f_5+f_6)
!--------------------------------------------------------------------------

  subroutine prolong3(wei,llbc,uubc,extc,func,&
                      llbf,uubf,extf,funf,&
                      llbp,uubp,SoA,Symmetry)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: wei
!                                       coarse    fine       coarse
  real*8,dimension(3),   intent(in) :: llbc,uubc,llbf,uubf,llbp,uubp
  integer,dimension(3),  intent(in) :: extc,extf
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in)   :: func
  real*8, dimension(extf(1),extf(2),extf(3)),intent(inout):: funf
  real*8, dimension(1:3), intent(in) :: SoA
  integer,intent(in)::Symmetry

!~~~~~~> local variables

  real*8, dimension(1:3) :: base
  integer,dimension(3) :: lbc,ubc,lbf,ubf,lbp,ubp,lbpc,ubpc
  integer,dimension(3) :: cxB,cxT,cxI
  integer :: i,j,k,ii,jj,kk
  real*8,  dimension(10,10,10) :: ya
  real*8, dimension(10,10) :: tmp2
  real*8, dimension(10) :: tmp1

  real*8, parameter :: C1=3.5d1/6.5536d4,C2=-4.05d2/6.5536d4,C3=5.67d2/1.6384d4
  real*8, parameter :: C4=-2.205d3/1.6384d4,C5=1.9845d4/3.2768d4

  integer::imini,imaxi,jmini,jmaxi,kmini,kmaxi
  integer::imino,imaxo,jmino,jmaxo,kmino,kmaxo
  logical::decide3d

  real*8,dimension(3) :: CD,FD

  if(wei.ne.3)then
     write(*,*)"prolongrestrict.f90::prolong3: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

! it's possible a iolated point for target but not for source
  CD = (uubc-llbc)/(extc-1)
  FD = CD/2

!take care the mismatch of the two segments of grid
  do i=1,3
     if(llbc(i) <= llbf(i))then
        base(i) = llbc(i)
     else
        j=idint((llbc(i)-llbf(i))/FD(i)+0.4)
        if(j/2*2 == j)then
           base(i) = llbf(i)
        else
           base(i) = llbf(i) - CD(i)/2
        endif
     endif
  enddo

!!! function idint:
!If A is of type REAL and |A| < 1, INT(A) equals 0. If |A| \geq 1, 
!then INT(A) equals the largest integer that does not exceed the range of A 
!and whose sign is the same as the sign of A.

    lbf = idint((llbf-base)/FD+0.4)+1
    ubf = idint((uubf-base)/FD+0.4)+1
    lbc = idint((llbc-base)/CD+0.4)+1
    ubc = idint((uubc-base)/CD+0.4)+1
    lbp = idint((llbp-base)/FD+0.4)+1
    lbpc = idint((llbp-base)/CD+0.4)+1
    ubp = idint((uubp-base)/FD+0.4)+1
    ubpc = idint((uubp-base)/CD+0.4)+1
!sanity check
  imino=lbp(1)-lbf(1) + 1
  imaxo=ubp(1)-lbf(1) + 1
  jmino=lbp(2)-lbf(2) + 1
  jmaxo=ubp(2)-lbf(2) + 1
  kmino=lbp(3)-lbf(3) + 1
  kmaxo=ubp(3)-lbf(3) + 1

  imini=lbpc(1)-lbc(1) + 1
  imaxi=ubpc(1)-lbc(1) + 1
  jmini=lbpc(2)-lbc(2) + 1
  jmaxi=ubpc(2)-lbc(2) + 1
  kmini=lbpc(3)-lbc(3) + 1
  kmaxi=ubpc(3)-lbc(3) + 1

  if(imino.lt.1.or.jmino.lt.1.or.kmino.lt.1.or.&
     imini.lt.1.or.jmini.lt.1.or.kmini.lt.1.or.&
     imaxo.gt.extf(1).or.jmaxo.gt.extf(2).or.kmaxo.gt.extf(3).or.&
     imaxi.gt.extc(1)-4.or.jmaxi.gt.extc(2)-4.or.kmaxi.gt.extc(3)-4)then
          write(*,*)"error in prolongation for"
          write(*,*)"from"
          write(*,*)llbc,uubc
          write(*,*)lbc,ubc
          write(*,*)"to"
          write(*,*)llbf,uubf
          write(*,*)lbf,ubf
          write(*,*)"want"
          write(*,*)llbp,uubp
          write(*,*)lbp,ubp,lbpc,ubpc
          if(imini.lt.1) write(*,*)"imini = ",imini
          if(jmini.lt.1) write(*,*)"jmini = ",jmini       
          if(kmini.lt.1) write(*,*)"kmini = ",kmini      
          if(imino.lt.1) write(*,*)"imino = ",imino       
          if(jmino.lt.1) write(*,*)"jmino = ",jmino       
          if(kmino.lt.1) write(*,*)"kmino = ",kmino   
          if(imaxi.gt.extc(1)) write(*,*)"imaxi = ",imaxi,"extc(1) = ",extc(1) 
          if(jmaxi.gt.extc(2)) write(*,*)"jmaxi = ",jmaxi,"extc(2) = ",extc(2) 
          if(kmaxi.gt.extc(3)) write(*,*)"kmaxi = ",kmaxi,"extc(3) = ",extc(3) 
          if(imaxo.gt.extf(1)) write(*,*)"imaxo = ",imaxo,"extf(1) = ",extf(1) 
          if(jmaxo.gt.extf(2)) write(*,*)"jmaxo = ",jmaxo,"extf(2) = ",extf(2) 
          if(kmaxo.gt.extf(3)) write(*,*)"kmaxo = ",kmaxo,"extf(3) = ",extf(3) 
          return
  endif
!~~~~~~> prolongation start...
  do k = kmino,kmaxo
   do j = jmino,jmaxo
    do i = imino,imaxo
! change to coarse level reference
!|*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*| 
!|x===============x===============x===============x========|
       cxI(1) = i
       cxI(2) = j
       cxI(3) = k
       cxI = (cxI+lbf)/2
! change to array index      
       cxI = cxI - lbc + 1

       ii=i+lbf(1)-1
       jj=j+lbf(2)-1
       kk=k+lbf(3)-1

       if(any(cxI+5 > extc)) write(*,*)"error in prolong"
       if(ii/2*2==ii)then
         if(jj/2*2==jj)then
           if(kk/2*2==kk)then
!  due to ghost zone, we can deal with symmetry boundary like this                                   
             if(cxI(1)>4.and.cxI(2)>4.and.cxI(3)>4)then
               tmp2= C1*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-4)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+5))&
                    +C2*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-3)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+4))&
                    +C3*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-2)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+3))&
                    +C4*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-1)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+2))&
                    +C5*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)  )+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+1))
             else
               cxB=cxI-4
               cxT=cxI+5
               if(decide3d(extc,func,func,cxB,cxT,SoA,ya,10,Symmetry))then
                write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                stop
               endif
               tmp2=C1*(ya(:,:,1)+ya(:,:,10))+C2*(ya(:,:,2)+ya(:,:,9))+C3*(ya(:,:,3)+ya(:,:,8)) &
                   +C4*(ya(:,:,4)+ya(:,:, 7))+C5*(ya(:,:,5)+ya(:,:,6))
             endif
             tmp1= C1*(tmp2(:,1)+tmp2(:,10))+C2*(tmp2(:,2)+tmp2(:,9))+C3*(tmp2(:,3)+tmp2(:,8)) &
                  +C4*(tmp2(:,4)+tmp2(:, 7))+C5*(tmp2(:,5)+tmp2(:,6))
             funf(i,j,k)= C1*(tmp1(1)+tmp1(10))+C2*(tmp1(2)+tmp1(9))+C3*(tmp1(3)+tmp1(8)) &
                         +C4*(tmp1(4)+tmp1( 7))+C5*(tmp1(5)+tmp1(6))
           else
             if(cxI(1)>4.and.cxI(2)>4.and.cxI(3)>4)then
                tmp2= func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3))
             else
                cxB=cxI-4
                cxT=cxI+5
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,10,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= ya(:,:,5)
             endif
             tmp1= C1*(tmp2(:,1)+tmp2(:,10))+C2*(tmp2(:,2)+tmp2(:,9))+C3*(tmp2(:,3)+tmp2(:,8)) &
                  +C4*(tmp2(:,4)+tmp2(:, 7))+C5*(tmp2(:,5)+tmp2(:,6))
             funf(i,j,k)= C1*(tmp1(1)+tmp1(10))+C2*(tmp1(2)+tmp1(9))+C3*(tmp1(3)+tmp1(8)) &
                         +C4*(tmp1(4)+tmp1( 7))+C5*(tmp1(5)+tmp1(6))
           endif
         else
           if(kk/2*2==kk)then
             if(cxI(1)>4.and.cxI(2)>4.and.cxI(3)>4)then
               tmp2= C1*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-4)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+5))&
                    +C2*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-3)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+4))&
                    +C3*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-2)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+3))&
                    +C4*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-1)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+2))&
                    +C5*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)  )+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+1))
             else
                cxB=cxI-4
                cxT=cxI+5
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,10,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
               tmp2=C1*(ya(:,:,1)+ya(:,:,10))+C2*(ya(:,:,2)+ya(:,:,9))+C3*(ya(:,:,3)+ya(:,:,8)) &
                   +C4*(ya(:,:,4)+ya(:,:, 7))+C5*(ya(:,:,5)+ya(:,:,6))
             endif
             tmp1= tmp2(:,5)
             funf(i,j,k)= C1*(tmp1(1)+tmp1(10))+C2*(tmp1(2)+tmp1(9))+C3*(tmp1(3)+tmp1(8)) &
                         +C4*(tmp1(4)+tmp1( 7))+C5*(tmp1(5)+tmp1(6))
           else
             if(cxI(1)>4.and.cxI(2)>4.and.cxI(3)>4)then
                tmp2= func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3))
             else
                cxB=cxI-4
                cxT=cxI+5
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,10,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= ya(:,:,5)
             endif
             tmp1= tmp2(:,5)
             funf(i,j,k)= C1*(tmp1(1)+tmp1(10))+C2*(tmp1(2)+tmp1(9))+C3*(tmp1(3)+tmp1(8)) &
                         +C4*(tmp1(4)+tmp1( 7))+C5*(tmp1(5)+tmp1(6))
           endif
         endif
       else
         if(jj/2*2==jj)then
           if(kk/2*2==kk)then               
             if(cxI(1)>4.and.cxI(2)>4.and.cxI(3)>4)then
               tmp2= C1*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-4)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+5))&
                    +C2*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-3)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+4))&
                    +C3*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-2)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+3))&
                    +C4*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-1)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+2))&
                    +C5*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)  )+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+1))
             else
                cxB=cxI-4
                cxT=cxI+5
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,10,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
               tmp2=C1*(ya(:,:,1)+ya(:,:,10))+C2*(ya(:,:,2)+ya(:,:,9))+C3*(ya(:,:,3)+ya(:,:,8)) &
                   +C4*(ya(:,:,4)+ya(:,:, 7))+C5*(ya(:,:,5)+ya(:,:,6))
             endif
             tmp1= C1*(tmp2(:,1)+tmp2(:,10))+C2*(tmp2(:,2)+tmp2(:,9))+C3*(tmp2(:,3)+tmp2(:,8)) &
                  +C4*(tmp2(:,4)+tmp2(:, 7))+C5*(tmp2(:,5)+tmp2(:,6))
             funf(i,j,k)= tmp1(5)
           else
             if(cxI(1)>4.and.cxI(2)>4.and.cxI(3)>4)then
                tmp2= func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3))
             else
                cxB=cxI-4
                cxT=cxI+5
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,10,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= ya(:,:,5)
             endif
             tmp1= C1*(tmp2(:,1)+tmp2(:,10))+C2*(tmp2(:,2)+tmp2(:,9))+C3*(tmp2(:,3)+tmp2(:,8)) &
                  +C4*(tmp2(:,4)+tmp2(:, 7))+C5*(tmp2(:,5)+tmp2(:,6))
             funf(i,j,k)= tmp1(5)
           endif
         else
           if(kk/2*2==kk)then
             if(cxI(1)>4.and.cxI(2)>4.and.cxI(3)>4)then
               tmp2= C1*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-4)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+5))&
                    +C2*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-3)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+4))&
                    +C3*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-2)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+3))&
                    +C4*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-1)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+2))&
                    +C5*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)  )+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+1))
             else
                cxB=cxI-4
                cxT=cxI+5
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,10,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
               tmp2=C1*(ya(:,:,1)+ya(:,:,10))+C2*(ya(:,:,2)+ya(:,:,9))+C3*(ya(:,:,3)+ya(:,:,8)) &
                   +C4*(ya(:,:,4)+ya(:,:, 7))+C5*(ya(:,:,5)+ya(:,:,6))
             endif
             tmp1= tmp2(:,5)
             funf(i,j,k)= tmp1(5)
           else
             funf(i,j,k)= func(cxI(1),cxI(2),cxI(3))
           endif
         endif
       endif
    enddo
   enddo
  enddo

  return

  end subroutine prolong3

#else
#ifdef Cell
!---------------------------------------------------------------------------------------------------------------
!
! Restrict from finner grids to coarser grids ignore the boundary point
!
! 10 points, 9th order interpolation
! 1   2   3   4   5   6   7   8   9   10
! *---*---*---*---*---*---*---*---*---*
!                   ^
! f=35/65536(f_1+f_10)-405/65536*(f_2+f_9) + 567/16384*(f_3+f_8) - 2205/16384*(f_4+f_7) + 19845/32768*(f_5+f_6)
!---------------------------------------------------------------------------------------------------------------

  subroutine restrict3(wei,llbc,uubc,extc,func,&
                       llbf,uubf,extf,funf,&
                       llbr,uubr,SoA,Symmetry)
  implicit none

!~~~~~~> input arguments
  integer,intent(in)::wei
!                                       coarse    fine       coarse
  real*8,dimension(3),   intent(in) :: llbc,uubc,llbf,uubf,llbr,uubr
  integer,dimension(3),  intent(in) :: extc,extf
  real*8, dimension(extc(1),extc(2),extc(3)),intent(inout):: func
  real*8, dimension(extf(1),extf(2),extf(3)),intent(in):: funf
  real*8, dimension(1:3), intent(in) :: SoA
  integer,intent(in)::Symmetry

!~~~~~~> local variables

  real*8, dimension(1:3) :: base
  integer,dimension(3) :: lbc,ubc,lbf,ubf,lbr,ubr,lbrf,ubrf
  integer,dimension(3) :: cxB,cxT,cxI
  integer :: i,j,k
  real*8,  dimension(10,10,10) :: ya
  real*8, dimension(10,10) :: tmp2
  real*8, dimension(10) :: tmp1
  real*8, parameter :: C1=3.5d1/6.5536d4,C2=-4.05d2/6.5536d4,C3=5.67d2/1.6384d4
  real*8, parameter :: C4=-2.205d3/1.6384d4,C5=1.9845d4/3.2768d4

  integer::imini,imaxi,jmini,jmaxi,kmini,kmaxi
  integer::imino,imaxo,jmino,jmaxo,kmino,kmaxo
  logical::decide3d

  real*8,dimension(3) :: CD,FD
  
  if(wei.ne.3)then
     write(*,*)"prolongrestrict.f90::restrict3: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

  CD = (uubc-llbc)/extc
  FD = (uubf-llbf)/extf

!take care the mismatch of the two segments of grid
  do i=1,3
     if(llbc(i) <= llbf(i))then
        base(i) = llbc(i)
     else
        j=idint((llbc(i)-llbf(i))/FD(i)+0.4)
        if(j/2*2 == j)then
           base(i) = llbf(i)
        else
           base(i) = llbf(i) - CD(i)/2
        endif
     endif
  enddo
!!! function idint:
!If A is of type REAL and |A| < 1, INT(A) equals 0. If |A| \geq 1, 
!then INT(A) equals the largest integer that does not exceed the range of A 
!and whose sign is the same as the sign of A.

! note say base = 0, llbf = 0, uubf = 2
! llbf->1 and uubf->2
    lbf = idint((llbf-base)/FD+0.4)+1
    ubf = idint((uubf-base)/FD+0.4)
    lbc = idint((llbc-base)/CD+0.4)+1
    ubc = idint((uubc-base)/CD+0.4)
    lbr = idint((llbr-base)/CD+0.4)+1
    lbrf = idint((llbr-base)/FD+0.4)+1
    ubr = idint((uubr-base)/CD+0.4)
    ubrf = idint((uubr-base)/FD+0.4)

!sanity check
  imino=lbr(1)-lbc(1) + 1
  imaxo=ubr(1)-lbc(1) + 1
  jmino=lbr(2)-lbc(2) + 1
  jmaxo=ubr(2)-lbc(2) + 1
  kmino=lbr(3)-lbc(3) + 1
  kmaxo=ubr(3)-lbc(3) + 1

  imini=lbrf(1)-lbf(1) + 1
  imaxi=ubrf(1)-lbf(1) + 1
  jmini=lbrf(2)-lbf(2) + 1
  jmaxi=ubrf(2)-lbf(2) + 1
  kmini=lbrf(3)-lbf(3) + 1
  kmaxi=ubrf(3)-lbf(3) + 1

  if(imino.lt.1.or.jmino.lt.1.or.kmino.lt.1.or.&
     imini.lt.1.or.jmini.lt.1.or.kmini.lt.1.or.&
     imaxo.gt.extc(1).or.jmaxo.gt.extc(2).or.kmaxo.gt.extc(3).or.&
     imaxi.gt.extf(1)-4.or.jmaxi.gt.extf(2)-4.or.kmaxi.gt.extf(3)-4)then
          write(*,*)"error in restrict for"
          write(*,*)"from"
          write(*,*)lbf,ubf
          write(*,*)"to"
          write(*,*)lbc,ubc
          write(*,*)"want"
          write(*,*)lbr,ubr,lbrf,ubrf
          write(*,*)"llbf = ",llbf
          write(*,*)"uubf = ",uubf
          write(*,*)"llbc = ",llbc
          write(*,*)"uubc = ",uubc
          write(*,*)"llbr = ",llbr
          write(*,*)"uubr = ",uubr
          write(*,*)"base = ",base
          stop
  endif

!~~~~~~> restriction start...
  do k = kmino,kmaxo
   do j = jmino,jmaxo
    do i = imino,imaxo

       cxI(1) = i
       cxI(2) = j
       cxI(3) = k
! change to fine level reference
!|---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*---| 
!|=======x===============x===============x===============x=======|
       cxI = 2*(cxI+lbc-1) - 1
! change to array index      
       cxI = cxI - lbf + 1

       if(any(cxI+5 > extf)) write(*,*)"error in restrict"
!  due to ghost zone, we can deal with symmetry boundary like this                   
       if(cxI(1)>4.and.cxI(2)>4.and.cxI(3)>4)then
          tmp2= C1*(funf(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-4)+funf(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+5))&
               +C2*(funf(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-3)+funf(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+4))&
               +C3*(funf(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-2)+funf(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+3))&
               +C4*(funf(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-1)+funf(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+2))&
               +C5*(funf(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)  )+funf(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+1))
       else
          cxB=cxI-4
          cxT=cxI+5
          if(decide3d(extf,funf,funf,cxB,cxT,SoA,ya,10,Symmetry))then
             write(*,*)"restrict3 position index: ",i+lbc(1)-1,j+lbc(2)-1,k+lbc(3)-1
             stop
          endif
          tmp2=C1*(ya(:,:,1)+ya(:,:,10))+C2*(ya(:,:,2)+ya(:,:,9))+C3*(ya(:,:,3)+ya(:,:,8)) &
              +C4*(ya(:,:,4)+ya(:,:, 7))+C5*(ya(:,:,5)+ya(:,:,6))
       endif
       tmp1= C1*(tmp2(:,1)+tmp2(:,10))+C2*(tmp2(:,2)+tmp2(:,9))+C3*(tmp2(:,3)+tmp2(:,8)) &
            +C4*(tmp2(:,4)+tmp2(:, 7))+C5*(tmp2(:,5)+tmp2(:,6))
       func(i,j,k)= C1*(tmp1(1)+tmp1(10))+C2*(tmp1(2)+tmp1(9))+C3*(tmp1(3)+tmp1(8)) &
                   +C4*(tmp1(4)+tmp1( 7))+C5*(tmp1(5)+tmp1(6))
    enddo
   enddo
  enddo
  
  return

  end subroutine restrict3
!--------------------------------------------------------------------------
!
! Prolongation from coarser grids to finer grids
! 10 points, 9th order interpolation
! 1   2   3   4   5   6   7   8   9   10
! *---*---*---*---*---*---*---*---*---*
!                  ^
!f=  13585/33554432*f_1-159885/33554432*f_2+230945/8388608*f_3- 969969/8388608*f_4+14549535/16777216*f_5
! +4849845/16777216*f_6- 692835/8388608*f_7+188955/8388608*f_8-138567/33554432*f_9+   12155/33554432*f_10
!--------------------------------------------------------------------------

  subroutine prolong3(wei,llbc,uubc,extc,func,&
                      llbf,uubf,extf,funf,&
                      llbp,uubp,SoA,Symmetry)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: wei
!                                       coarse    fine       coarse
  real*8,dimension(3),   intent(in) :: llbc,uubc,llbf,uubf,llbp,uubp
  integer,dimension(3),  intent(in) :: extc,extf
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in)   :: func
  real*8, dimension(extf(1),extf(2),extf(3)),intent(inout):: funf
  real*8, dimension(1:3), intent(in) :: SoA
  integer,intent(in)::Symmetry

!~~~~~~> local variables

  real*8, dimension(1:3) :: base
  integer,dimension(3) :: lbc,ubc,lbf,ubf,lbp,ubp,lbpc,ubpc
  integer,dimension(3) :: cxB,cxT,cxI
  integer :: i,j,k,ii,jj,kk
  real*8,  dimension(10,10,10) :: ya
  real*8, dimension(10,10) :: tmp2
  real*8, dimension(10) :: tmp1

  real*8, parameter :: C1=1.3585d4/3.3554432d7,C2=-1.59885d5/3.3554432d7,C3=2.30945d5/8.388608d6
  real*8, parameter :: C4=-9.69969d5/8.388608d6,C5=1.4549535d7/1.6777216d7,C6=4.849845d6/1.6777216d7
  real*8, parameter :: C7=-6.92835d5/8.388608d6,C8=1.88955d5/8.388608d6,C9=-1.38567d5/3.3554432d7
  real*8, parameter :: C10=1.2155d4/3.3554432d7

  integer::imini,imaxi,jmini,jmaxi,kmini,kmaxi
  integer::imino,imaxo,jmino,jmaxo,kmino,kmaxo
  logical::decide3d

  real*8,dimension(3) :: CD,FD
  
  if(wei.ne.3)then
     write(*,*)"prolongrestrict.f90::prolong3: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

  CD = (uubc-llbc)/extc
  FD = (uubf-llbf)/extf

!take care the mismatch of the two segments of grid
  do i=1,3
     if(llbc(i) <= llbf(i))then
        base(i) = llbc(i)
     else
        j=idint((llbc(i)-llbf(i))/FD(i)+0.4)
        if(j/2*2 == j)then
           base(i) = llbf(i)
        else
           base(i) = llbf(i) - CD(i)/2
        endif
     endif
  enddo

!!! function idint:
!If A is of type REAL and |A| < 1, INT(A) equals 0. If |A| \geq 1, 
!then INT(A) equals the largest integer that does not exceed the range of A 
!and whose sign is the same as the sign of A.

    lbf = idint((llbf-base)/FD+0.4)+1
    ubf = idint((uubf-base)/FD+0.4)
    lbc = idint((llbc-base)/CD+0.4)+1
    ubc = idint((uubc-base)/CD+0.4)
    lbp = idint((llbp-base)/FD+0.4)+1
    lbpc = idint((llbp-base)/CD+0.4)+1
    ubp = idint((uubp-base)/FD+0.4)
    ubpc = idint((uubp-base)/CD+0.4)

!sanity check
  imino=lbp(1)-lbf(1) + 1
  imaxo=ubp(1)-lbf(1) + 1
  jmino=lbp(2)-lbf(2) + 1
  jmaxo=ubp(2)-lbf(2) + 1
  kmino=lbp(3)-lbf(3) + 1
  kmaxo=ubp(3)-lbf(3) + 1

  imini=lbpc(1)-lbc(1) + 1
  imaxi=ubpc(1)-lbc(1) + 1
  jmini=lbpc(2)-lbc(2) + 1
  jmaxi=ubpc(2)-lbc(2) + 1
  kmini=lbpc(3)-lbc(3) + 1
  kmaxi=ubpc(3)-lbc(3) + 1

  if(imino.lt.1.or.jmino.lt.1.or.kmino.lt.1.or.&
     imini.lt.1.or.jmini.lt.1.or.kmini.lt.1.or.&
     imaxo.gt.extf(1).or.jmaxo.gt.extf(2).or.kmaxo.gt.extf(3).or.&
     imaxi.gt.extc(1)-4.or.jmaxi.gt.extc(2)-4.or.kmaxi.gt.extc(3)-4)then
          write(*,*)"error in prolongation for"
          write(*,*)"from"
          write(*,*)llbc,uubc
          write(*,*)lbc,ubc
          write(*,*)"to"
          write(*,*)llbf,uubf
          write(*,*)lbf,ubf
          write(*,*)"want"
          write(*,*)llbp,uubp
          write(*,*)lbp,ubp,lbpc,ubpc
          return
  endif

!~~~~~~> prolongation start...
  do k = kmino,kmaxo
   do j = jmino,jmaxo
    do i = imino,imaxo
       cxI(1) = i
       cxI(2) = j
       cxI(3) = k
! change to coarse level reference
!|---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*---| 
!|=======x===============x===============x===============x=======|
       cxI = (cxI+lbf-1)/2
! change to array index      
       cxI = cxI - lbc + 1

       if(any(cxI+5 > extc)) write(*,*)"error in prolong"
       ii=i+lbf(1)-1
       jj=j+lbf(2)-1
       kk=k+lbf(3)-1
       if(ii/2*2==ii)then
         if(jj/2*2==jj)then
           if(kk/2*2==kk)then
!  due to ghost zone, we can deal with symmetry boundary like this                   
             if(cxI(1)>4.and.cxI(2)>4.and.cxI(3)>4)then
                tmp2= C1 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-4)+&
                      C2 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-3)+&
                      C3 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-2)+&
                      C4 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-1)+&
                      C5 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)  )+&
                      C6 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+1)+&
                      C7 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+2)+&
                      C8 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+3)+&
                      C9 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+4)+&
                      C10*func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+5)
             else
                cxB=cxI-4
                cxT=cxI+5
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,10,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= C1*ya(:,:,1)+C2*ya(:,:,2)+C3*ya(:,:,3)+C4*ya(:,:,4)+C5 *ya(:,:, 5)+& 
                      C6*ya(:,:,6)+C7*ya(:,:,7)+C8*ya(:,:,8)+C9*ya(:,:,9)+C10*ya(:,:,10)
             endif
             tmp1= C1*tmp2(:,1)+C2*tmp2(:,2)+C3*tmp2(:,3)+C4*tmp2(:,4)+C5 *tmp2(:, 5)+&
                   C6*tmp2(:,6)+C7*tmp2(:,7)+C8*tmp2(:,8)+C9*tmp2(:,9)+C10*tmp2(:,10)
             funf(i,j,k)= C1*tmp1(1)+C2*tmp1(2)+C3*tmp1(3)+C4*tmp1(4)+ C5*tmp1( 5)+&
                          C6*tmp1(6)+C7*tmp1(7)+C8*tmp1(8)+C9*tmp1(9)+C10*tmp1(10)
           else
             if(cxI(1)>4.and.cxI(2)>4.and.cxI(3)>4)then
                tmp2= C10*func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-4)+&
                      C9 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-3)+&
                      C8 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-2)+&
                      C7 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-1)+&
                      C6 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)  )+&
                      C5 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+1)+&
                      C4 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+2)+&
                      C3 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+3)+&
                      C2 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+4)+&
                      C1 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+5)
             else
                cxB=cxI-4
                cxT=cxI+5
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,10,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= C10*ya(:,:,1)+C9*ya(:,:,2)+C8*ya(:,:,3)+C7*ya(:,:,4)+C6*ya(:,:, 5)+&
                      C5 *ya(:,:,6)+C4*ya(:,:,7)+C3*ya(:,:,8)+C2*ya(:,:,9)+C1*ya(:,:,10)
             endif
             tmp1= C1*tmp2(:,1)+C2*tmp2(:,2)+C3*tmp2(:,3)+C4*tmp2(:,4)+C5 *tmp2(:, 5)+&
                   C6*tmp2(:,6)+C7*tmp2(:,7)+C8*tmp2(:,8)+C9*tmp2(:,9)+C10*tmp2(:,10)
             funf(i,j,k)=  C1*tmp1(1)+C2*tmp1(2)+C3*tmp1(3)+C4*tmp1(4)+C5 *tmp1( 5)+&
                           C6*tmp1(6)+C7*tmp1(7)+C8*tmp1(8)+C9*tmp1(9)+C10*tmp1(10)
           endif
         else
           if(kk/2*2==kk)then
             if(cxI(1)>4.and.cxI(2)>4.and.cxI(3)>4)then
                tmp2= C1 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-4)+&
                      C2 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-3)+&
                      C3 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-2)+&
                      C4 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-1)+&
                      C5 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)  )+&
                      C6 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+1)+&
                      C7 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+2)+&
                      C8 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+3)+&
                      C9 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+4)+&
                      C10*func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+5)
             else
                cxB=cxI-4
                cxT=cxI+5
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,10,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= C1*ya(:,:,1)+C2*ya(:,:,2)+C3*ya(:,:,3)+C4*ya(:,:,4)+C5 *ya(:,:, 5)+&
                      C6*ya(:,:,6)+C7*ya(:,:,7)+C8*ya(:,:,8)+C9*ya(:,:,9)+C10*ya(:,:,10)
             endif
             tmp1= C10*tmp2(:,1)+C9*tmp2(:,2)+C8*tmp2(:,3)+C7*tmp2(:,4)+C6*tmp2(:, 5)+&
                   C5 *tmp2(:,6)+C4*tmp2(:,7)+C3*tmp2(:,8)+C2*tmp2(:,9)+C1*tmp2(:,10)
             funf(i,j,k)= C1*tmp1(1)+C2*tmp1(2)+C3*tmp1(3)+C4*tmp1(4)+C5 *tmp1( 5)+&
                          C6*tmp1(6)+C7*tmp1(7)+C8*tmp1(8)+C9*tmp1(9)+C10*tmp1(10)
           else
             if(cxI(1)>4.and.cxI(2)>4.and.cxI(3)>4)then
                tmp2= C10*func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-4)+&
                      C9 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-3)+&
                      C8 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-2)+&
                      C7 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-1)+&
                      C6 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)  )+&
                      C5 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+1)+&
                      C4 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+2)+&
                      C3 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+3)+&
                      C2 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+4)+&
                      C1 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+5)
             else
                cxB=cxI-4
                cxT=cxI+5
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,10,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= C10*ya(:,:,1)+C9*ya(:,:,2)+C8*ya(:,:,3)+C7*ya(:,:,4)+C6*ya(:,:, 5)+&
                      C5 *ya(:,:,6)+C4*ya(:,:,7)+C3*ya(:,:,8)+C2*ya(:,:,9)+C1*ya(:,:,10)
             endif
             tmp1= C10*tmp2(:,1)+C9*tmp2(:,2)+C8*tmp2(:,3)+C7*tmp2(:,4)+C6*tmp2(:, 5)+&
                   C5 *tmp2(:,6)+C4*tmp2(:,7)+C3*tmp2(:,8)+C2*tmp2(:,9)+C1*tmp2(:,10)
             funf(i,j,k)=  C1*tmp1(1)+C2*tmp1(2)+C3*tmp1(3)+C4*tmp1(4)+C5 *tmp1( 5)+&
                           C6*tmp1(6)+C7*tmp1(7)+C8*tmp1(8)+C9*tmp1(9)+C10*tmp1(10)
           endif
         endif
       else
         if(jj/2*2==jj)then
           if(kk/2*2==kk)then               
             if(cxI(1)>4.and.cxI(2)>4.and.cxI(3)>4)then
                tmp2= C1 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-4)+&
                      C2 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-3)+&
                      C3 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-2)+&
                      C4 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-1)+&
                      C5 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)  )+&
                      C6 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+1)+&
                      C7 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+2)+&
                      C8 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+3)+&
                      C9 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+4)+&
                      C10*func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+5)
             else
                cxB=cxI-4
                cxT=cxI+5
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,10,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= C1*ya(:,:,1)+C2*ya(:,:,2)+C3*ya(:,:,3)+C4*ya(:,:,4)+C5 *ya(:,:, 5)+&
                      C6*ya(:,:,6)+C7*ya(:,:,7)+C8*ya(:,:,8)+C9*ya(:,:,9)+C10*ya(:,:,10)
             endif
             tmp1= C1*tmp2(:,1)+C2*tmp2(:,2)+C3*tmp2(:,3)+C4*tmp2(:,4)+C5 *tmp2(:, 5)+&
                   C6*tmp2(:,6)+C7*tmp2(:,7)+C8*tmp2(:,8)+C9*tmp2(:,9)+C10*tmp2(:,10)
             funf(i,j,k)= C10*tmp1(1)+C9*tmp1(2)+C8*tmp1(3)+C7*tmp1(4)+C6*tmp1( 5)+&
                          C5 *tmp1(6)+C4*tmp1(7)+C3*tmp1(8)+C2*tmp1(9)+C1*tmp1(10)
           else
             if(cxI(1)>4.and.cxI(2)>4.and.cxI(3)>4)then
                tmp2= C10*func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-4)+&
                      C9 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-3)+&
                      C8 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-2)+&
                      C7 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-1)+&
                      C6 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)  )+&
                      C5 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+1)+&
                      C4 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+2)+&
                      C3 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+3)+&
                      C2 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+4)+&
                      C1 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+5)
             else
                cxB=cxI-4
                cxT=cxI+5
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,10,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= C10*ya(:,:,1)+C9*ya(:,:,2)+C8*ya(:,:,3)+C7*ya(:,:,4)+C6*ya(:,:, 5)+&
                      C5 *ya(:,:,6)+C4*ya(:,:,7)+C3*ya(:,:,8)+C2*ya(:,:,9)+C1*ya(:,:,10)
             endif
             tmp1= C1*tmp2(:,1)+C2*tmp2(:,2)+C3*tmp2(:,3)+C4*tmp2(:,4)+C5 *tmp2(:, 5)+&
                   C6*tmp2(:,6)+C7*tmp2(:,7)+C8*tmp2(:,8)+C9*tmp2(:,9)+C10*tmp2(:,10)
             funf(i,j,k)=  C10*tmp1(1)+C9*tmp1(2)+C8*tmp1(3)+C7*tmp1(4)+C6*tmp1( 5)+&
                           C5 *tmp1(6)+C4*tmp1(7)+C3*tmp1(8)+C2*tmp1(9)+C1*tmp1(10)
           endif
         else
           if(kk/2*2==kk)then
             if(cxI(1)>4.and.cxI(2)>4.and.cxI(3)>4)then
                tmp2= C1 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-4)+&
                      C2 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-3)+&
                      C3 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-2)+&
                      C4 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-1)+&
                      C5 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)  )+&
                      C6 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+1)+&
                      C7 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+2)+&
                      C8 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+3)+&
                      C9 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+4)+&
                      C10*func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+5)
             else
                cxB=cxI-4
                cxT=cxI+5
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,10,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= C1*ya(:,:,1)+C2*ya(:,:,2)+C3*ya(:,:,3)+C4*ya(:,:,4)+C5 *ya(:,:, 5)+&
                      C6*ya(:,:,6)+C7*ya(:,:,7)+C8*ya(:,:,8)+C9*ya(:,:,9)+C10*ya(:,:,10)
             endif
             tmp1= C10*tmp2(:,1)+C9*tmp2(:,2)+C8*tmp2(:,3)+C7*tmp2(:,4)+C6*tmp2(:, 5)+&
                   C5 *tmp2(:,6)+C4*tmp2(:,7)+C3*tmp2(:,8)+C2*tmp2(:,9)+C1*tmp2(:,10)
             funf(i,j,k)= C10*tmp1(1)+C9*tmp1(2)+C8*tmp1(3)+C7*tmp1(4)+C6*tmp1( 5)+&
                          C5 *tmp1(6)+C4*tmp1(7)+C3*tmp1(8)+C2*tmp1(9)+C1*tmp1(10)
           else
             if(cxI(1)>4.and.cxI(2)>4.and.cxI(3)>4)then
                tmp2= C10*func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-4)+&
                      C9 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-3)+&
                      C8 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-2)+&
                      C7 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-1)+&
                      C6 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)  )+&
                      C5 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+1)+&
                      C4 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+2)+&
                      C3 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+3)+&
                      C2 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+4)+&
                      C1 *func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+5)
             else
                cxB=cxI-4
                cxT=cxI+5
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,10,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= C10*ya(:,:,1)+C9*ya(:,:,2)+C8*ya(:,:,3)+C7*ya(:,:,4)+C6*ya(:,:, 5)+&
                      C5 *ya(:,:,6)+C4*ya(:,:,7)+C3*ya(:,:,8)+C2*ya(:,:,9)+C1*ya(:,:,10)
             endif
             tmp1= C10*tmp2(:,1)+C9*tmp2(:,2)+C8*tmp2(:,3)+C7*tmp2(:,4)+C6*tmp2(:, 5)+&
                   C5 *tmp2(:,6)+C4*tmp2(:,7)+C3*tmp2(:,8)+C2*tmp2(:,9)+C1*tmp2(:,10)
             funf(i,j,k)=  C10*tmp1(1)+C9*tmp1(2)+C8*tmp1(3)+C7*tmp1(4)+C6*tmp1( 5)+&
                           C5 *tmp1(6)+C4*tmp1(7)+C3*tmp1(8)+C2*tmp1(9)+C1*tmp1(10)
           endif
         endif
       endif
    enddo
   enddo
  enddo

  return

  end subroutine prolong3
#else
#error Not define Vertex nor Cell
#endif  
#endif

#endif
