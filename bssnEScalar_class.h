//$Id: bssn_class.h,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
#ifndef BSSNESCALAR_CLASS_H
#define BSSNESCALAR_CLASS_H

#ifdef newc
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
using namespace std;
#else
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#endif

#include <mpi.h>

#include "cgh.h"
#include "ShellPatch.h"
#include "misc.h"
#include "var.h"
#include "MyList.h"
#include "monitor.h"
#include "surface_integral.h"

   class bssn_class
{
  protected:
       int myrank;
       cgh *GH;
       ShellPatch *SH;
       double PhysTime;

       int checkrun; 
       char checkfilename[50];
       int Steps;
       double StartTime,TotalTime;
       double AnasTime,DumpTime,CheckTime;
       double LastAnas; 
       double Courant;
       double numepss,numepsb;
       int Symmetry;
       int maxl,decn;
       double maxrex,drex;
       int trfls,a_lev;

       double dT;
       double chitiny;

       double **Porg0,**Porgbr,**Porg,**Porg1,**Porg_rhs;
       int BH_num;

       var *phio,*trKo;
       var *gxxo,*gxyo,*gxzo,*gyyo,*gyzo,*gzzo;
       var *Axxo,*Axyo,*Axzo,*Ayyo,*Ayzo,*Azzo;
       var *Gmxo,*Gmyo,*Gmzo;
       var *Lapo,*Sfxo,*Sfyo,*Sfzo;
       var *dtSfxo,*dtSfyo,*dtSfzo;
       var *Sphio,*Spio;

       var *phi0,*trK0;
       var *gxx0,*gxy0,*gxz0,*gyy0,*gyz0,*gzz0;
       var *Axx0,*Axy0,*Axz0,*Ayy0,*Ayz0,*Azz0;
       var *Gmx0,*Gmy0,*Gmz0;
       var *Lap0,*Sfx0,*Sfy0,*Sfz0;
       var *dtSfx0,*dtSfy0,*dtSfz0;
       var *Sphi0,*Spi0;

       var *phi,*trK;
       var *gxx,*gxy,*gxz,*gyy,*gyz,*gzz;
       var *Axx,*Axy,*Axz,*Ayy,*Ayz,*Azz;
       var *Gmx,*Gmy,*Gmz;
       var *Lap,*Sfx,*Sfy,*Sfz;
       var *dtSfx,*dtSfy,*dtSfz;
       var *Sphi,*Spi;

       var *phi1,*trK1;
       var *gxx1,*gxy1,*gxz1,*gyy1,*gyz1,*gzz1;
       var *Axx1,*Axy1,*Axz1,*Ayy1,*Ayz1,*Azz1;
       var *Gmx1,*Gmy1,*Gmz1;
       var *Lap1,*Sfx1,*Sfy1,*Sfz1;
       var *dtSfx1,*dtSfy1,*dtSfz1;
       var *Sphi1,*Spi1;

       var *phi_rhs,*trK_rhs;
       var *gxx_rhs,*gxy_rhs,*gxz_rhs,*gyy_rhs,*gyz_rhs,*gzz_rhs;
       var *Axx_rhs,*Axy_rhs,*Axz_rhs,*Ayy_rhs,*Ayz_rhs,*Azz_rhs;
       var *Gmx_rhs,*Gmy_rhs,*Gmz_rhs;
       var *Lap_rhs,*Sfx_rhs,*Sfy_rhs,*Sfz_rhs;
       var *dtSfx_rhs,*dtSfy_rhs,*dtSfz_rhs;
       var *Sphi_rhs,*Spi_rhs;

       var *rho,*Sx,*Sy,*Sz,*Sxx,*Sxy,*Sxz,*Syy,*Syz,*Szz;

       var *Gamxxx,*Gamxxy,*Gamxxz,*Gamxyy,*Gamxyz,*Gamxzz;
       var *Gamyxx,*Gamyxy,*Gamyxz,*Gamyyy,*Gamyyz,*Gamyzz;
       var *Gamzxx,*Gamzxy,*Gamzxz,*Gamzyy,*Gamzyz,*Gamzzz;

       var *Rxx,*Rxy,*Rxz,*Ryy,*Ryz,*Rzz;
       
       var *Rpsi4,*Ipsi4;
       var *t1Rpsi4,*t1Ipsi4,*t2Rpsi4,*t2Ipsi4;

       MyList<var> *StateList,*SynchList_pre,*SynchList_cor,*RHSList;
       MyList<var> *OldStateList,*DumpList,*CheckList;

       monitor *ErrorMonitor,*Psi4Monitor,*BHMonitor;
       surface_integral *Waveshell;

  public:
       bssn_class(double Couranti,double StartTimei,double TotalTimei,double DumpTimei,double CheckTimei,double AnasTimei,
		  int Symmetryi,int checkruni,char *checkfilenamei,double numepssi,double numepsbi,
		  int a_levi,int maxli,int decni,double maxrexi,double drexi);
       ~bssn_class();
       void Read_Ansorg();
       void Evolve(int Steps);
       void RecursiveStep(int lev);
       void Step(int lev,int YN);
       void RestrictProlong(int lev,int YN, bool BB);
       void ProlongRestrict(int lev,int YN, bool BB);
       void Compute_Psi4(int lev);
       void Setup_Black_Hole_position();
       void compute_Porg_rhs(double **BH_PS,double **BH_RHS,var *forx,var *fory,var *forz,int lev);
};
#endif    /* BSSNESCALAR_CLASS_H */
