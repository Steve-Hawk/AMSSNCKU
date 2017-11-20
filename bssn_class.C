//$Id: bssn_class.C,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
#ifdef newc
#include <sstream>
#include <cstdio>
using namespace std;
#else
#include <stdio.h>
#endif

#include <time.h>

#include "microdef.h"
#include "misc.h"
#include "Ansorg.h"
#include "fmisc.h"
#include "Parallel.h"
#include "bssn_class.h"
#include "bssn_rhs.h"
#include "initial_puncture.h"
#include "enforce_algebra.h"
#include "rungekutta4_rout.h"
#include "sommerfeld_rout.h"
#include "getnp4.h"
#include "shellfunctions.h"

bssn_class::bssn_class(double Couranti,double StartTimei,double TotalTimei,double DumpTimei,double CheckTimei,double AnasTimei,
		  int Symmetryi,int checkruni,char *checkfilenamei,double numepssi,double numepsbi,
		  int a_levi,int maxli,int decni,double maxrexi,double drexi):
Courant(Couranti),StartTime(StartTimei),TotalTime(TotalTimei),DumpTime(DumpTimei),CheckTime(CheckTimei),AnasTime(AnasTimei),
Symmetry(Symmetryi),checkrun(checkruni),numepss(numepssi),numepsb(numepsbi),
a_lev(a_levi),maxl(maxli),decn(decni),maxrex(maxrexi),drex(drexi)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  
    if(checkrun)
    {
    }
    else
    {
      PhysTime = StartTime;
      Setup_Black_Hole_position();
    }
// setup Monitors
   {
     stringstream a_stream;
     a_stream.setf(ios::left);
     a_stream<<"# Error log information";
     ErrorMonitor=new monitor("Error.log",myrank,a_stream.str());

     a_stream.clear();
     a_stream.str("");
     a_stream<< setw(15) << "# time";
     char str[50];
     for(int pl=2;pl<maxl+1;pl++)
        for(int pm=-pl;pm<pl+1;pm++)
	{
	  sprintf(str,"R%02dm%03d",pl,pm);
	  a_stream<< setw(16) << str;
	  sprintf(str,"I%02dm%03d",pl,pm);
	  a_stream<< setw(16) << str;
	}
     Psi4Monitor=new monitor("bssn_psi4.dat",myrank,a_stream.str());

     a_stream.clear();
     a_stream.str("");
     a_stream<< setw(15) << "# time";
     BHMonitor=new monitor("bssn_BH.dat",myrank,a_stream.str());
   }
// setup sphere integration engine 
   Waveshell = new surface_integral(Symmetry); 

    trfls = 0;
    chitiny = 0;
// read parameter from file
    {
      char filename[50];
      strcpy(filename,"input.par");
      const int LEN = 256;   
      char pline[LEN];
      string str, sgrp, skey, sval;
      int sind;
      ifstream inf(filename,  ifstream::in);
      if ( !inf.good() && myrank==0 ) 
      {
        if(ErrorMonitor->outfile) ErrorMonitor->outfile<<"Can not open parameter file "<<filename<<" for inputing information of black holes"<<endl;
       	MPI_Abort(MPI_COMM_WORLD,1);
      }

      for( int i=1; inf.good(); i++)  
      {
        inf.getline(pline, LEN);
        str = pline;
     
        int status = misc::parse_parts(str, sgrp, skey, sval, sind);
        if ( status == -1 ) 
	{
	  if(ErrorMonitor->outfile) ErrorMonitor->outfile<<"error reading parameter file "<<filename<<" in line "<<i<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	}
        else if( status == 0 ) continue;
   
	if(sgrp == "BSSN" && skey == "chitiny") chitiny = atof(sval.c_str());
	else if(sgrp == "BSSN" && skey == "time refinement start from level") trfls = atoi(sval.c_str());
      }
      inf.close();
    }
// echo information of lower bound of chi 
      if(myrank==0)
      {
	cout<<"chitiny = "<<chitiny<<endl;
	cout<<"time refinement start from level #"<<trfls<<endl;
      }

    chitiny = chitiny-1; //because we have subtracted one from chi

    strcpy(checkfilename, checkfilenamei); 

    int ngfs=0;
    phio=new var("phio",ngfs++, 1, 1, 1); trKo=new var("trKo",ngfs++, 1, 1, 1);
    gxxo=new var("gxxo",ngfs++, 1, 1, 1); gxyo=new var("gxyo",ngfs++,-1,-1, 1); gxzo=new var("gxzo",ngfs++,-1, 1,-1);
    gyyo=new var("gyyo",ngfs++, 1, 1, 1); gyzo=new var("gyzo",ngfs++, 1,-1,-1); gzzo=new var("gzzo",ngfs++, 1, 1, 1);
    Axxo=new var("Axxo",ngfs++, 1, 1, 1); Axyo=new var("Axyo",ngfs++,-1,-1, 1); Axzo=new var("Axzo",ngfs++,-1, 1,-1);
    Ayyo=new var("Ayyo",ngfs++, 1, 1, 1); Ayzo=new var("Ayzo",ngfs++, 1,-1,-1); Azzo=new var("Azzo",ngfs++, 1, 1, 1);
    Gmxo=new var("Gmxo",ngfs++,-1, 1, 1); Gmyo=new var("Gmyo",ngfs++, 1,-1, 1); Gmzo=new var("Gmzo",ngfs++, 1, 1,-1);
    Lapo=new var("Lapo",ngfs++, 1, 1, 1); 
      Sfxo=new var(  "Sfxo",ngfs++,-1, 1, 1);   Sfyo=new var(  "Sfyo",ngfs++, 1,-1, 1);   Sfzo=new var(  "Sfzo",ngfs++, 1, 1,-1);
    dtSfxo=new var("dtSfxo",ngfs++,-1, 1, 1); dtSfyo=new var("dtSfyo",ngfs++, 1,-1, 1); dtSfzo=new var("dtSfzo",ngfs++, 1, 1,-1);

    phi0=new var("phi0",ngfs++, 1, 1, 1); trK0=new var("trK0",ngfs++, 1, 1, 1);
    gxx0=new var("gxx0",ngfs++, 1, 1, 1); gxy0=new var("gxy0",ngfs++,-1,-1, 1); gxz0=new var("gxz0",ngfs++,-1, 1,-1);
    gyy0=new var("gyy0",ngfs++, 1, 1, 1); gyz0=new var("gyz0",ngfs++, 1,-1,-1); gzz0=new var("gzz0",ngfs++, 1, 1, 1);
    Axx0=new var("Axx0",ngfs++, 1, 1, 1); Axy0=new var("Axy0",ngfs++,-1,-1, 1); Axz0=new var("Axz0",ngfs++,-1, 1,-1);
    Ayy0=new var("Ayy0",ngfs++, 1, 1, 1); Ayz0=new var("Ayz0",ngfs++, 1,-1,-1); Azz0=new var("Azz0",ngfs++, 1, 1, 1);
    Gmx0=new var("Gmx0",ngfs++,-1, 1, 1); Gmy0=new var("Gmy0",ngfs++, 1,-1, 1); Gmz0=new var("Gmz0",ngfs++, 1, 1,-1);
    Lap0=new var("Lap0",ngfs++, 1, 1, 1); 
      Sfx0=new var(  "Sfx0",ngfs++,-1, 1, 1);   Sfy0=new var(  "Sfy0",ngfs++, 1,-1, 1);   Sfz0=new var(  "Sfz0",ngfs++, 1, 1,-1);
    dtSfx0=new var("dtSfx0",ngfs++,-1, 1, 1); dtSfy0=new var("dtSfy0",ngfs++, 1,-1, 1); dtSfz0=new var("dtSfz0",ngfs++, 1, 1,-1);

    phi=new var("phi",ngfs++, 1, 1, 1); trK=new var("trK",ngfs++, 1, 1, 1);
    gxx=new var("gxx",ngfs++, 1, 1, 1); gxy=new var("gxy",ngfs++,-1,-1, 1); gxz=new var("gxz",ngfs++,-1, 1,-1);
    gyy=new var("gyy",ngfs++, 1, 1, 1); gyz=new var("gyz",ngfs++, 1,-1,-1); gzz=new var("gzz",ngfs++, 1, 1, 1);
    Axx=new var("Axx",ngfs++, 1, 1, 1); Axy=new var("Axy",ngfs++,-1,-1, 1); Axz=new var("Axz",ngfs++,-1, 1,-1);
    Ayy=new var("Ayy",ngfs++, 1, 1, 1); Ayz=new var("Ayz",ngfs++, 1,-1,-1); Azz=new var("Azz",ngfs++, 1, 1, 1);
    Gmx=new var("Gmx",ngfs++,-1, 1, 1); Gmy=new var("Gmy",ngfs++, 1,-1, 1); Gmz=new var("Gmz",ngfs++, 1, 1,-1);
    Lap=new var("Lap",ngfs++, 1, 1, 1); 
      Sfx=new var(  "Sfx",ngfs++,-1, 1, 1);   Sfy=new var(  "Sfy",ngfs++, 1,-1, 1);   Sfz=new var(  "Sfz",ngfs++, 1, 1,-1);
    dtSfx=new var("dtSfx",ngfs++,-1, 1, 1); dtSfy=new var("dtSfy",ngfs++, 1,-1, 1); dtSfz=new var("dtSfz",ngfs++, 1, 1,-1);

    phi1=new var("phi1",ngfs++, 1, 1, 1); trK1=new var("trK1",ngfs++, 1, 1, 1);
    gxx1=new var("gxx1",ngfs++, 1, 1, 1); gxy1=new var("gxy1",ngfs++,-1,-1, 1); gxz1=new var("gxz1",ngfs++,-1, 1,-1);
    gyy1=new var("gyy1",ngfs++, 1, 1, 1); gyz1=new var("gyz1",ngfs++, 1,-1,-1); gzz1=new var("gzz1",ngfs++, 1, 1, 1);
    Axx1=new var("Axx1",ngfs++, 1, 1, 1); Axy1=new var("Axy1",ngfs++,-1,-1, 1); Axz1=new var("Axz1",ngfs++,-1, 1,-1);
    Ayy1=new var("Ayy1",ngfs++, 1, 1, 1); Ayz1=new var("Ayz1",ngfs++, 1,-1,-1); Azz1=new var("Azz1",ngfs++, 1, 1, 1);
    Gmx1=new var("Gmx1",ngfs++,-1, 1, 1); Gmy1=new var("Gmy1",ngfs++, 1,-1, 1); Gmz1=new var("Gmz1",ngfs++, 1, 1,-1);
    Lap1=new var("Lap1",ngfs++, 1, 1, 1); 
      Sfx1=new var(  "Sfx1",ngfs++,-1, 1, 1);   Sfy1=new var(  "Sfy1",ngfs++, 1,-1, 1);   Sfz1=new var(  "Sfz1",ngfs++, 1, 1,-1);
    dtSfx1=new var("dtSfx1",ngfs++,-1, 1, 1); dtSfy1=new var("dtSfy1",ngfs++, 1,-1, 1); dtSfz1=new var("dtSfz1",ngfs++, 1, 1,-1);

    phi_rhs=new var("phi_rhs",ngfs++, 1, 1, 1); trK_rhs=new var("trK_rhs",ngfs++, 1, 1, 1);
    gxx_rhs=new var("gxx_rhs",ngfs++, 1, 1, 1); gxy_rhs=new var("gxy_rhs",ngfs++,-1,-1, 1); gxz_rhs=new var("gxz_rhs",ngfs++,-1, 1,-1);
    gyy_rhs=new var("gyy_rhs",ngfs++, 1, 1, 1); gyz_rhs=new var("gyz_rhs",ngfs++, 1,-1,-1); gzz_rhs=new var("gzz_rhs",ngfs++, 1, 1, 1);
    Axx_rhs=new var("Axx_rhs",ngfs++, 1, 1, 1); Axy_rhs=new var("Axy_rhs",ngfs++,-1,-1, 1); Axz_rhs=new var("Axz_rhs",ngfs++,-1, 1,-1);
    Ayy_rhs=new var("Ayy_rhs",ngfs++, 1, 1, 1); Ayz_rhs=new var("Ayz_rhs",ngfs++, 1,-1,-1); Azz_rhs=new var("Azz_rhs",ngfs++, 1, 1, 1);
    Gmx_rhs=new var("Gmx_rhs",ngfs++,-1, 1, 1); Gmy_rhs=new var("Gmy_rhs",ngfs++, 1,-1, 1); Gmz_rhs=new var("Gmz_rhs",ngfs++, 1, 1,-1);
    Lap_rhs=new var("Lap_rhs",ngfs++, 1, 1, 1); 
      Sfx_rhs=new var(  "Sfx_rhs",ngfs++,-1, 1, 1);   Sfy_rhs=new var(  "Sfy_rhs",ngfs++, 1,-1, 1);   Sfz_rhs=new var(  "Sfz_rhs",ngfs++, 1, 1,-1);
    dtSfx_rhs=new var("dtSfx_rhs",ngfs++,-1, 1, 1); dtSfy_rhs=new var("dtSfy_rhs",ngfs++, 1,-1, 1); dtSfz_rhs=new var("dtSfz_rhs",ngfs++, 1, 1,-1);

    rho=new var("rho",ngfs++, 1, 1, 1); Sx=new var("Sx",ngfs++,-1, 1, 1); Sy=new var("Sy",ngfs++, 1,-1, 1); Sz=new var("Sz",ngfs++, 1, 1,-1);
    Sxx=new var("Sxx",ngfs++, 1, 1, 1); Sxy=new var("Sxy",ngfs++,-1,-1, 1); Sxz=new var("Sxz",ngfs++,-1, 1,-1); 
    Syy=new var("Syy",ngfs++, 1, 1, 1); Syz=new var("Syz",ngfs++, 1,-1,-1); Szz=new var("Szz",ngfs++, 1, 1, 1);

    Gamxxx=new var("Gamxxx",ngfs++,-1, 1, 1); Gamxxy=new var("Gamxxy",ngfs++, 1,-1, 1); Gamxxz=new var("Gamxxz",ngfs++, 1, 1,-1);
    Gamxyy=new var("Gamxyy",ngfs++,-1, 1, 1); Gamxyz=new var("Gamxyz",ngfs++,-1,-1,-1); Gamxzz=new var("Gamxzz",ngfs++,-1, 1, 1);
    Gamyxx=new var("Gamyxx",ngfs++, 1,-1, 1); Gamyxy=new var("Gamyxy",ngfs++,-1, 1, 1); Gamyxz=new var("Gamyxz",ngfs++,-1,-1,-1);
    Gamyyy=new var("Gamyyy",ngfs++, 1,-1, 1); Gamyyz=new var("Gamyyz",ngfs++, 1, 1,-1); Gamyzz=new var("Gamyzz",ngfs++, 1,-1, 1);
    Gamzxx=new var("Gamzxx",ngfs++, 1, 1,-1); Gamzxy=new var("Gamzxy",ngfs++,-1,-1,-1); Gamzxz=new var("Gamzxz",ngfs++,-1, 1, 1);
    Gamzyy=new var("Gamzyy",ngfs++, 1, 1,-1); Gamzyz=new var("Gamzyz",ngfs++, 1,-1, 1); Gamzzz=new var("Gamzzz",ngfs++, 1, 1,-1);

    Rxx=new var("Rxx",ngfs++, 1, 1, 1); Rxy=new var("Rxy",ngfs++,-1,-1, 1); Rxz=new var("Rxz",ngfs++,-1, 1,-1);
    Ryy=new var("Ryy",ngfs++, 1, 1, 1); Ryz=new var("Ryz",ngfs++, 1,-1,-1); Rzz=new var("Rzz",ngfs++, 1, 1, 1);

// refer to PRD, 77, 024027 (2008)
      Rpsi4=new var(  "Rpsi4",ngfs++, 1, 1, 1);   Ipsi4=new var(  "Ipsi4",ngfs++,-1,-1,-1);
    t1Rpsi4=new var("t1Rpsi4",ngfs++, 1, 1, 1); t1Ipsi4=new var("t1Ipsi4",ngfs++,-1,-1,-1);
    t2Rpsi4=new var("t2Rpsi4",ngfs++, 1, 1, 1); t2Ipsi4=new var("t2Ipsi4",ngfs++,-1,-1,-1);

    if(myrank==0)cout<<"you have setted "<<ngfs<<" grid functions."<<endl;

   OldStateList=new MyList<var>(phio); 
   OldStateList->insert(trKo);
   OldStateList->insert(gxxo); OldStateList->insert(gxyo); OldStateList->insert(gxzo);
   OldStateList->insert(gyyo); OldStateList->insert(gyzo); OldStateList->insert(gzzo);
   OldStateList->insert(Axxo); OldStateList->insert(Axyo); OldStateList->insert(Axzo);
   OldStateList->insert(Ayyo); OldStateList->insert(Ayzo); OldStateList->insert(Azzo);
   OldStateList->insert(Gmxo); OldStateList->insert(Gmyo); OldStateList->insert(Gmzo);
   OldStateList->insert(Lapo);
   OldStateList->insert(Sfxo); OldStateList->insert(Sfyo); OldStateList->insert(Sfzo);
   OldStateList->insert(dtSfxo); OldStateList->insert(dtSfyo); OldStateList->insert(dtSfzo);
   
   StateList=new MyList<var>(phi0); 
   StateList->insert(trK0);
   StateList->insert(gxx0); StateList->insert(gxy0); StateList->insert(gxz0);
   StateList->insert(gyy0); StateList->insert(gyz0); StateList->insert(gzz0);
   StateList->insert(Axx0); StateList->insert(Axy0); StateList->insert(Axz0);
   StateList->insert(Ayy0); StateList->insert(Ayz0); StateList->insert(Azz0);
   StateList->insert(Gmx0); StateList->insert(Gmy0); StateList->insert(Gmz0);
   StateList->insert(Lap0);
   StateList->insert(Sfx0); StateList->insert(Sfy0); StateList->insert(Sfz0);
   StateList->insert(dtSfx0); StateList->insert(dtSfy0); StateList->insert(dtSfz0);

   RHSList=new MyList<var>(phi_rhs); 
   RHSList->insert(trK_rhs);
   RHSList->insert(gxx_rhs); RHSList->insert(gxy_rhs); RHSList->insert(gxz_rhs);
   RHSList->insert(gyy_rhs); RHSList->insert(gyz_rhs); RHSList->insert(gzz_rhs);
   RHSList->insert(Axx_rhs); RHSList->insert(Axy_rhs); RHSList->insert(Axz_rhs);
   RHSList->insert(Ayy_rhs); RHSList->insert(Ayz_rhs); RHSList->insert(Azz_rhs);
   RHSList->insert(Gmx_rhs); RHSList->insert(Gmy_rhs); RHSList->insert(Gmz_rhs);
   RHSList->insert(Lap_rhs);
   RHSList->insert(Sfx_rhs); RHSList->insert(Sfy_rhs); RHSList->insert(Sfz_rhs);
   RHSList->insert(dtSfx_rhs); RHSList->insert(dtSfy_rhs); RHSList->insert(dtSfz_rhs);
   
   SynchList_pre=new MyList<var>(phi); 
   SynchList_pre->insert(trK);
   SynchList_pre->insert(gxx); SynchList_pre->insert(gxy); SynchList_pre->insert(gxz);
   SynchList_pre->insert(gyy); SynchList_pre->insert(gyz); SynchList_pre->insert(gzz);
   SynchList_pre->insert(Axx); SynchList_pre->insert(Axy); SynchList_pre->insert(Axz);
   SynchList_pre->insert(Ayy); SynchList_pre->insert(Ayz); SynchList_pre->insert(Azz);
   SynchList_pre->insert(Gmx); SynchList_pre->insert(Gmy); SynchList_pre->insert(Gmz);
   SynchList_pre->insert(Lap);
   SynchList_pre->insert(Sfx); SynchList_pre->insert(Sfy); SynchList_pre->insert(Sfz);
   SynchList_pre->insert(dtSfx); SynchList_pre->insert(dtSfy); SynchList_pre->insert(dtSfz);

   SynchList_cor=new MyList<var>(phi1); 
   SynchList_cor->insert(trK1);
   SynchList_cor->insert(gxx1); SynchList_cor->insert(gxy1); SynchList_cor->insert(gxz1);
   SynchList_cor->insert(gyy1); SynchList_cor->insert(gyz1); SynchList_cor->insert(gzz1);
   SynchList_cor->insert(Axx1); SynchList_cor->insert(Axy1); SynchList_cor->insert(Axz1);
   SynchList_cor->insert(Ayy1); SynchList_cor->insert(Ayz1); SynchList_cor->insert(Azz1);
   SynchList_cor->insert(Gmx1); SynchList_cor->insert(Gmy1); SynchList_cor->insert(Gmz1);
   SynchList_cor->insert(Lap1);
   SynchList_cor->insert(Sfx1); SynchList_cor->insert(Sfy1); SynchList_cor->insert(Sfz1);
   SynchList_cor->insert(dtSfx1); SynchList_cor->insert(dtSfy1); SynchList_cor->insert(dtSfz1);
   
   DumpList=new MyList<var>(phi0);
   DumpList->insert(trK0);
   DumpList->insert(gxx0); DumpList->insert(gxy0); DumpList->insert(gxz0);
   DumpList->insert(gyy0); DumpList->insert(gyz0); DumpList->insert(gzz0);
   DumpList->insert(Axx0); DumpList->insert(Axy0); DumpList->insert(Axz0);
   DumpList->insert(Ayy0); DumpList->insert(Ayz0); DumpList->insert(Azz0);
   DumpList->insert(Gmx0); DumpList->insert(Gmy0); DumpList->insert(Gmz0);
   DumpList->insert(Lap0);
   DumpList->insert(Sfx0); DumpList->insert(Sfy0); DumpList->insert(Sfz0);
   DumpList->insert(dtSfx0); DumpList->insert(dtSfy0); DumpList->insert(dtSfz0);
   DumpList->insert(Rpsi4); DumpList->insert(Ipsi4);

   GH=new cgh(0,ngfs,Symmetry,"input.par",checkrun,ErrorMonitor);
   GH->compose_cgh(nprocs);

#ifdef WithShell   
   {
   int shapehh[dim];
   double Rrange[2];
// read parameter from file
    {
      char filename[50];
      strcpy(filename,"input.par");
      const int LEN = 256;   
      char pline[LEN];
      string str, sgrp, skey, sval;
      int sind;
      ifstream inf(filename,  ifstream::in);
      if ( !inf.good() && ErrorMonitor->outfile) 
      {
        ErrorMonitor->outfile<<"Can not open parameter file "<<filename
		             <<" for inputing information of Shell patches"<<endl;
       	MPI_Abort(MPI_COMM_WORLD,1);
      }

      for( int i=1; inf.good(); i++)  
      {
        inf.getline(pline, LEN);
        str = pline;
     
        int status = misc::parse_parts(str, sgrp, skey, sval, sind);
        if ( status == -1 ) 
	{
	  if(ErrorMonitor->outfile) ErrorMonitor->outfile<<"error reading parameter file "<<filename<<" in line "<<i<<endl; 
	  MPI_Abort(MPI_COMM_WORLD,1);
	}
        else if( status == 0 ) continue;
   
	if(sgrp == "BSSN")
	{
	 if ( skey == "Shell shape")  shapehh[sind] = atof(sval.c_str());
	 else if ( skey == "Shell R range")  Rrange[sind] = atof(sval.c_str());
	}
      }
      inf.close();
    }
   SH=new ShellPatch(0,ngfs,shapehh,Rrange,Symmetry,myrank);
   SH->matchcheck(GH->PatL[0]);
   }
   SH->compose_sh(nprocs);
   SH->setupcordtrans();
   SH->Dump_xyz(0,0,1);
   SH->setupintintstuff(nprocs,GH->PatL[0],Symmetry);
#endif   

   double h = GH->PatL[0]->data->blb->data->getdX(0);
   for(int i=1;i<dim;i++) h = min(h,GH->PatL[0]->data->blb->data->getdX(i));
   dT = Courant*h;
}
 bssn_class::~bssn_class()
{
  StateList->clearList();
  RHSList->clearList();
  OldStateList->clearList();
  SynchList_pre->clearList();
  SynchList_cor->clearList();
  DumpList->clearList();

  delete phio; delete trKo;
  delete gxxo; delete gxyo; delete gxzo; delete gyyo; delete gyzo; delete gzzo;
  delete Axxo; delete Axyo; delete Axzo; delete Ayyo; delete Ayzo; delete Azzo;
  delete Gmxo; delete Gmyo; delete Gmzo;
  delete Lapo; delete Sfxo; delete Sfyo; delete Sfzo;
  delete dtSfxo; delete dtSfyo; delete dtSfzo;

  delete phi0; delete trK0;
  delete gxx0; delete gxy0; delete gxz0; delete gyy0; delete gyz0; delete gzz0;
  delete Axx0; delete Axy0; delete Axz0; delete Ayy0; delete Ayz0; delete Azz0;
  delete Gmx0; delete Gmy0; delete Gmz0;
  delete Lap0; delete Sfx0; delete Sfy0; delete Sfz0;
  delete dtSfx0; delete dtSfy0; delete dtSfz0;

  delete phi; delete trK;
  delete gxx; delete gxy; delete gxz; delete gyy; delete gyz; delete gzz;
  delete Axx; delete Axy; delete Axz; delete Ayy; delete Ayz; delete Azz;
  delete Gmx; delete Gmy; delete Gmz;
  delete Lap; delete Sfx; delete Sfy; delete Sfz;
  delete dtSfx; delete dtSfy; delete dtSfz;

  delete phi1; delete trK1;
  delete gxx1; delete gxy1; delete gxz1; delete gyy1; delete gyz1; delete gzz1;
  delete Axx1; delete Axy1; delete Axz1; delete Ayy1; delete Ayz1; delete Azz1;
  delete Gmx1; delete Gmy1; delete Gmz1;
  delete Lap1; delete Sfx1; delete Sfy1; delete Sfz1;
  delete dtSfx1; delete dtSfy1; delete dtSfz1;

  delete phi_rhs; delete trK_rhs;
  delete gxx_rhs; delete gxy_rhs; delete gxz_rhs; delete gyy_rhs; delete gyz_rhs; delete gzz_rhs;
  delete Axx_rhs; delete Axy_rhs; delete Axz_rhs; delete Ayy_rhs; delete Ayz_rhs; delete Azz_rhs;
  delete Gmx_rhs; delete Gmy_rhs; delete Gmz_rhs;
  delete Lap_rhs; delete Sfx_rhs; delete Sfy_rhs; delete Sfz_rhs;
  delete dtSfx_rhs; delete dtSfy_rhs; delete dtSfz_rhs;

  delete rho; delete Sx; delete Sy; delete Sz; delete Sxx; delete Sxy; delete Sxz; delete Syy; delete Syz; delete Szz;

  delete Gamxxx; delete Gamxxy; delete Gamxxz; delete Gamxyy; delete Gamxyz; delete Gamxzz;
  delete Gamyxx; delete Gamyxy; delete Gamyxz; delete Gamyyy; delete Gamyyz; delete Gamyzz;
  delete Gamzxx; delete Gamzxy; delete Gamzxz; delete Gamzyy; delete Gamzyz; delete Gamzzz;

  delete Rxx; delete Rxy; delete Rxz; delete Ryy; delete Ryz; delete Rzz;
  
  delete Rpsi4; delete Ipsi4;
  delete t1Rpsi4; delete t1Ipsi4;
  delete t2Rpsi4; delete t2Ipsi4;

  delete GH;
#ifdef WithShell   
  delete SH;
#endif  

    for(int i=0;i<BH_num;i++)
    {
     delete[] Porg0[i];
     delete[] Porgbr[i];
     delete[] Porg[i];
     delete[] Porg1[i];
     delete[] Porg_rhs[i];
    }

     delete[] Porg0;
     delete[] Porgbr;
     delete[] Porg;
     delete[] Porg1;
     delete[] Porg_rhs;

   delete ErrorMonitor;
   delete Psi4Monitor;
   delete BHMonitor;
   delete Waveshell;
}
void bssn_class::Setup_Initial_Data()
{ 
  if(checkrun)
  {
  }
  else
  {
   char filename[50];
   strcpy(filename,"input.par");
   int BH_NM;
   double *Porg_here,*Pmom,*Spin,*Mass;
// read parameter from file
    {
      const int LEN = 256;   
      char pline[LEN];
      string str, sgrp, skey, sval;
      int sind;
      ifstream inf(filename,  ifstream::in);
      if ( !inf.good() && myrank==0 ) 
      {
        if(ErrorMonitor->outfile) ErrorMonitor->outfile<<"Can not open parameter file "<<filename<<" for inputing information of black holes"<<endl;
       	MPI_Abort(MPI_COMM_WORLD,1);
      }

      for( int i=1; inf.good(); i++)  
      {
        inf.getline(pline, LEN);
        str = pline;
     
        int status = misc::parse_parts(str, sgrp, skey, sval, sind);
        if ( status == -1 ) 
	{
	  if(ErrorMonitor->outfile) ErrorMonitor->outfile<<"error reading parameter file "<<filename<<" in line "<<i<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	}
        else if( status == 0 ) continue;
   
	if(sgrp == "BSSN" && skey == "BH_num") {BH_NM = atoi(sval.c_str()); break;}
      }
      inf.close();
    }
	      
    Porg_here = new double[3*BH_NM];
    Pmom = new double[3*BH_NM];
    Spin = new double[3*BH_NM];
    Mass = new double[BH_NM];
// read parameter from file
    {
      const int LEN = 256;   
      char pline[LEN];
      string str, sgrp, skey, sval;
      int sind;
      ifstream inf(filename,  ifstream::in);
      if ( !inf.good() && myrank==0 ) 
      {
        if(ErrorMonitor->outfile) ErrorMonitor->outfile<<"Can not open parameter file "<<filename
		                      <<" for inputing information of black holes"<<endl;
       	MPI_Abort(MPI_COMM_WORLD,1);
      }

      for( int i=1; inf.good(); i++)  
      {
        inf.getline(pline, LEN);
        str = pline;
     
        int status = misc::parse_parts(str, sgrp, skey, sval, sind);
        if ( status == -1 ) 
	{
	  if(ErrorMonitor->outfile) ErrorMonitor->outfile<<"error reading parameter file "<<filename<<" in line "<<i<<endl; 
	  MPI_Abort(MPI_COMM_WORLD,1);
	}
        else if( status == 0 ) continue;
   
	if(sgrp == "BSSN" && sind < BH_NM)
	{
              if ( skey == "Mass" )  Mass[sind] = atof(sval.c_str());
	 else if ( skey == "Porgx")  Porg_here[sind*3  ] = atof(sval.c_str());
	 else if ( skey == "Porgy")  Porg_here[sind*3+1] = atof(sval.c_str());
	 else if ( skey == "Porgz")  Porg_here[sind*3+2]= atof(sval.c_str());
	 else if ( skey == "Spinx")  Spin[sind*3  ] = atof(sval.c_str());
	 else if ( skey == "Spiny")  Spin[sind*3+1] = atof(sval.c_str());
	 else if ( skey == "Spinz")  Spin[sind*3+2] = atof(sval.c_str());
	 else if ( skey == "Pmomx")  Pmom[sind*3  ] = atof(sval.c_str());
	 else if ( skey == "Pmomy")  Pmom[sind*3+1] = atof(sval.c_str());
	 else if ( skey == "Pmomz")  Pmom[sind*3+2] = atof(sval.c_str());
	}
      }
      inf.close();
    }
// set initial data
    for(int lev=0;lev<GH->levels;lev++)
    {
     MyList<Patch> *Pp=GH->PatL[lev];
     while(Pp)
     {
      MyList<Block> *BL=Pp->data->blb;
      while(BL)
      {
       Block *cg=BL->data;
       if(myrank == cg->rank) 
       {
	   f_get_initial_nbhs(cg->shape,cg->X[0],cg->X[1],cg->X[2],
                     cg->fgfs[phi0->sgfn],cg->fgfs[trK0->sgfn], 
                     cg->fgfs[gxx0->sgfn],cg->fgfs[gxy0->sgfn],cg->fgfs[gxz0->sgfn],cg->fgfs[gyy0->sgfn],cg->fgfs[gyz0->sgfn],cg->fgfs[gzz0->sgfn],
                     cg->fgfs[Axx0->sgfn],cg->fgfs[Axy0->sgfn],cg->fgfs[Axz0->sgfn],cg->fgfs[Ayy0->sgfn],cg->fgfs[Ayz0->sgfn],cg->fgfs[Azz0->sgfn],
                     cg->fgfs[Gmx0->sgfn],cg->fgfs[Gmy0->sgfn],cg->fgfs[Gmz0->sgfn],
                     cg->fgfs[Lap0->sgfn],cg->fgfs[Sfx0->sgfn],cg->fgfs[Sfy0->sgfn],cg->fgfs[Sfz0->sgfn],
                     cg->fgfs[dtSfx0->sgfn],cg->fgfs[dtSfy0->sgfn],cg->fgfs[dtSfz0->sgfn],Mass,Porg_here,Pmom,Spin,BH_NM);
       }
       if(BL == Pp->data->ble) break;
       BL=BL->next;
      }
      Pp=Pp->next;
     }
    }
#ifdef WithShell   
// ShellPatch part
     MyList<ss_patch> *Pp=SH->PatL;
     while(Pp)
     {
      MyList<Block> *BL=Pp->data->blb;
      while(BL)
      {
       Block *cg=BL->data;
       if(myrank == cg->rank) 
       {
	   f_get_initial_nbhs_sh(cg->shape,cg->fgfs[Pp->data->fngfs+ShellPatch::gx],cg->fgfs[Pp->data->fngfs+ShellPatch::gy],
		     cg->fgfs[Pp->data->fngfs+ShellPatch::gz],
                     cg->fgfs[phi0->sgfn],cg->fgfs[trK0->sgfn], 
                     cg->fgfs[gxx0->sgfn],cg->fgfs[gxy0->sgfn],cg->fgfs[gxz0->sgfn],cg->fgfs[gyy0->sgfn],cg->fgfs[gyz0->sgfn],cg->fgfs[gzz0->sgfn],
                     cg->fgfs[Axx0->sgfn],cg->fgfs[Axy0->sgfn],cg->fgfs[Axz0->sgfn],cg->fgfs[Ayy0->sgfn],cg->fgfs[Ayz0->sgfn],cg->fgfs[Azz0->sgfn],
                     cg->fgfs[Gmx0->sgfn],cg->fgfs[Gmy0->sgfn],cg->fgfs[Gmz0->sgfn],
                     cg->fgfs[Lap0->sgfn],cg->fgfs[Sfx0->sgfn],cg->fgfs[Sfy0->sgfn],cg->fgfs[Sfz0->sgfn],
                     cg->fgfs[dtSfx0->sgfn],cg->fgfs[dtSfy0->sgfn],cg->fgfs[dtSfz0->sgfn],Mass,Porg_here,Pmom,Spin,BH_NM);
       }
       if(BL == Pp->data->ble) break;
       BL=BL->next;
      }
      Pp=Pp->next;
     }
#endif

   delete[] Porg_here;
   delete[] Pmom;
   delete[] Spin;
   delete[] Mass;
// dump read_in initial data 
//   SH->Synch(GH->PatL[0],StateList,Symmetry);   
//   for(int lev=0;lev<GH->levels;lev++) Parallel::Dump_Data(GH->PatL[lev],StateList,0,PhysTime,dT);
//   SH->Dump_Data(StateList,0,PhysTime,dT);
//   exit(0);
  }
}
// Read initial data solved by Ansorg, PRD 70, 064011 (2004)
void bssn_class::Read_Ansorg()
{
  if(checkrun)
  {
  }
  else
  {
   if(myrank==0)cout<<"Read initial data from Ansorg's solver, please be sure the input parameters for black holes are puncture parameters!!"<<endl;
   char filename[50];
   strcpy(filename,"input.par");
   int BH_NM;
   double *Porg_here,*Pmom,*Spin,*Mass;
// read parameter from file
    {
      const int LEN = 256;   
      char pline[LEN];
      string str, sgrp, skey, sval;
      int sind;
      ifstream inf(filename,  ifstream::in);
      if ( !inf.good() && myrank==0 ) 
      {
        if(ErrorMonitor->outfile) ErrorMonitor->outfile<<"Can not open parameter file "<<filename<<" for inputing information of black holes"<<endl;
       	MPI_Abort(MPI_COMM_WORLD,1);
      }

      for( int i=1; inf.good(); i++)  
      {
        inf.getline(pline, LEN);
        str = pline;
     
        int status = misc::parse_parts(str, sgrp, skey, sval, sind);
        if ( status == -1 ) 
	{
	  if(ErrorMonitor->outfile) ErrorMonitor->outfile<<"error reading parameter file "<<filename<<" in line "<<i<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	}
        else if( status == 0 ) continue;
   
	if(sgrp == "BSSN" && skey == "BH_num") {BH_NM = atoi(sval.c_str()); break;}
      }
      inf.close();
    }
	      
    Porg_here = new double[3*BH_NM];
    Pmom = new double[3*BH_NM];
    Spin = new double[3*BH_NM];
    Mass = new double[BH_NM];
// read parameter from file
    {
      const int LEN = 256;   
      char pline[LEN];
      string str, sgrp, skey, sval;
      int sind;
      ifstream inf(filename,  ifstream::in);
      if ( !inf.good() && myrank==0 ) 
      {
        if(ErrorMonitor->outfile) ErrorMonitor->outfile<<"Can not open parameter file "<<filename
		                      <<" for inputing information of black holes"<<endl;
       	MPI_Abort(MPI_COMM_WORLD,1);
      }

      for( int i=1; inf.good(); i++)  
      {
        inf.getline(pline, LEN);
        str = pline;
     
        int status = misc::parse_parts(str, sgrp, skey, sval, sind);
        if ( status == -1 ) 
	{
	  if(ErrorMonitor->outfile) ErrorMonitor->outfile<<"error reading parameter file "<<filename<<" in line "<<i<<endl; 
	  MPI_Abort(MPI_COMM_WORLD,1);
	}
        else if( status == 0 ) continue;
   
	if(sgrp == "BSSN" && sind < BH_NM)
	{
              if ( skey == "Mass" )  Mass[sind] = atof(sval.c_str());
	 else if ( skey == "Porgx")  Porg_here[sind*3  ] = atof(sval.c_str());
	 else if ( skey == "Porgy")  Porg_here[sind*3+1] = atof(sval.c_str());
	 else if ( skey == "Porgz")  Porg_here[sind*3+2]= atof(sval.c_str());
	 else if ( skey == "Spinx")  Spin[sind*3  ] = atof(sval.c_str());
	 else if ( skey == "Spiny")  Spin[sind*3+1] = atof(sval.c_str());
	 else if ( skey == "Spinz")  Spin[sind*3+2] = atof(sval.c_str());
	 else if ( skey == "Pmomx")  Pmom[sind*3  ] = atof(sval.c_str());
	 else if ( skey == "Pmomy")  Pmom[sind*3+1] = atof(sval.c_str());
	 else if ( skey == "Pmomz")  Pmom[sind*3+2] = atof(sval.c_str());
	}
      }
      inf.close();
    }
   int order = 6;
   Ansorg read_ansorg("Ansorg.psid",order);
// set initial data
    for(int lev=0;lev<GH->levels;lev++)
    {
     MyList<Patch> *Pp=GH->PatL[lev];
     while(Pp)
     {
      MyList<Block> *BL=Pp->data->blb;
      while(BL)
      {
       Block *cg=BL->data;
       if(myrank == cg->rank) 
       {
           for(int k=0;k<cg->shape[2];k++)
             for(int j=0;j<cg->shape[1];j++)
               for(int i=0;i<cg->shape[0];i++)
                  cg->fgfs[phi0->sgfn][i+j*cg->shape[0]+k*cg->shape[0]*cg->shape[1]] = 
			  read_ansorg.ps_u_at_xyz(cg->X[0][i],cg->X[1][j],cg->X[2][k]);

	   f_get_ansorg_nbhs(cg->shape,cg->X[0],cg->X[1],cg->X[2],
                     cg->fgfs[phi0->sgfn],cg->fgfs[trK0->sgfn], 
                     cg->fgfs[gxx0->sgfn],cg->fgfs[gxy0->sgfn],cg->fgfs[gxz0->sgfn],cg->fgfs[gyy0->sgfn],cg->fgfs[gyz0->sgfn],cg->fgfs[gzz0->sgfn],
                     cg->fgfs[Axx0->sgfn],cg->fgfs[Axy0->sgfn],cg->fgfs[Axz0->sgfn],cg->fgfs[Ayy0->sgfn],cg->fgfs[Ayz0->sgfn],cg->fgfs[Azz0->sgfn],
                     cg->fgfs[Gmx0->sgfn],cg->fgfs[Gmy0->sgfn],cg->fgfs[Gmz0->sgfn],
                     cg->fgfs[Lap0->sgfn],cg->fgfs[Sfx0->sgfn],cg->fgfs[Sfy0->sgfn],cg->fgfs[Sfz0->sgfn],
                     cg->fgfs[dtSfx0->sgfn],cg->fgfs[dtSfy0->sgfn],cg->fgfs[dtSfz0->sgfn],
                     Mass,Porg_here,Pmom,Spin,BH_NM);
       }
       if(BL == Pp->data->ble) break;
       BL=BL->next;
      }
      Pp=Pp->next;
     }
    }
#ifdef WithShell   
// ShellPatch part
     MyList<ss_patch> *Pp=SH->PatL;
     while(Pp)
     {
      MyList<Block> *BL=Pp->data->blb;
      while(BL)
      {
       Block *cg=BL->data;
       if(myrank == cg->rank) 
       {
           for(int k=0;k<cg->shape[2];k++)
             for(int j=0;j<cg->shape[1];j++)
               for(int i=0;i<cg->shape[0];i++)
                  cg->fgfs[phi0->sgfn][i+j*cg->shape[0]+k*cg->shape[0]*cg->shape[1]] = 
			  read_ansorg.ps_u_at_xyz(cg->fgfs[Pp->data->fngfs+ShellPatch::gx][i+j*cg->shape[0]+k*cg->shape[0]*cg->shape[1]],
					          cg->fgfs[Pp->data->fngfs+ShellPatch::gy][i+j*cg->shape[0]+k*cg->shape[0]*cg->shape[1]],
		                                  cg->fgfs[Pp->data->fngfs+ShellPatch::gz][i+j*cg->shape[0]+k*cg->shape[0]*cg->shape[1]]);

	   f_get_ansorg_nbhs_ss(cg->shape,cg->fgfs[Pp->data->fngfs+ShellPatch::gx],cg->fgfs[Pp->data->fngfs+ShellPatch::gy],
		     cg->fgfs[Pp->data->fngfs+ShellPatch::gz],
                     cg->fgfs[phi0->sgfn],cg->fgfs[trK0->sgfn], 
                     cg->fgfs[gxx0->sgfn],cg->fgfs[gxy0->sgfn],cg->fgfs[gxz0->sgfn],cg->fgfs[gyy0->sgfn],cg->fgfs[gyz0->sgfn],cg->fgfs[gzz0->sgfn],
                     cg->fgfs[Axx0->sgfn],cg->fgfs[Axy0->sgfn],cg->fgfs[Axz0->sgfn],cg->fgfs[Ayy0->sgfn],cg->fgfs[Ayz0->sgfn],cg->fgfs[Azz0->sgfn],
                     cg->fgfs[Gmx0->sgfn],cg->fgfs[Gmy0->sgfn],cg->fgfs[Gmz0->sgfn],
                     cg->fgfs[Lap0->sgfn],cg->fgfs[Sfx0->sgfn],cg->fgfs[Sfy0->sgfn],cg->fgfs[Sfz0->sgfn],
                     cg->fgfs[dtSfx0->sgfn],cg->fgfs[dtSfy0->sgfn],cg->fgfs[dtSfz0->sgfn],
                     Mass,Porg_here,Pmom,Spin,BH_NM);
       }
       if(BL == Pp->data->ble) break;
       BL=BL->next;
      }
      Pp=Pp->next;
     }
#endif

   delete[] Porg_here;
   delete[] Pmom;
   delete[] Spin;
   delete[] Mass;
// dump read_in initial data  
//   for(int lev=0;lev<GH->levels;lev++) Parallel::Dump_Data(GH->PatL[lev],StateList,0,PhysTime,dT);
  }
}
void bssn_class::Evolve(int Steps)
{ 
   clock_t prev_clock,curr_clock; 
   double LastDump = 0.0,LastCheck = 0.0;
   LastAnas = 0;

   double dT_mon = dT*pow(0.5,max(0,trfls));

   for(int ncount=1;ncount<Steps+1;ncount++)
   {
     if(myrank==0)curr_clock=clock();
     RecursiveStep(0);
     GH->Regrid(Symmetry,BH_num,Porgbr,Porg0,
		SynchList_cor,OldStateList,StateList,SynchList_pre,
		fgt(PhysTime-dT_mon,StartTime,dT_mon/2),ErrorMonitor);

     LastDump   += dT_mon;
     LastCheck  += dT_mon;

     if( LastDump >= DumpTime )
     {
       for(int lev=0;lev<GH->levels;lev++)
          Parallel::Dump_Data(GH->PatL[lev],DumpList,0,PhysTime,dT_mon);
#ifdef WithShell   
       SH->Dump_Data(DumpList,0,PhysTime,dT_mon);
#endif     

       LastDump = 0;
     }
     if(myrank==0)
     {
      prev_clock=curr_clock;
      curr_clock=clock();
      cout<<"Timestep # "<<ncount<<": integrating to time: "<<PhysTime<<endl;
      cout<<"used "<<(double)(curr_clock-prev_clock)/((double)CLOCKS_PER_SEC)<<" seconds!"<<endl;
     }
     if(PhysTime >= TotalTime) break;
   }
}
void bssn_class::RecursiveStep(int lev)
{
  int NoIterations = 1,YN;
  if( lev <= trfls ) NoIterations = 1;
  else               NoIterations = 2;

  for( int i = 0; i < NoIterations; i++)
  {
//     if(myrank==0) cout<<"level now = "<<lev<<" NoIteration = "<<i<<endl;
     YN = (i==NoIterations-1)?1:0; //1: same time level for coarse level and fine level
     Step(lev,YN);
     if( lev < GH->levels - 1 )
     {
       int lf = lev + 1;
       RecursiveStep(lf);
     }
     else 
       PhysTime   += dT*pow(0.5,lev);
  }
}
void bssn_class::Step(int lev,int YN)
{
    double dT_lev = dT*pow(0.5,max(lev,trfls));
    bool BB=fgt(PhysTime,StartTime,dT_lev/2);
    double ndeps=numepss;
    if(lev < GH->movls) ndeps=numepsb;
    double TRK4=PhysTime;    
    int iter_count = 0; //count RK4 substeps
    int pre=0,cor=1;
    int ERROR = 0;

  MyList<ss_patch> *sPp;
// Predictor     
  MyList<Patch> *Pp=GH->PatL[lev];
  while(Pp)
  {
    MyList<Block> *BP=Pp->data->blb;
    while(BP)
    {
      Block *cg=BP->data;
      if(myrank == cg->rank) 
      {
        f_enforce_ga(cg->shape,
                     cg->fgfs[gxx0->sgfn],cg->fgfs[gxy0->sgfn],cg->fgfs[gxz0->sgfn],cg->fgfs[gyy0->sgfn],cg->fgfs[gyz0->sgfn],cg->fgfs[gzz0->sgfn],
                     cg->fgfs[Axx0->sgfn],cg->fgfs[Axy0->sgfn],cg->fgfs[Axz0->sgfn],cg->fgfs[Ayy0->sgfn],cg->fgfs[Ayz0->sgfn],cg->fgfs[Azz0->sgfn]);

        if(f_compute_rhs_bssn(cg->shape,TRK4,cg->X[0],cg->X[1],cg->X[2], 
                     cg->fgfs[phi0->sgfn],cg->fgfs[trK0->sgfn], 
                     cg->fgfs[gxx0->sgfn],cg->fgfs[gxy0->sgfn],cg->fgfs[gxz0->sgfn],cg->fgfs[gyy0->sgfn],cg->fgfs[gyz0->sgfn],cg->fgfs[gzz0->sgfn],
                     cg->fgfs[Axx0->sgfn],cg->fgfs[Axy0->sgfn],cg->fgfs[Axz0->sgfn],cg->fgfs[Ayy0->sgfn],cg->fgfs[Ayz0->sgfn],cg->fgfs[Azz0->sgfn],
                     cg->fgfs[Gmx0->sgfn],cg->fgfs[Gmy0->sgfn],cg->fgfs[Gmz0->sgfn],
                     cg->fgfs[Lap0->sgfn],cg->fgfs[Sfx0->sgfn],cg->fgfs[Sfy0->sgfn],cg->fgfs[Sfz0->sgfn],
                     cg->fgfs[dtSfx0->sgfn],cg->fgfs[dtSfy0->sgfn],cg->fgfs[dtSfz0->sgfn],
		     cg->fgfs[phi_rhs->sgfn],cg->fgfs[trK_rhs->sgfn],
                     cg->fgfs[gxx_rhs->sgfn],cg->fgfs[gxy_rhs->sgfn],cg->fgfs[gxz_rhs->sgfn],
		     cg->fgfs[gyy_rhs->sgfn],cg->fgfs[gyz_rhs->sgfn],cg->fgfs[gzz_rhs->sgfn],
                     cg->fgfs[Axx_rhs->sgfn],cg->fgfs[Axy_rhs->sgfn],cg->fgfs[Axz_rhs->sgfn],
		     cg->fgfs[Ayy_rhs->sgfn],cg->fgfs[Ayz_rhs->sgfn],cg->fgfs[Azz_rhs->sgfn],
                     cg->fgfs[Gmx_rhs->sgfn],cg->fgfs[Gmy_rhs->sgfn],cg->fgfs[Gmz_rhs->sgfn],
                     cg->fgfs[Lap_rhs->sgfn],cg->fgfs[Sfx_rhs->sgfn],cg->fgfs[Sfy_rhs->sgfn],cg->fgfs[Sfz_rhs->sgfn],
                     cg->fgfs[dtSfx_rhs->sgfn],cg->fgfs[dtSfy_rhs->sgfn],cg->fgfs[dtSfz_rhs->sgfn],
		     cg->fgfs[rho->sgfn],cg->fgfs[Sx->sgfn],cg->fgfs[Sy->sgfn],cg->fgfs[Sz->sgfn],
		     cg->fgfs[Sxx->sgfn],cg->fgfs[Sxy->sgfn],cg->fgfs[Sxz->sgfn],cg->fgfs[Syy->sgfn],cg->fgfs[Syz->sgfn],cg->fgfs[Szz->sgfn],
                     cg->fgfs[Gamxxx->sgfn],cg->fgfs[Gamxxy->sgfn],cg->fgfs[Gamxxz->sgfn],
		     cg->fgfs[Gamxyy->sgfn],cg->fgfs[Gamxyz->sgfn],cg->fgfs[Gamxzz->sgfn],
                     cg->fgfs[Gamyxx->sgfn],cg->fgfs[Gamyxy->sgfn],cg->fgfs[Gamyxz->sgfn],
		     cg->fgfs[Gamyyy->sgfn],cg->fgfs[Gamyyz->sgfn],cg->fgfs[Gamyzz->sgfn],
                     cg->fgfs[Gamzxx->sgfn],cg->fgfs[Gamzxy->sgfn],cg->fgfs[Gamzxz->sgfn],
		     cg->fgfs[Gamzyy->sgfn],cg->fgfs[Gamzyz->sgfn],cg->fgfs[Gamzzz->sgfn],
                     cg->fgfs[Rxx->sgfn],cg->fgfs[Rxy->sgfn],cg->fgfs[Rxz->sgfn],cg->fgfs[Ryy->sgfn],cg->fgfs[Ryz->sgfn],cg->fgfs[Rzz->sgfn],
      		     Symmetry,lev,ndeps))
		     {
			cout<<"find NaN in domain: ("<<cg->bbox[0]<<":"<<cg->bbox[3]<<","<<cg->bbox[1]<<":"<<cg->bbox[4]<<","
			    <<cg->bbox[2]<<":"<<cg->bbox[5]<<")"<<endl;
			ERROR = 1;
		     }

// rk4 substep and boundary
       {
	    MyList<var> *varl0=StateList,*varl=SynchList_pre,*varlrhs=RHSList; // we do not check the correspondence here
            while(varl0)
            {
#ifndef WithShell   
		double vl=1;
                if(lev==0) //sommerfeld indeed
                   f_sommerfeld_routbam(cg->shape,cg->X[0],cg->X[1],cg->X[2],
                     Pp->data->bbox[0],Pp->data->bbox[1],Pp->data->bbox[2],Pp->data->bbox[3],Pp->data->bbox[4],Pp->data->bbox[5],
		     cg->fgfs[varlrhs->data->sgfn],
                     cg->fgfs[varl0->data->sgfn],vl,varl0->data->SoA,
	             Symmetry);

#endif
	        f_rungekutta4_rout(cg->shape, dT_lev,cg->fgfs[varl0->data->sgfn],cg->fgfs[varl->data->sgfn],cg->fgfs[varlrhs->data->sgfn],
		                   iter_count);
#ifndef WithShell  
                if(lev > 0) //fix BD point
#endif	
	        f_sommerfeld_rout(cg->shape,cg->X[0],cg->X[1],cg->X[2],
                     Pp->data->bbox[0],Pp->data->bbox[1],Pp->data->bbox[2],Pp->data->bbox[3],Pp->data->bbox[4],Pp->data->bbox[5],
		     dT_lev,cg->fgfs[phi0->sgfn],
                     cg->fgfs[Lap0->sgfn],cg->fgfs[varl0->data->sgfn],cg->fgfs[varl->data->sgfn],varl0->data->SoA,
	             Symmetry,cor);

                varl0=varl0->next;
                varl=varl->next;
                varlrhs=varlrhs->next;
            }
	}
	f_lowerboundset(cg->shape,cg->fgfs[phi->sgfn],chitiny);
      }
      if(BP==Pp->data->ble) break;
      BP=BP->next;
    }
    Pp=Pp->next;
  }
//check error information
  {int erh=ERROR;MPI_Allreduce(&erh,&ERROR,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD); }
  if(ERROR)
  {
         Parallel::Dump_Data(GH->PatL[lev],StateList,0,PhysTime,dT_lev);
         if(myrank == 0)
	 {
            if(ErrorMonitor->outfile) 
	       ErrorMonitor->outfile<<"find NaN in state variables at t = "<<PhysTime<<", lev = "<<lev<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	 }
  }

#ifdef WithShell  
// evolve Shell Patches
  if(lev==0)
  {
    sPp=SH->PatL;
    while(sPp)
    {
      MyList<Block> *BP=sPp->data->blb;
      int fngfs = sPp->data->fngfs;
      while(BP)
      {
        Block *cg=BP->data;
        if(myrank == cg->rank) 
        {
          f_enforce_ga(cg->shape,
                     cg->fgfs[gxx0->sgfn],cg->fgfs[gxy0->sgfn],cg->fgfs[gxz0->sgfn],cg->fgfs[gyy0->sgfn],cg->fgfs[gyz0->sgfn],cg->fgfs[gzz0->sgfn],
                     cg->fgfs[Axx0->sgfn],cg->fgfs[Axy0->sgfn],cg->fgfs[Axz0->sgfn],cg->fgfs[Ayy0->sgfn],cg->fgfs[Ayz0->sgfn],cg->fgfs[Azz0->sgfn]);

          if(f_compute_rhs_bssn_ss(cg->shape,TRK4,cg->X[0],cg->X[1],cg->X[2], 
                     cg->fgfs[fngfs+ShellPatch::gx],cg->fgfs[fngfs+ShellPatch::gy],cg->fgfs[fngfs+ShellPatch::gz],
                     cg->fgfs[fngfs+ShellPatch::drhodx],cg->fgfs[fngfs+ShellPatch::drhody],cg->fgfs[fngfs+ShellPatch::drhodz],
                     cg->fgfs[fngfs+ShellPatch::dsigmadx],cg->fgfs[fngfs+ShellPatch::dsigmady],cg->fgfs[fngfs+ShellPatch::dsigmadz],
                     cg->fgfs[fngfs+ShellPatch::dRdx],cg->fgfs[fngfs+ShellPatch::dRdy],cg->fgfs[fngfs+ShellPatch::dRdz],
                     cg->fgfs[fngfs+ShellPatch::drhodxx],cg->fgfs[fngfs+ShellPatch::drhodxy],cg->fgfs[fngfs+ShellPatch::drhodxz],
                     cg->fgfs[fngfs+ShellPatch::drhodyy],cg->fgfs[fngfs+ShellPatch::drhodyz],cg->fgfs[fngfs+ShellPatch::drhodzz],
                     cg->fgfs[fngfs+ShellPatch::dsigmadxx],cg->fgfs[fngfs+ShellPatch::dsigmadxy],cg->fgfs[fngfs+ShellPatch::dsigmadxz],
                     cg->fgfs[fngfs+ShellPatch::dsigmadyy],cg->fgfs[fngfs+ShellPatch::dsigmadyz],cg->fgfs[fngfs+ShellPatch::dsigmadzz],
                     cg->fgfs[fngfs+ShellPatch::dRdxx],cg->fgfs[fngfs+ShellPatch::dRdxy],cg->fgfs[fngfs+ShellPatch::dRdxz],
                     cg->fgfs[fngfs+ShellPatch::dRdyy],cg->fgfs[fngfs+ShellPatch::dRdyz],cg->fgfs[fngfs+ShellPatch::dRdzz],
                     cg->fgfs[phi0->sgfn],cg->fgfs[trK0->sgfn], 
                     cg->fgfs[gxx0->sgfn],cg->fgfs[gxy0->sgfn],cg->fgfs[gxz0->sgfn],cg->fgfs[gyy0->sgfn],cg->fgfs[gyz0->sgfn],cg->fgfs[gzz0->sgfn],
                     cg->fgfs[Axx0->sgfn],cg->fgfs[Axy0->sgfn],cg->fgfs[Axz0->sgfn],cg->fgfs[Ayy0->sgfn],cg->fgfs[Ayz0->sgfn],cg->fgfs[Azz0->sgfn],
                     cg->fgfs[Gmx0->sgfn],cg->fgfs[Gmy0->sgfn],cg->fgfs[Gmz0->sgfn],
                     cg->fgfs[Lap0->sgfn],cg->fgfs[Sfx0->sgfn],cg->fgfs[Sfy0->sgfn],cg->fgfs[Sfz0->sgfn],
                     cg->fgfs[dtSfx0->sgfn],cg->fgfs[dtSfy0->sgfn],cg->fgfs[dtSfz0->sgfn],
		     cg->fgfs[phi_rhs->sgfn],cg->fgfs[trK_rhs->sgfn],
                     cg->fgfs[gxx_rhs->sgfn],cg->fgfs[gxy_rhs->sgfn],cg->fgfs[gxz_rhs->sgfn],
		     cg->fgfs[gyy_rhs->sgfn],cg->fgfs[gyz_rhs->sgfn],cg->fgfs[gzz_rhs->sgfn],
                     cg->fgfs[Axx_rhs->sgfn],cg->fgfs[Axy_rhs->sgfn],cg->fgfs[Axz_rhs->sgfn],
		     cg->fgfs[Ayy_rhs->sgfn],cg->fgfs[Ayz_rhs->sgfn],cg->fgfs[Azz_rhs->sgfn],
                     cg->fgfs[Gmx_rhs->sgfn],cg->fgfs[Gmy_rhs->sgfn],cg->fgfs[Gmz_rhs->sgfn],
                     cg->fgfs[Lap_rhs->sgfn],cg->fgfs[Sfx_rhs->sgfn],cg->fgfs[Sfy_rhs->sgfn],cg->fgfs[Sfz_rhs->sgfn],
                     cg->fgfs[dtSfx_rhs->sgfn],cg->fgfs[dtSfy_rhs->sgfn],cg->fgfs[dtSfz_rhs->sgfn],
		     cg->fgfs[rho->sgfn],cg->fgfs[Sx->sgfn],cg->fgfs[Sy->sgfn],cg->fgfs[Sz->sgfn],
		     cg->fgfs[Sxx->sgfn],cg->fgfs[Sxy->sgfn],cg->fgfs[Sxz->sgfn],cg->fgfs[Syy->sgfn],cg->fgfs[Syz->sgfn],cg->fgfs[Szz->sgfn],
                     cg->fgfs[Gamxxx->sgfn],cg->fgfs[Gamxxy->sgfn],cg->fgfs[Gamxxz->sgfn],
		     cg->fgfs[Gamxyy->sgfn],cg->fgfs[Gamxyz->sgfn],cg->fgfs[Gamxzz->sgfn],
                     cg->fgfs[Gamyxx->sgfn],cg->fgfs[Gamyxy->sgfn],cg->fgfs[Gamyxz->sgfn],
		     cg->fgfs[Gamyyy->sgfn],cg->fgfs[Gamyyz->sgfn],cg->fgfs[Gamyzz->sgfn],
                     cg->fgfs[Gamzxx->sgfn],cg->fgfs[Gamzxy->sgfn],cg->fgfs[Gamzxz->sgfn],
		     cg->fgfs[Gamzyy->sgfn],cg->fgfs[Gamzyz->sgfn],cg->fgfs[Gamzzz->sgfn],
                     cg->fgfs[Rxx->sgfn],cg->fgfs[Rxy->sgfn],cg->fgfs[Rxz->sgfn],cg->fgfs[Ryy->sgfn],cg->fgfs[Ryz->sgfn],cg->fgfs[Rzz->sgfn],
      		     Symmetry,lev,ndeps,sPp->data->sst))
		     {
			cout<<"find NaN in Shell domain: sst = "<<sPp->data->sst<<", ("<<cg->bbox[0]<<":"<<cg->bbox[3]<<","
			<<cg->bbox[1]<<":"<<cg->bbox[4]<<","<<cg->bbox[2]<<":"<<cg->bbox[5]<<")"<<endl;
			ERROR = 1;
		     }

// rk4 substep and boundary
         {
	    MyList<var> *varl0=StateList,*varl=SynchList_pre,*varlrhs=RHSList; // we do not check the correspondence here
            while(varl0)
            {
		double vl=1;
                //sommerfeld indeed for outter boudary while fix BD for inner boundary
                f_sommerfeld_routbam_ss(cg->shape,cg->X[0],cg->X[1],cg->X[2],
                     sPp->data->bbox[0],sPp->data->bbox[1],sPp->data->bbox[2],sPp->data->bbox[3],sPp->data->bbox[4],sPp->data->bbox[5],
		     cg->fgfs[varlrhs->data->sgfn],
                     cg->fgfs[varl0->data->sgfn],vl,varl0->data->SoA,
	             Symmetry);

	        f_rungekutta4_rout(cg->shape, dT_lev,cg->fgfs[varl0->data->sgfn],cg->fgfs[varl->data->sgfn],cg->fgfs[varlrhs->data->sgfn],
		                   iter_count);

                varl0=varl0->next;
                varl=varl->next;
                varlrhs=varlrhs->next;
            }
	  }
 	  f_lowerboundset(cg->shape,cg->fgfs[phi->sgfn],chitiny);
        }
        if(BP==sPp->data->ble) break;
        BP=BP->next;
      }
      sPp=sPp->next;
    }
  }
//check error information
  {int erh=ERROR;MPI_Allreduce(&erh,&ERROR,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD); }
  if(ERROR)
  {
         SH->Dump_Data(StateList,0,PhysTime,dT_lev);
         if(myrank == 0)
	 {
            if(ErrorMonitor->outfile) 
	       ErrorMonitor->outfile<<"find NaN in state variables on Shell Patches at t = "<<PhysTime<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	 }
  }
#endif

    Parallel::Sync(GH->PatL[lev],SynchList_pre,Symmetry);

#ifdef WithShell
    if(lev==0) 
    {
      clock_t prev_clock,curr_clock;
      if(myrank==0)curr_clock=clock();
      SH->Synch(SynchList_pre,Symmetry);     
      if(myrank==0)
     {
      prev_clock=curr_clock;
      curr_clock=clock();
      cout<<"Shell stuff synchronization used "<<(double)(curr_clock-prev_clock)/((double)CLOCKS_PER_SEC)<<" seconds!"<<endl;
     }
    }
#endif

// for black hole position
  if(BH_num>0 && lev==GH->levels-1)
  {
    compute_Porg_rhs(Porg0,Porg_rhs,Sfx0,Sfy0,Sfz0,lev);
    for(int ithBH=0;ithBH<BH_num;ithBH++)
    {
       f_rungekutta4_scalar(dT_lev,Porg0[ithBH][0],Porg[ithBH][0],Porg_rhs[ithBH][0],iter_count);
       f_rungekutta4_scalar(dT_lev,Porg0[ithBH][1],Porg[ithBH][1],Porg_rhs[ithBH][1],iter_count);
       f_rungekutta4_scalar(dT_lev,Porg0[ithBH][2],Porg[ithBH][2],Porg_rhs[ithBH][2],iter_count);
       if(Symmetry > 0) Porg[ithBH][2]=fabs(Porg[ithBH][2]);
       if(Symmetry == 2)
       {
         Porg[ithBH][0]=fabs(Porg[ithBH][0]);
         Porg[ithBH][1]=fabs(Porg[ithBH][1]);
       }
       if(!finite(Porg[ithBH][0]) || !finite(Porg[ithBH][1]) || !finite(Porg[ithBH][2]))
       {
         if(ErrorMonitor->outfile) 
	    ErrorMonitor->outfile<<"predictor step finds NaN for BH's position from ("
		                 <<Porg0[ithBH][0]<<","<<Porg0[ithBH][1]<<","<<Porg0[ithBH][2]<<")"<<endl;

         MyList<var> * DG_List=new MyList<var>(Sfx0);
       	 DG_List->insert(Sfx0); DG_List->insert(Sfy0); DG_List->insert(Sfz0);
         Parallel::Dump_Data(GH->PatL[lev],DG_List,0,PhysTime,dT_lev);
	 DG_List->clearList();
       }
    }
  }
// data analysis part
// Warning NOTE: the variables1 are used as temp storege room
   if(lev==a_lev)
   {
     if(LastAnas >= AnasTime)
     {
      Compute_Psi4(lev);
      double *RP,*IP;
      int NN=0;
      for(int pl=2;pl<maxl+1;pl++)
          for(int pm=-pl;pm<pl+1;pm++) NN++;
      RP = new double[NN];
      IP = new double[NN];
      double Rex = maxrex;
      for(int i=0;i<decn;i++)
      {
	if(lev>0)                    Waveshell->surf_Wave(Rex,lev,GH, Rpsi4, Ipsi4,2,maxl,NN,RP,IP,ErrorMonitor);
	else
	{
	   if(Rex<GH->bbox[0][0][3]) Waveshell->surf_Wave(Rex,lev,GH, Rpsi4, Ipsi4,2,maxl,NN,RP,IP,ErrorMonitor);
	   else                      Waveshell->surf_Wave(Rex,lev,SH, Rpsi4, Ipsi4,2,maxl,NN,RP,IP,ErrorMonitor);
	}
	Psi4Monitor->writefile(PhysTime,NN,RP,IP);
	Rex = Rex-drex;
      }
      delete[] RP; delete[] IP;

      double *pox;
      pox=new double[dim*BH_num];
      for(int bhi=0;bhi<BH_num;bhi++)
	for(int i=0;i<dim;i++)
	   pox[dim*bhi+i] = Porg0[bhi][i];
      BHMonitor->writefile(PhysTime,dim*BH_num,pox);
      delete[] pox;
      LastAnas = 0;
     }
     LastAnas   += dT_lev;
   }
// corrector
for(iter_count = 1; iter_count < 4; iter_count++)
 {	   
// for RK4: t0, t0+dt/2, t0+dt/2, t0+dt;	 
  if(iter_count==1 || iter_count==3) TRK4 += dT_lev/2;	 
  Pp=GH->PatL[lev];
  while(Pp)
  {
    MyList<Block> *BP=Pp->data->blb;
    while(BP)
    {
      Block *cg=BP->data;
      if(myrank == cg->rank) 
      {
        f_enforce_ga(cg->shape,
                     cg->fgfs[gxx->sgfn],cg->fgfs[gxy->sgfn],cg->fgfs[gxz->sgfn],cg->fgfs[gyy->sgfn],cg->fgfs[gyz->sgfn],cg->fgfs[gzz->sgfn],
                     cg->fgfs[Axx->sgfn],cg->fgfs[Axy->sgfn],cg->fgfs[Axz->sgfn],cg->fgfs[Ayy->sgfn],cg->fgfs[Ayz->sgfn],cg->fgfs[Azz->sgfn]);

        if(f_compute_rhs_bssn(cg->shape,TRK4,cg->X[0],cg->X[1],cg->X[2], 
                     cg->fgfs[phi->sgfn],cg->fgfs[trK->sgfn], 
                     cg->fgfs[gxx->sgfn],cg->fgfs[gxy->sgfn],cg->fgfs[gxz->sgfn],cg->fgfs[gyy->sgfn],cg->fgfs[gyz->sgfn],cg->fgfs[gzz->sgfn],
                     cg->fgfs[Axx->sgfn],cg->fgfs[Axy->sgfn],cg->fgfs[Axz->sgfn],cg->fgfs[Ayy->sgfn],cg->fgfs[Ayz->sgfn],cg->fgfs[Azz->sgfn],
                     cg->fgfs[Gmx->sgfn],cg->fgfs[Gmy->sgfn],cg->fgfs[Gmz->sgfn],
                     cg->fgfs[Lap->sgfn],cg->fgfs[Sfx->sgfn],cg->fgfs[Sfy->sgfn],cg->fgfs[Sfz->sgfn],
                     cg->fgfs[dtSfx->sgfn],cg->fgfs[dtSfy->sgfn],cg->fgfs[dtSfz->sgfn],
		     cg->fgfs[phi1->sgfn],cg->fgfs[trK1->sgfn],
                     cg->fgfs[gxx1->sgfn],cg->fgfs[gxy1->sgfn],cg->fgfs[gxz1->sgfn],
		     cg->fgfs[gyy1->sgfn],cg->fgfs[gyz1->sgfn],cg->fgfs[gzz1->sgfn],
                     cg->fgfs[Axx1->sgfn],cg->fgfs[Axy1->sgfn],cg->fgfs[Axz1->sgfn],
		     cg->fgfs[Ayy1->sgfn],cg->fgfs[Ayz1->sgfn],cg->fgfs[Azz1->sgfn],
                     cg->fgfs[Gmx1->sgfn],cg->fgfs[Gmy1->sgfn],cg->fgfs[Gmz1->sgfn],
                     cg->fgfs[Lap1->sgfn],cg->fgfs[Sfx1->sgfn],cg->fgfs[Sfy1->sgfn],cg->fgfs[Sfz1->sgfn],
                     cg->fgfs[dtSfx1->sgfn],cg->fgfs[dtSfy1->sgfn],cg->fgfs[dtSfz1->sgfn],
		     cg->fgfs[rho->sgfn],cg->fgfs[Sx->sgfn],cg->fgfs[Sy->sgfn],cg->fgfs[Sz->sgfn],
		     cg->fgfs[Sxx->sgfn],cg->fgfs[Sxy->sgfn],cg->fgfs[Sxz->sgfn],cg->fgfs[Syy->sgfn],cg->fgfs[Syz->sgfn],cg->fgfs[Szz->sgfn],
                     cg->fgfs[Gamxxx->sgfn],cg->fgfs[Gamxxy->sgfn],cg->fgfs[Gamxxz->sgfn],
		     cg->fgfs[Gamxyy->sgfn],cg->fgfs[Gamxyz->sgfn],cg->fgfs[Gamxzz->sgfn],
                     cg->fgfs[Gamyxx->sgfn],cg->fgfs[Gamyxy->sgfn],cg->fgfs[Gamyxz->sgfn],
		     cg->fgfs[Gamyyy->sgfn],cg->fgfs[Gamyyz->sgfn],cg->fgfs[Gamyzz->sgfn],
                     cg->fgfs[Gamzxx->sgfn],cg->fgfs[Gamzxy->sgfn],cg->fgfs[Gamzxz->sgfn],
		     cg->fgfs[Gamzyy->sgfn],cg->fgfs[Gamzyz->sgfn],cg->fgfs[Gamzzz->sgfn],
                     cg->fgfs[Rxx->sgfn],cg->fgfs[Rxy->sgfn],cg->fgfs[Rxz->sgfn],cg->fgfs[Ryy->sgfn],cg->fgfs[Ryz->sgfn],cg->fgfs[Rzz->sgfn],
      		     Symmetry,lev,ndeps))
		     {
                        cout<<"find NaN in domain: ("<<cg->bbox[0]<<":"<<cg->bbox[3]<<","<<cg->bbox[1]<<":"<<cg->bbox[4]<<","
		            <<cg->bbox[2]<<":"<<cg->bbox[5]<<")"<<endl;
			ERROR = 1;
		     }
// rk4 substep and boundary				
	{
	    MyList<var> *varl0=StateList,*varl=SynchList_pre,*varl1=SynchList_cor,*varlrhs=RHSList; // we do not check the correspondence here
            while(varl0)
            {
#ifndef WithShell  		    
		double vl=1;
                if(lev==0) //sommerfeld indeed
                   f_sommerfeld_routbam(cg->shape,cg->X[0],cg->X[1],cg->X[2],
                     Pp->data->bbox[0],Pp->data->bbox[1],Pp->data->bbox[2],Pp->data->bbox[3],Pp->data->bbox[4],Pp->data->bbox[5],
		     cg->fgfs[varl1->data->sgfn],
                     cg->fgfs[varl->data->sgfn],vl,varl0->data->SoA,
	             Symmetry);
#endif
	        f_rungekutta4_rout(cg->shape, dT_lev,cg->fgfs[varl0->data->sgfn],cg->fgfs[varl1->data->sgfn],cg->fgfs[varlrhs->data->sgfn],
		                   iter_count);

#ifndef WithShell  
                if(lev > 0) //fix BD point
#endif	
	        f_sommerfeld_rout(cg->shape,cg->X[0],cg->X[1],cg->X[2],
                     Pp->data->bbox[0],Pp->data->bbox[1],Pp->data->bbox[2],Pp->data->bbox[3],Pp->data->bbox[4],Pp->data->bbox[5],
		     dT_lev,cg->fgfs[phi0->sgfn],
                     cg->fgfs[Lap0->sgfn],cg->fgfs[varl0->data->sgfn],cg->fgfs[varl1->data->sgfn],varl0->data->SoA,
	             Symmetry,cor);

                varl0=varl0->next;
                varl = varl->next;
                varl1=varl1->next;
                varlrhs=varlrhs->next;
            }
	}
	f_lowerboundset(cg->shape,cg->fgfs[phi1->sgfn],chitiny);
      }
      if(BP==Pp->data->ble) break;
      BP=BP->next;
    }
    Pp=Pp->next;
  }

//check error information
  {int erh=ERROR;MPI_Allreduce(&erh,&ERROR,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD); }
    if(ERROR)
    {
         Parallel::Dump_Data(GH->PatL[lev],SynchList_pre,0,PhysTime,dT_lev);
         if(myrank == 0)
	 {
            if(ErrorMonitor->outfile) 
	       ErrorMonitor->outfile<<"find NaN in RK4 substep#"<<iter_count<<" variables at t = "<<PhysTime<<", lev = "<<lev<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	 }
    }

#ifdef WithShell  
// evolve Shell Patches
  if(lev==0)
  {
    sPp=SH->PatL;
    while(sPp)
    {
      MyList<Block> *BP=sPp->data->blb;
      int fngfs = sPp->data->fngfs;
      while(BP)
      {
        Block *cg=BP->data;
        if(myrank == cg->rank) 
        {
          f_enforce_ga(cg->shape,
                     cg->fgfs[gxx->sgfn],cg->fgfs[gxy->sgfn],cg->fgfs[gxz->sgfn],cg->fgfs[gyy->sgfn],cg->fgfs[gyz->sgfn],cg->fgfs[gzz->sgfn],
                     cg->fgfs[Axx->sgfn],cg->fgfs[Axy->sgfn],cg->fgfs[Axz->sgfn],cg->fgfs[Ayy->sgfn],cg->fgfs[Ayz->sgfn],cg->fgfs[Azz->sgfn]);

          if(f_compute_rhs_bssn_ss(cg->shape,TRK4,cg->X[0],cg->X[1],cg->X[2], 
                     cg->fgfs[fngfs+ShellPatch::gx],cg->fgfs[fngfs+ShellPatch::gy],cg->fgfs[fngfs+ShellPatch::gz],
                     cg->fgfs[fngfs+ShellPatch::drhodx],cg->fgfs[fngfs+ShellPatch::drhody],cg->fgfs[fngfs+ShellPatch::drhodz],
                     cg->fgfs[fngfs+ShellPatch::dsigmadx],cg->fgfs[fngfs+ShellPatch::dsigmady],cg->fgfs[fngfs+ShellPatch::dsigmadz],
                     cg->fgfs[fngfs+ShellPatch::dRdx],cg->fgfs[fngfs+ShellPatch::dRdy],cg->fgfs[fngfs+ShellPatch::dRdz],
                     cg->fgfs[fngfs+ShellPatch::drhodxx],cg->fgfs[fngfs+ShellPatch::drhodxy],cg->fgfs[fngfs+ShellPatch::drhodxz],
                     cg->fgfs[fngfs+ShellPatch::drhodyy],cg->fgfs[fngfs+ShellPatch::drhodyz],cg->fgfs[fngfs+ShellPatch::drhodzz],
                     cg->fgfs[fngfs+ShellPatch::dsigmadxx],cg->fgfs[fngfs+ShellPatch::dsigmadxy],cg->fgfs[fngfs+ShellPatch::dsigmadxz],
                     cg->fgfs[fngfs+ShellPatch::dsigmadyy],cg->fgfs[fngfs+ShellPatch::dsigmadyz],cg->fgfs[fngfs+ShellPatch::dsigmadzz],
                     cg->fgfs[fngfs+ShellPatch::dRdxx],cg->fgfs[fngfs+ShellPatch::dRdxy],cg->fgfs[fngfs+ShellPatch::dRdxz],
                     cg->fgfs[fngfs+ShellPatch::dRdyy],cg->fgfs[fngfs+ShellPatch::dRdyz],cg->fgfs[fngfs+ShellPatch::dRdzz],
                     cg->fgfs[phi->sgfn],cg->fgfs[trK->sgfn], 
                     cg->fgfs[gxx->sgfn],cg->fgfs[gxy->sgfn],cg->fgfs[gxz->sgfn],cg->fgfs[gyy->sgfn],cg->fgfs[gyz->sgfn],cg->fgfs[gzz->sgfn],
                     cg->fgfs[Axx->sgfn],cg->fgfs[Axy->sgfn],cg->fgfs[Axz->sgfn],cg->fgfs[Ayy->sgfn],cg->fgfs[Ayz->sgfn],cg->fgfs[Azz->sgfn],
                     cg->fgfs[Gmx->sgfn],cg->fgfs[Gmy->sgfn],cg->fgfs[Gmz->sgfn],
                     cg->fgfs[Lap->sgfn],cg->fgfs[Sfx->sgfn],cg->fgfs[Sfy->sgfn],cg->fgfs[Sfz->sgfn],
                     cg->fgfs[dtSfx->sgfn],cg->fgfs[dtSfy->sgfn],cg->fgfs[dtSfz->sgfn],
		     cg->fgfs[phi1->sgfn],cg->fgfs[trK1->sgfn],
                     cg->fgfs[gxx1->sgfn],cg->fgfs[gxy1->sgfn],cg->fgfs[gxz1->sgfn],
		     cg->fgfs[gyy1->sgfn],cg->fgfs[gyz1->sgfn],cg->fgfs[gzz1->sgfn],
                     cg->fgfs[Axx1->sgfn],cg->fgfs[Axy1->sgfn],cg->fgfs[Axz1->sgfn],
		     cg->fgfs[Ayy1->sgfn],cg->fgfs[Ayz1->sgfn],cg->fgfs[Azz1->sgfn],
                     cg->fgfs[Gmx1->sgfn],cg->fgfs[Gmy1->sgfn],cg->fgfs[Gmz1->sgfn],
                     cg->fgfs[Lap1->sgfn],cg->fgfs[Sfx1->sgfn],cg->fgfs[Sfy1->sgfn],cg->fgfs[Sfz1->sgfn],
                     cg->fgfs[dtSfx1->sgfn],cg->fgfs[dtSfy1->sgfn],cg->fgfs[dtSfz1->sgfn],
		     cg->fgfs[rho->sgfn],cg->fgfs[Sx->sgfn],cg->fgfs[Sy->sgfn],cg->fgfs[Sz->sgfn],
		     cg->fgfs[Sxx->sgfn],cg->fgfs[Sxy->sgfn],cg->fgfs[Sxz->sgfn],cg->fgfs[Syy->sgfn],cg->fgfs[Syz->sgfn],cg->fgfs[Szz->sgfn],
                     cg->fgfs[Gamxxx->sgfn],cg->fgfs[Gamxxy->sgfn],cg->fgfs[Gamxxz->sgfn],
		     cg->fgfs[Gamxyy->sgfn],cg->fgfs[Gamxyz->sgfn],cg->fgfs[Gamxzz->sgfn],
                     cg->fgfs[Gamyxx->sgfn],cg->fgfs[Gamyxy->sgfn],cg->fgfs[Gamyxz->sgfn],
		     cg->fgfs[Gamyyy->sgfn],cg->fgfs[Gamyyz->sgfn],cg->fgfs[Gamyzz->sgfn],
                     cg->fgfs[Gamzxx->sgfn],cg->fgfs[Gamzxy->sgfn],cg->fgfs[Gamzxz->sgfn],
		     cg->fgfs[Gamzyy->sgfn],cg->fgfs[Gamzyz->sgfn],cg->fgfs[Gamzzz->sgfn],
                     cg->fgfs[Rxx->sgfn],cg->fgfs[Rxy->sgfn],cg->fgfs[Rxz->sgfn],cg->fgfs[Ryy->sgfn],cg->fgfs[Ryz->sgfn],cg->fgfs[Rzz->sgfn],
      		     Symmetry,lev,ndeps,sPp->data->sst))
		     {
			cout<<"find NaN in Shell domain: sst = "<<sPp->data->sst<<", ("<<cg->bbox[0]<<":"<<cg->bbox[3]<<","
			<<cg->bbox[1]<<":"<<cg->bbox[4]<<","<<cg->bbox[2]<<":"<<cg->bbox[5]<<")"<<endl;
			ERROR = 1;
		     }
// rk4 substep and boundary	
	  {
	    MyList<var> *varl0=StateList,*varl=SynchList_pre,*varl1=SynchList_cor,*varlrhs=RHSList; // we do not check the correspondence here
            while(varl0)
            {
		double vl=1;
                //sommerfeld indeed for outter boudary while fix BD for inner boundary
                f_sommerfeld_routbam_ss(cg->shape,cg->X[0],cg->X[1],cg->X[2],
                     sPp->data->bbox[0],sPp->data->bbox[1],sPp->data->bbox[2],sPp->data->bbox[3],sPp->data->bbox[4],sPp->data->bbox[5],
		     cg->fgfs[varl1->data->sgfn],
                     cg->fgfs[varl->data->sgfn],vl,varl0->data->SoA,
	             Symmetry);

	        f_rungekutta4_rout(cg->shape, dT_lev,cg->fgfs[varl0->data->sgfn],cg->fgfs[varl1->data->sgfn],cg->fgfs[varlrhs->data->sgfn],
		                   iter_count);

                varl0=varl0->next;
                varl = varl->next;
                varl1=varl1->next;
                varlrhs=varlrhs->next;
            }
	  }
	  f_lowerboundset(cg->shape,cg->fgfs[phi1->sgfn],chitiny);
        }
        if(BP==sPp->data->ble) break;
        BP=BP->next;
      }
      sPp=sPp->next;
    }
  }
//check error information
  {int erh=ERROR;MPI_Allreduce(&erh,&ERROR,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD); }
    if(ERROR)
    {
         SH->Dump_Data(SynchList_pre,0,PhysTime,dT_lev);
         if(myrank == 0)
	 {
            if(ErrorMonitor->outfile) 
	       ErrorMonitor->outfile<<"find NaN on Shell Patches in RK4 substep#"<<iter_count<<" variables at t = "<<PhysTime<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	 }
    }
#endif

    Parallel::Sync(GH->PatL[lev],SynchList_cor,Symmetry);

#ifdef WithShell   
    if(lev==0) 
    {
      clock_t prev_clock,curr_clock;
      if(myrank==0)curr_clock=clock();
      SH->Synch(SynchList_cor,Symmetry);  
      if(myrank==0)
     {
      prev_clock=curr_clock;
      curr_clock=clock();
      cout<<"Shell stuff synchronization used "<<(double)(curr_clock-prev_clock)/((double)CLOCKS_PER_SEC)<<" seconds!"<<endl;
     }
    }
#endif
// for black hole position
  if(BH_num>0 && lev==GH->levels-1)
  {
    compute_Porg_rhs(Porg,Porg1,Sfx,Sfy,Sfz,lev);
    for(int ithBH=0;ithBH<BH_num;ithBH++)
    {
       f_rungekutta4_scalar(dT_lev,Porg0[ithBH][0],Porg1[ithBH][0],Porg_rhs[ithBH][0],iter_count);
       f_rungekutta4_scalar(dT_lev,Porg0[ithBH][1],Porg1[ithBH][1],Porg_rhs[ithBH][1],iter_count);
       f_rungekutta4_scalar(dT_lev,Porg0[ithBH][2],Porg1[ithBH][2],Porg_rhs[ithBH][2],iter_count);
       if(Symmetry > 0) Porg1[ithBH][2]=fabs(Porg1[ithBH][2]);
       if(Symmetry == 2)
       {
         Porg1[ithBH][0]=fabs(Porg1[ithBH][0]);
         Porg1[ithBH][1]=fabs(Porg1[ithBH][1]);
       }
       if(!finite(Porg1[ithBH][0]) || !finite(Porg1[ithBH][1]) || !finite(Porg1[ithBH][2]))
       {
         if(ErrorMonitor->outfile) 
	    ErrorMonitor->outfile<<iter_count<<" corrector step finds NaN for BH's position from ("
		                 <<Porg[ithBH][0]<<","<<Porg[ithBH][1]<<","<<Porg[ithBH][2]<<")"<<endl;

         MyList<var> * DG_List=new MyList<var>(Sfx0);
       	 DG_List->insert(Sfx0); DG_List->insert(Sfy0); DG_List->insert(Sfz0);
         Parallel::Dump_Data(GH->PatL[lev],DG_List,0,PhysTime,dT_lev);
	 DG_List->clearList();
       }
    }
  }
//swap time level
    if( iter_count < 3 ) 
    {
     Pp=GH->PatL[lev];
     while(Pp)
     {
       MyList<Block> *BP=Pp->data->blb;
       while(BP)
       {
         Block *cg=BP->data;
         cg->swapList(SynchList_pre,SynchList_cor,myrank);
         if(BP==Pp->data->ble) break;
         BP=BP->next;
       }
       Pp=Pp->next;
     }
#ifdef WithShell    
     if(lev==0)
     {
      sPp=SH->PatL;
      while(sPp)
      {
        MyList<Block> *BP=sPp->data->blb;
        while(BP)
        {
          Block *cg=BP->data;
          cg->swapList(SynchList_pre,SynchList_cor,myrank);
          if(BP==sPp->data->ble) break;
          BP=BP->next;
        }
        sPp=sPp->next;
      }
     }
#endif   
// for black hole position
     if(BH_num>0 && lev==GH->levels-1)
     {
       for(int ithBH=0;ithBH<BH_num;ithBH++)
       {
         Porg[ithBH][0]=Porg1[ithBH][0];
         Porg[ithBH][1]=Porg1[ithBH][1];
         Porg[ithBH][2]=Porg1[ithBH][2];
       }
     }
    }
 }
// mesh refinement boundary part
   RestrictProlong(lev,YN, BB);

#ifdef WithShell       
   if(lev==0) 
   {
      clock_t prev_clock,curr_clock;
      if(myrank==0)curr_clock=clock();
      SH->CS_Inter(SynchList_cor,Symmetry);  
      if(myrank==0)
     {
      prev_clock=curr_clock;
      curr_clock=clock();
      cout<<"CS_Inter used "<<(double)(curr_clock-prev_clock)/((double)CLOCKS_PER_SEC)<<" seconds!"<<endl;
     }
   }
#endif    
// note the data structure before update
// SynchList_cor 1   -----------
//                       
// StateList     0   -----------
//
// OldStateList  old -----------
// update   
     Pp=GH->PatL[lev];
     while(Pp)
     {
       MyList<Block> *BP=Pp->data->blb;
       while(BP)
       {
         Block *cg=BP->data;
         cg->swapList(StateList,SynchList_cor,myrank);
         cg->swapList(OldStateList,SynchList_cor,myrank);
         if(BP==Pp->data->ble) break;
         BP=BP->next;
       }
       Pp=Pp->next;
     }
#ifdef WithShell  
     if(lev==0)
     {
      sPp=SH->PatL;
      while(sPp)
      {
       MyList<Block> *BP=sPp->data->blb;
       while(BP)
       {
         Block *cg=BP->data;
         cg->swapList(StateList,SynchList_cor,myrank);
         cg->swapList(OldStateList,SynchList_cor,myrank);
         if(BP==sPp->data->ble) break;
         BP=BP->next;
       }
       sPp=sPp->next;
      }
     }
#endif    
// for black hole position
     if(BH_num>0 && lev==GH->levels-1)
     {
       for(int ithBH=0;ithBH<BH_num;ithBH++)
       {
         Porg0[ithBH][0]=Porg1[ithBH][0];
         Porg0[ithBH][1]=Porg1[ithBH][1];
         Porg0[ithBH][2]=Porg1[ithBH][2];
       }
     }
}
void bssn_class::RestrictProlong(int lev,int YN, bool BB)
{
  if(lev>0)
  {
    MyList<Patch> *Pp,*Ppc;
    if(lev>trfls && YN==0) // time refinement levels and for intermediat time level
    {
       Pp=GH->PatL[lev-1];
       while(Pp)
       {
         if(BB)  Parallel::prepare_inter_time_level(Pp->data,StateList,OldStateList,SynchList_cor,
		                                    SynchList_pre, 0);  // use SynchList_pre as temporal storage space
         else    Parallel::prepare_inter_time_level(Pp->data,StateList,OldStateList,
		                                    SynchList_pre, 0);  // use SynchList_pre as temporal storage space
	 Pp=Pp->next;
       }

       Parallel::Restrict(GH->PatL[lev-1],GH->PatL[lev],SynchList_cor,SynchList_pre,Symmetry);
      
       Parallel::Sync(GH->PatL[lev-1],SynchList_pre,Symmetry);

       Ppc=GH->PatL[lev-1];
       while(Ppc)
       {
         Pp=GH->PatL[lev];
	 while(Pp)
	 {
           Parallel::OutBdLow2Hi(Ppc->data,Pp->data,SynchList_pre,SynchList_cor,Symmetry);
	   Pp=Pp->next;
	 }
	 Ppc=Ppc->next;
       }
    }   
    else           // no time refinement levels and for all same time levels
    {
       Parallel::Restrict(GH->PatL[lev-1],GH->PatL[lev],SynchList_cor,StateList,Symmetry);
	
       Parallel::Sync(GH->PatL[lev-1],StateList,Symmetry);

       Ppc=GH->PatL[lev-1];
       while(Ppc)
       {
         Pp=GH->PatL[lev];
	 while(Pp)
	 {
           Parallel::OutBdLow2Hi(Ppc->data,Pp->data,StateList,SynchList_cor,Symmetry);
	   Pp=Pp->next;
	 }
	 Ppc=Ppc->next;
       }
    }

    Parallel::Sync(GH->PatL[lev],SynchList_cor,Symmetry);
  }
}
void bssn_class::ProlongRestrict(int lev,int YN, bool BB)
{
  if(lev>0)
  {
    MyList<Patch> *Pp,*Ppc;
    if(lev>trfls && YN==0) // time refinement levels and for intermediat time level
    {
       Pp=GH->PatL[lev-1];
       while(Pp)
       {
         if(BB)  Parallel::prepare_inter_time_level(Pp->data,StateList,OldStateList,SynchList_cor,
		                                    SynchList_pre, 0);  // use SynchList_pre as temporal storage space
         else    Parallel::prepare_inter_time_level(Pp->data,StateList,OldStateList,
		                                    SynchList_pre, 0);  // use SynchList_pre as temporal storage space
	 Pp=Pp->next;
       }

       Ppc=GH->PatL[lev-1];
       while(Ppc)
       {
         Pp=GH->PatL[lev];
	 while(Pp)
	 {
           Parallel::OutBdLow2Hi(Ppc->data,Pp->data,SynchList_pre,SynchList_cor,Symmetry);
	   Pp=Pp->next;
	 }
	 Ppc=Ppc->next;
       }
    }   
    else           // no time refinement levels and for all same time levels
    {
       Ppc=GH->PatL[lev-1];
       while(Ppc)
       {
         Pp=GH->PatL[lev];
	 while(Pp)
	 {
           Parallel::OutBdLow2Hi(Ppc->data,Pp->data,StateList,SynchList_cor,Symmetry);
	   Pp=Pp->next;
	 }
	 Ppc=Ppc->next;
       }
            
       Parallel::Restrict(GH->PatL[lev-1],GH->PatL[lev],SynchList_cor,StateList,Symmetry);

       Parallel::Sync(GH->PatL[lev-1],StateList,Symmetry);
    }

    Parallel::Sync(GH->PatL[lev],SynchList_cor,Symmetry);
  }
}
void bssn_class::Compute_Psi4(int lev)
{
  MyList<Patch> *Pp=GH->PatL[lev];
  while(Pp)
  {
    MyList<Block> *BP=Pp->data->blb;
    while(BP)
    {
      Block *cg=BP->data;
      if(myrank == cg->rank) 
      {      
// the input arguments Gamma^i_jk and R_ij do not need synch, because we do not need to derivate them	       
         f_getnp4(cg->shape,cg->X[0],cg->X[1],cg->X[2],
                  cg->fgfs[phi0->sgfn],cg->fgfs[trK0->sgfn], 
                  cg->fgfs[gxx0->sgfn],cg->fgfs[gxy0->sgfn],cg->fgfs[gxz0->sgfn],cg->fgfs[gyy0->sgfn],cg->fgfs[gyz0->sgfn],cg->fgfs[gzz0->sgfn],
                  cg->fgfs[Axx0->sgfn],cg->fgfs[Axy0->sgfn],cg->fgfs[Axz0->sgfn],cg->fgfs[Ayy0->sgfn],cg->fgfs[Ayz0->sgfn],cg->fgfs[Azz0->sgfn],
                  cg->fgfs[Gamxxx->sgfn],cg->fgfs[Gamxxy->sgfn],cg->fgfs[Gamxxz->sgfn],
                  cg->fgfs[Gamxyy->sgfn],cg->fgfs[Gamxyz->sgfn],cg->fgfs[Gamxzz->sgfn],
                  cg->fgfs[Gamyxx->sgfn],cg->fgfs[Gamyxy->sgfn],cg->fgfs[Gamyxz->sgfn],
	          cg->fgfs[Gamyyy->sgfn],cg->fgfs[Gamyyz->sgfn],cg->fgfs[Gamyzz->sgfn],
                  cg->fgfs[Gamzxx->sgfn],cg->fgfs[Gamzxy->sgfn],cg->fgfs[Gamzxz->sgfn],
	          cg->fgfs[Gamzyy->sgfn],cg->fgfs[Gamzyz->sgfn],cg->fgfs[Gamzzz->sgfn],
                  cg->fgfs[Rxx->sgfn],cg->fgfs[Rxy->sgfn],cg->fgfs[Rxz->sgfn],cg->fgfs[Ryy->sgfn],cg->fgfs[Ryz->sgfn],cg->fgfs[Rzz->sgfn],
	          cg->fgfs[Rpsi4->sgfn],cg->fgfs[Ipsi4->sgfn],
	          Symmetry);	  
      }
      if(BP==Pp->data->ble) break;
      BP=BP->next;
    }
    Pp=Pp->next;
  }

#ifdef WithShell
// ShellPatch part    
    if(lev==0) 
    {
     MyList<ss_patch> *Pp=SH->PatL;
     while(Pp)
     {
      MyList<Block> *BL=Pp->data->blb;
      int fngfs = Pp->data->fngfs;
      while(BL)
      {
       Block *cg=BL->data;
       if(myrank == cg->rank) 
       {
         f_getnp4_ss(cg->shape,cg->X[0],cg->X[1],cg->X[2], 
                  cg->fgfs[fngfs+ShellPatch::gx],cg->fgfs[fngfs+ShellPatch::gy],cg->fgfs[fngfs+ShellPatch::gz],
                  cg->fgfs[fngfs+ShellPatch::drhodx],cg->fgfs[fngfs+ShellPatch::drhody],cg->fgfs[fngfs+ShellPatch::drhodz],
                  cg->fgfs[fngfs+ShellPatch::dsigmadx],cg->fgfs[fngfs+ShellPatch::dsigmady],cg->fgfs[fngfs+ShellPatch::dsigmadz],
                  cg->fgfs[fngfs+ShellPatch::dRdx],cg->fgfs[fngfs+ShellPatch::dRdy],cg->fgfs[fngfs+ShellPatch::dRdz],
                  cg->fgfs[fngfs+ShellPatch::drhodxx],cg->fgfs[fngfs+ShellPatch::drhodxy],cg->fgfs[fngfs+ShellPatch::drhodxz],
                  cg->fgfs[fngfs+ShellPatch::drhodyy],cg->fgfs[fngfs+ShellPatch::drhodyz],cg->fgfs[fngfs+ShellPatch::drhodzz],
                  cg->fgfs[fngfs+ShellPatch::dsigmadxx],cg->fgfs[fngfs+ShellPatch::dsigmadxy],cg->fgfs[fngfs+ShellPatch::dsigmadxz],
                  cg->fgfs[fngfs+ShellPatch::dsigmadyy],cg->fgfs[fngfs+ShellPatch::dsigmadyz],cg->fgfs[fngfs+ShellPatch::dsigmadzz],
                  cg->fgfs[fngfs+ShellPatch::dRdxx],cg->fgfs[fngfs+ShellPatch::dRdxy],cg->fgfs[fngfs+ShellPatch::dRdxz],
                  cg->fgfs[fngfs+ShellPatch::dRdyy],cg->fgfs[fngfs+ShellPatch::dRdyz],cg->fgfs[fngfs+ShellPatch::dRdzz],
                  cg->fgfs[phi0->sgfn],cg->fgfs[trK0->sgfn], 
                  cg->fgfs[gxx0->sgfn],cg->fgfs[gxy0->sgfn],cg->fgfs[gxz0->sgfn],cg->fgfs[gyy0->sgfn],cg->fgfs[gyz0->sgfn],cg->fgfs[gzz0->sgfn],
                  cg->fgfs[Axx0->sgfn],cg->fgfs[Axy0->sgfn],cg->fgfs[Axz0->sgfn],cg->fgfs[Ayy0->sgfn],cg->fgfs[Ayz0->sgfn],cg->fgfs[Azz0->sgfn],
                  cg->fgfs[Gamxxx->sgfn],cg->fgfs[Gamxxy->sgfn],cg->fgfs[Gamxxz->sgfn],
                  cg->fgfs[Gamxyy->sgfn],cg->fgfs[Gamxyz->sgfn],cg->fgfs[Gamxzz->sgfn],
                  cg->fgfs[Gamyxx->sgfn],cg->fgfs[Gamyxy->sgfn],cg->fgfs[Gamyxz->sgfn],
	          cg->fgfs[Gamyyy->sgfn],cg->fgfs[Gamyyz->sgfn],cg->fgfs[Gamyzz->sgfn],
                  cg->fgfs[Gamzxx->sgfn],cg->fgfs[Gamzxy->sgfn],cg->fgfs[Gamzxz->sgfn],
	          cg->fgfs[Gamzyy->sgfn],cg->fgfs[Gamzyz->sgfn],cg->fgfs[Gamzzz->sgfn],
                  cg->fgfs[Rxx->sgfn],cg->fgfs[Rxy->sgfn],cg->fgfs[Rxz->sgfn],cg->fgfs[Ryy->sgfn],cg->fgfs[Ryz->sgfn],cg->fgfs[Rzz->sgfn],
	          cg->fgfs[Rpsi4->sgfn],cg->fgfs[Ipsi4->sgfn],
	          Symmetry,Pp->data->sst);
       }
       if(BL == Pp->data->ble) break;
       BL=BL->next;
      }
      Pp=Pp->next;
     }
    }
#endif    

    MyList<var> * DG_List=new MyList<var>(Rpsi4);
    DG_List->insert(Ipsi4);
    Parallel::Sync(GH->PatL[lev],DG_List,Symmetry);
#ifdef WithShell   
    if(lev==0) 
    {
      SH->Synch(DG_List,Symmetry);  
    }
#endif
    DG_List->clearList();
}
void bssn_class::Setup_Black_Hole_position()
{ 
   char filename[50];
   strcpy(filename,"input.par");
   double *Porg_here,*Pmom,*Spin,*Mass;
// read parameter from file
    {
      const int LEN = 256;   
      char pline[LEN];
      string str, sgrp, skey, sval;
      int sind;
      ifstream inf(filename,  ifstream::in);
      if ( !inf.good() && myrank==0 ) 
      {
        if(ErrorMonitor->outfile) ErrorMonitor->outfile<<"Can not open parameter file "<<filename<<" for inputing information of black holes"<<endl;
       	MPI_Abort(MPI_COMM_WORLD,1);
      }

      for( int i=1; inf.good(); i++)  
      {
        inf.getline(pline, LEN);
        str = pline;
     
        int status = misc::parse_parts(str, sgrp, skey, sval, sind);
        if ( status == -1 ) 
	{
	  if(ErrorMonitor->outfile) ErrorMonitor->outfile<<"error reading parameter file "<<filename<<" in line "<<i<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	}
        else if( status == 0 ) continue;
   
	if(sgrp == "BSSN" && skey == "BH_num") {BH_num = atoi(sval.c_str()); break;}
      }
      inf.close();
    }
	      
    Porg_here = new double[3*BH_num];
    Pmom = new double[3*BH_num];
    Spin = new double[3*BH_num];
    Mass = new double[BH_num];
// read parameter from file
    {
      const int LEN = 256;   
      char pline[LEN];
      string str, sgrp, skey, sval;
      int sind;
      ifstream inf(filename,  ifstream::in);
      if ( !inf.good() && myrank==0 ) 
      {
        if(ErrorMonitor->outfile) ErrorMonitor->outfile<<"Can not open parameter file "<<filename
		                      <<" for inputing information of black holes"<<endl;
       	MPI_Abort(MPI_COMM_WORLD,1);
      }

      for( int i=1; inf.good(); i++)  
      {
        inf.getline(pline, LEN);
        str = pline;
     
        int status = misc::parse_parts(str, sgrp, skey, sval, sind);
        if ( status == -1 ) 
	{
	  if(ErrorMonitor->outfile) ErrorMonitor->outfile<<"error reading parameter file "<<filename<<" in line "<<i<<endl; 
	  MPI_Abort(MPI_COMM_WORLD,1);
	}
        else if( status == 0 ) continue;
   
	if(sgrp == "BSSN" && sind < BH_num)
	{
              if ( skey == "Mass" )  Mass[sind] = atof(sval.c_str());
	 else if ( skey == "Porgx")  Porg_here[sind*3  ] = atof(sval.c_str());
	 else if ( skey == "Porgy")  Porg_here[sind*3+1] = atof(sval.c_str());
	 else if ( skey == "Porgz")  Porg_here[sind*3+2]= atof(sval.c_str());
	 else if ( skey == "Spinx")  Spin[sind*3  ] = atof(sval.c_str());
	 else if ( skey == "Spiny")  Spin[sind*3+1] = atof(sval.c_str());
	 else if ( skey == "Spinz")  Spin[sind*3+2] = atof(sval.c_str());
	 else if ( skey == "Pmomx")  Pmom[sind*3  ] = atof(sval.c_str());
	 else if ( skey == "Pmomy")  Pmom[sind*3+1] = atof(sval.c_str());
	 else if ( skey == "Pmomz")  Pmom[sind*3+2] = atof(sval.c_str());
	}
      }
      inf.close();
    }
// echo information of Black holes
      if(myrank==0)
      {
	cout<<"initial information of "<<BH_num<<" Black Hole(s)"<<endl;
        cout<<setw(16) << "Mass"
            << setw(16) << "x"
	    << setw(16) << "y"
	    << setw(16) << "z"
	    << setw(16) << "Px"
	    << setw(16) << "Py"
	    << setw(16) << "Pz"
	    << setw(16) << "Sx"
	    << setw(16) << "Sy"
	    << setw(16) << "Sz"<<endl;
	for(int i=0;i<BH_num;i++)
	{
        cout<<setw(16) << Mass[i]
            << setw(16) << Porg_here[i*3  ]
	    << setw(16) << Porg_here[i*3+1]
	    << setw(16) << Porg_here[i*3+2]
	    << setw(16) << Pmom[i*3  ]
	    << setw(16) << Pmom[i*3+1]
	    << setw(16) << Pmom[i*3+2]
	    << setw(16) << Spin[i*3  ]
	    << setw(16) << Spin[i*3+1]
	    << setw(16) << Spin[i*3+2]<<endl;
	}
      }
// set initial data
    Porg0 = new double*[BH_num];
    Porgbr = new double*[BH_num];
    Porg  = new double*[BH_num];
    Porg1 = new double*[BH_num];
    Porg_rhs=new double *[BH_num];
    for(int i=0;i<BH_num;i++)
    {
     Porg0[i] = new double[3];
     Porgbr[i] = new double[3];
     for(int j=0;j<dim;j++) Porgbr[i][j]=Porg0[i][j]=Porg_here[i*dim+j];
     Porg[i]  = new double[3];
     Porg1[i] = new double[3];
     Porg_rhs[i]=new double [3];
    }

   delete[] Porg_here;
   delete[] Pmom;
   delete[] Spin;
   delete[] Mass;
}
void bssn_class::compute_Porg_rhs(double **BH_PS,double **BH_RHS,var *forx,var *fory,var *forz,int lev)
{
  const int InList = 3;

  MyList<var> * DG_List=new MyList<var>(forx);
  DG_List->insert(fory); DG_List->insert(forz);

  int n;
  double *x1,*y1,*z1;
  double *shellf;
  shellf=new double[3*BH_num];
  double *pox[3];
  for(int i=0;i<3;i++) pox[i] = new double[BH_num];
  for( n = 0; n < BH_num; n++)
  {
    pox[0][n] = BH_PS[n][0];
    pox[1][n] = BH_PS[n][1];
    pox[2][n] = BH_PS[n][2];
  }

  if(!Parallel::PatList_Interp_Points(GH->PatL[lev],DG_List,BH_num,pox,shellf,Symmetry))
	  ErrorMonitor->outfile<<"fail to find black holes at t = "<<PhysTime<<endl; 
  
  for( n = 0; n < BH_num; n++)
  {
    BH_RHS[n][0]=-shellf[3*n  ];
    BH_RHS[n][1]=-shellf[3*n+1];
    BH_RHS[n][2]=-shellf[3*n+2];
  }
   
  DG_List->clearList();
  delete[] shellf;
  for(int i=0;i<3;i++) delete[] pox[i];
}
