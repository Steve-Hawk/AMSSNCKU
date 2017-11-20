//$Id: surface_integral.C,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
//----------------------------------------------------------------
// Using Gauss-Legendre quadrature in theta direction
// and   trapezoidal rule in phi direction (from Second Euler-Maclaurin summation formula, we can see that
// this method gives expolential convergence for periodic function)
//----------------------------------------------------------------
#ifdef newc
#include <iostream>
#include <iomanip>
#include <fstream>
#include <strstream>
#include <cmath>
using namespace std;
#else
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <string.h>
#include <math.h>
#endif
#include <mpi.h>

#include "misc.h"
#include "cgh.h"
#include "Parallel.h"
#include "surface_integral.h"

#define PI M_PI
//|============================================================================
//| Constructor
//|============================================================================

  surface_integral::surface_integral(int iSymmetry):
   Symmetry(iSymmetry)
  {
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    MPI_Comm_size(MPI_COMM_WORLD,&cpusize);
    int N=40;
// read parameter from file
    {
      const int LEN = 256;   
      char pline[LEN];
      string str, sgrp, skey, sval;
      int sind;
      const char pname[] = "input.par"; 
      ifstream inf(pname,  ifstream::in);
      if ( !inf.good() && myrank==0 ) { cout<<"Can not open parameter file "<<pname<<endl; MPI_Abort(MPI_COMM_WORLD,1);}

      for( int i=1; inf.good(); i++)  
      {
        inf.getline(pline, LEN);
        str = pline;
     
        int status = misc::parse_parts(str, sgrp, skey, sval, sind);
        if ( status == -1 ) { cout<<"error reading parameter file "<<pname<<" in line "<<i<<endl; MPI_Abort(MPI_COMM_WORLD,1);} 
        else if( status == 0 ) continue;
   
	if(sgrp == "SurfaceIntegral")
	{
              if ( skey == "number of points for quarter sphere" )  N = atoi(sval.c_str()); 
	}
      }
      inf.close();
    }
//|-----number of points for whole [0,pi] x [0,2pi]
    N_phi  = 4*N;   // for simplicity, we require this number must be 4*N
    N_theta= 2*N;   //                                                2*N
    
    if(myrank==0)
    {
	cout<<"-----------------------------------------------------------------------"<<endl;
#ifdef GaussInt
	cout<<"    spherical integration for wave form extraction with Gauss method   "<<endl;
#else
	cout<<" spherical integration for wave form extraction with mid point method  "<<endl;
#endif     
     	cout<<"N_phi = "<<N_phi<<endl;
        cout<<"N_theta = "<<N_theta<<endl;
	cout<<"-----------------------------------------------------------------------"<<endl;
    }
   
#ifdef GaussInt
//  weight function cover all of [0,pi]
   arcostheta = new double[N_theta];
   wtcostheta = new double[N_theta];

// note: theta in [0,pi/2], upper half sphere, corresponds to 1 < costheta < 0
   misc::gaulegf(-1.0, 1.0,arcostheta,wtcostheta,N_theta);
// due to symmetry, I need first half array corresponds to upper sphere, note these two arrays must match each other
   misc::inversearray(arcostheta,N_theta);
   misc::inversearray(wtcostheta,N_theta);
#endif   

  if(      Symmetry == 2 )
   {
         N_phi=N_phi/4;
         N_theta=N_theta/2;
         dphi =  PI / (2.0 * N_phi);
    dcostheta = 1.0 / N_theta;
   }
  else if( Symmetry == 1 )
   {
         N_theta=N_theta/2;
         dphi = 2.0 * PI / N_phi;
    dcostheta = 1.0 / N_theta;
   }
  else if( Symmetry == 0 )
   {
         dphi = 2.0 * PI / N_phi;
    dcostheta = 2.0 / N_theta;
   }
  else if(myrank==0)
  {
       cout<<"surface_integral::surface_integral: not supported Symmetry setting!"<<endl;
       MPI_Abort(MPI_COMM_WORLD,1);
  }

#ifndef GaussInt
//  weight function cover all of [0,pi]
   arcostheta = new double[N_theta];
#endif   
  n_tot = N_theta * N_phi;
  nx_g = new double[n_tot]; ny_g = new double[n_tot]; nz_g = new double[n_tot];

  int n = 0;
  double costheta, sintheta, ph;

  for(int i = 0; i < N_theta; ++i)
   {
#ifndef GaussInt	   
    arcostheta[i] = 1.0 - (i + 0.5) * dcostheta;
#endif    
    costheta = arcostheta[i];
    sintheta = sqrt(1.0 - costheta * costheta);

    for(int j = 0; j < N_phi; ++j)
     {
      ph = (j + 0.5) * dphi;
//normal vector respect to the constant R sphere      
      nx_g[n] = sintheta * cos(ph);
      ny_g[n] = sintheta * sin(ph);
      nz_g[n] = costheta;
      n++;
     }
   }
  }

//|============================================================================
//| Destructor
//|============================================================================
 surface_integral::~surface_integral()
{
  delete[] nx_g; delete[] ny_g; delete[] nz_g;
  delete[] arcostheta;
#ifdef GaussInt  
  delete[] wtcostheta;
#endif  
}
//|----------------------------------------------------------------
//  spin weighted spinw component of psi4, general routine
//  l takes from spinw to maxl; m takes from -l to l
//|----------------------------------------------------------------
  void surface_integral::surf_Wave(double rex,int lev,cgh *GH, var *Rpsi4, var *Ipsi4,
		                   int spinw,int maxl, int NN, double *RP, double *IP,
				   monitor *Monitor)  //NN is the length of RP and IP
   {
    if(myrank==0 && GH->grids[lev] != 1 && Monitor->outfile) Monitor->outfile<<"WARNING: surface integral on multipatches"<<endl;

    const int InList = 2;

    MyList<var> * DG_List=new MyList<var>(Rpsi4);
    DG_List->insert(Ipsi4);

    int n;
    double *pox[3];
    for(int i=0;i<3;i++) pox[i] = new double[n_tot];
    for( n = 0; n < n_tot; n++)
       {
        pox[0][n] = rex * nx_g[n];
        pox[1][n] = rex * ny_g[n];
        pox[2][n] = rex * nz_g[n];
       }

    double *shellf;
    shellf=new double[n_tot*InList];

    GH->PatL[lev]->data->Interp_Points(DG_List,n_tot,pox,shellf,Symmetry);    

    int mp, Lp, Nmin, Nmax;

    mp = n_tot / cpusize;
    Lp = n_tot - cpusize * mp;

    if( Lp > myrank ) { Nmin = myrank * mp + myrank; Nmax = Nmin + mp    ;}
    else              { Nmin = myrank * mp + Lp;     Nmax = Nmin + mp - 1;}

//|~~~~~> Integrate the dot product of Dphi with the surface normal.

    double *RP_out, *IP_out;
    RP_out = new double[NN];
    IP_out = new double[NN];    
    
    for(int ii=0;ii<NN;ii++)
    {
      RP_out[ii] = 0;
      IP_out[ii] = 0;
    }
// theta part    
    double costheta, thetap;
    double cosmphi,sinmphi;

    int i,j;
    int lpsy=0;
         if( Symmetry == 0 )     lpsy=1;
    else if( Symmetry == 1 )     lpsy=2;
    else if( Symmetry == 2 )     lpsy=8;

    double psi4RR,psi4II;
    for( n = Nmin; n <= Nmax; n++) 
     {
//       need round off always	     
        i = int(n/N_phi); // int(1.723) = 1, int(-1.732) = -1
        j = n - i * N_phi;
        
	int countlm=0;
	for(int pl=spinw;pl<maxl+1;pl++)
          for(int pm=-pl;pm<pl+1;pm++)
	  {
	for(int lp=0;lp<lpsy;lp++)
	{
 	 switch(lp)
	 {
	  case 0:  //+++ (theta, phi)
          costheta = arcostheta[i];
	  cosmphi = cos(pm * (j+0.5) * dphi);
	  sinmphi = sin(pm * (j+0.5) * dphi);
	  psi4RR = shellf[InList*n  ];
	  psi4II = shellf[InList*n+1];
	  break;
	  case 1:  //++- (pi-theta, phi)
          costheta = -arcostheta[i];
	  cosmphi = cos(pm * (j+0.5) * dphi);
	  sinmphi = sin(pm * (j+0.5) * dphi);
	  psi4RR = shellf[InList*n  ];
	  psi4II =-shellf[InList*n+1];
	  break;
	  case 2:  //+-+ (theta, 2*pi-phi)
          costheta = arcostheta[i];
	  cosmphi = cos(pm * (j+0.5) * dphi);
	  sinmphi =-sin(pm * (j+0.5) * dphi);
	  psi4RR = shellf[InList*n  ];
	  psi4II =-shellf[InList*n+1];
	  break;
	  case 3:  //+-- (pi-theta, 2*pi-phi)
          costheta = -arcostheta[i];
	  cosmphi = cos(pm * (j+0.5) * dphi);
	  sinmphi =-sin(pm * (j+0.5) * dphi);
	  psi4RR = shellf[InList*n  ];
	  psi4II = shellf[InList*n+1];
	  break;
	  case 4:  //-++ (theta, pi-phi)
          costheta = arcostheta[i];
	  cosmphi = cos(pm * (PI - (j+0.5) * dphi));
	  sinmphi = sin(pm * (PI - (j+0.5) * dphi));
	  psi4RR = shellf[InList*n  ];
	  psi4II =-shellf[InList*n+1];
	  break;
	  case 5:  //-+- (pi-theta, pi-phi)
          costheta = -arcostheta[i];
 	  cosmphi = cos(pm * (PI - (j+0.5) * dphi));
	  sinmphi = sin(pm * (PI - (j+0.5) * dphi));
	  psi4RR = shellf[InList*n  ];
	  psi4II = shellf[InList*n+1];
	  break;
	  case 6:  //--+ (theta, pi+phi)
          costheta = arcostheta[i];
	  cosmphi = cos(pm * (PI + (j+0.5) * dphi));
	  sinmphi = sin(pm * (PI + (j+0.5) * dphi));
	  psi4RR = shellf[InList*n  ];
	  psi4II = shellf[InList*n+1];
	  break;
	  case 7:  //--- (pi-theta, pi+phi)
          costheta = -arcostheta[i];
	  cosmphi = cos(pm * (PI + (j+0.5) * dphi));
	  sinmphi = sin(pm * (PI + (j+0.5) * dphi));
	  psi4RR = shellf[InList*n  ];
	  psi4II =-shellf[InList*n+1];
	 }
	  
	 thetap = sqrt((2*pl+1.0)/4.0/PI)*misc::Wigner_d_function(pl,pm,spinw,costheta); //note the variation from -2 to 2
#ifdef GaussInt
// wtcostheta is even function respect costheta
         RP_out[countlm] = RP_out[countlm] + thetap * (psi4RR * cosmphi + psi4II * sinmphi)*wtcostheta[i];
  	 IP_out[countlm] = IP_out[countlm] + thetap * (psi4II * cosmphi - psi4RR * sinmphi)*wtcostheta[i];
#else	 
         RP_out[countlm] = RP_out[countlm] + thetap * (psi4RR * cosmphi + psi4II * sinmphi);
  	 IP_out[countlm] = IP_out[countlm] + thetap * (psi4II * cosmphi - psi4RR * sinmphi);
#endif	 
	}
	countlm++;  //no sanity check for countlm and NN which should be noted in the input parameters
	  }
     }

    for(int ii=0;ii<NN;ii++)
    {
#ifdef GaussInt
      RP_out[ii] = RP_out[ii] * rex * dphi;
      IP_out[ii] = IP_out[ii] * rex * dphi;
#else	
      RP_out[ii] = RP_out[ii] * rex * dphi * dcostheta;
      IP_out[ii] = IP_out[ii] * rex * dphi * dcostheta;
#endif      
    }
//|------+  Communicate and sum the results from each processor.

    MPI_Allreduce(RP_out,RP,NN,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(IP_out,IP,NN,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

//|------= Free memory.

    delete[] pox[0];
    delete[] pox[1];
    delete[] pox[2];
    delete[] shellf;
    delete[] RP_out; delete[] IP_out;
    DG_List->clearList();
   }
//|----------------------------------------------------------------
//  for shell patch
//|----------------------------------------------------------------
  void surface_integral::surf_Wave(double rex,int lev,ShellPatch *GH, var *Rpsi4, var *Ipsi4,
		                   int spinw,int maxl, int NN, double *RP, double *IP,
				   monitor *Monitor)  //NN is the length of RP and IP
   {
    const int InList = 2;

    MyList<var> * DG_List=new MyList<var>(Rpsi4);
    DG_List->insert(Ipsi4);

    int n;
    double *pox[3];
    for(int i=0;i<3;i++) pox[i] = new double[n_tot];
    for( n = 0; n < n_tot; n++)
       {
        pox[0][n] = rex * nx_g[n];
        pox[1][n] = rex * ny_g[n];
        pox[2][n] = rex * nz_g[n];
       }

    double *shellf;
    shellf=new double[n_tot*InList];

    GH->Interp_Points(DG_List,n_tot,pox,shellf,Symmetry);    

    int mp, Lp, Nmin, Nmax;

    mp = n_tot / cpusize;
    Lp = n_tot - cpusize * mp;

    if( Lp > myrank ) { Nmin = myrank * mp + myrank; Nmax = Nmin + mp    ;}
    else              { Nmin = myrank * mp + Lp;     Nmax = Nmin + mp - 1;}

//|~~~~~> Integrate the dot product of Dphi with the surface normal.

    double *RP_out, *IP_out;
    RP_out = new double[NN];
    IP_out = new double[NN];    
    
    for(int ii=0;ii<NN;ii++)
    {
      RP_out[ii] = 0;
      IP_out[ii] = 0;
    }
// theta part    
    double costheta, thetap;
    double cosmphi,sinmphi;

    int i,j;
    int lpsy=0;
         if( Symmetry == 0 )     lpsy=1;
    else if( Symmetry == 1 )     lpsy=2;
    else if( Symmetry == 2 )     lpsy=8;

    double psi4RR,psi4II;
    for( n = Nmin; n <= Nmax; n++) 
     {
//       need round off always	     
        i = int(n/N_phi); // int(1.723) = 1, int(-1.732) = -1
        j = n - i * N_phi;
        
	int countlm=0;
	for(int pl=spinw;pl<maxl+1;pl++)
          for(int pm=-pl;pm<pl+1;pm++)
	  {
	for(int lp=0;lp<lpsy;lp++)
	{
 	 switch(lp)
	 {
	  case 0:  //+++ (theta, phi)
          costheta = arcostheta[i];
	  cosmphi = cos(pm * (j+0.5) * dphi);
	  sinmphi = sin(pm * (j+0.5) * dphi);
	  psi4RR = shellf[InList*n  ];
	  psi4II = shellf[InList*n+1];
	  break;
	  case 1:  //++- (pi-theta, phi)
          costheta = -arcostheta[i];
	  cosmphi = cos(pm * (j+0.5) * dphi);
	  sinmphi = sin(pm * (j+0.5) * dphi);
	  psi4RR = shellf[InList*n  ];
	  psi4II =-shellf[InList*n+1];
	  break;
	  case 2:  //+-+ (theta, 2*pi-phi)
          costheta = arcostheta[i];
	  cosmphi = cos(pm * (j+0.5) * dphi);
	  sinmphi =-sin(pm * (j+0.5) * dphi);
	  psi4RR = shellf[InList*n  ];
	  psi4II =-shellf[InList*n+1];
	  break;
	  case 3:  //+-- (pi-theta, 2*pi-phi)
          costheta = -arcostheta[i];
	  cosmphi = cos(pm * (j+0.5) * dphi);
	  sinmphi =-sin(pm * (j+0.5) * dphi);
	  psi4RR = shellf[InList*n  ];
	  psi4II = shellf[InList*n+1];
	  break;
	  case 4:  //-++ (theta, pi-phi)
          costheta = arcostheta[i];
	  cosmphi = cos(pm * (PI - (j+0.5) * dphi));
	  sinmphi = sin(pm * (PI - (j+0.5) * dphi));
	  psi4RR = shellf[InList*n  ];
	  psi4II =-shellf[InList*n+1];
	  break;
	  case 5:  //-+- (pi-theta, pi-phi)
          costheta = -arcostheta[i];
 	  cosmphi = cos(pm * (PI - (j+0.5) * dphi));
	  sinmphi = sin(pm * (PI - (j+0.5) * dphi));
	  psi4RR = shellf[InList*n  ];
	  psi4II = shellf[InList*n+1];
	  break;
	  case 6:  //--+ (theta, pi+phi)
          costheta = arcostheta[i];
	  cosmphi = cos(pm * (PI + (j+0.5) * dphi));
	  sinmphi = sin(pm * (PI + (j+0.5) * dphi));
	  psi4RR = shellf[InList*n  ];
	  psi4II = shellf[InList*n+1];
	  break;
	  case 7:  //--- (pi-theta, pi+phi)
          costheta = -arcostheta[i];
	  cosmphi = cos(pm * (PI + (j+0.5) * dphi));
	  sinmphi = sin(pm * (PI + (j+0.5) * dphi));
	  psi4RR = shellf[InList*n  ];
	  psi4II =-shellf[InList*n+1];
	 }
	  
	 thetap = sqrt((2*pl+1.0)/4.0/PI)*misc::Wigner_d_function(pl,pm,spinw,costheta); //note the variation from -2 to 2
#ifdef GaussInt
// wtcostheta is even function respect costheta
         RP_out[countlm] = RP_out[countlm] + thetap * (psi4RR * cosmphi + psi4II * sinmphi)*wtcostheta[i];
  	 IP_out[countlm] = IP_out[countlm] + thetap * (psi4II * cosmphi - psi4RR * sinmphi)*wtcostheta[i];
#else	 
         RP_out[countlm] = RP_out[countlm] + thetap * (psi4RR * cosmphi + psi4II * sinmphi);
  	 IP_out[countlm] = IP_out[countlm] + thetap * (psi4II * cosmphi - psi4RR * sinmphi);
#endif	 
	}
	countlm++;  //no sanity check for countlm and NN which should be noted in the input parameters
	  }
     }

    for(int ii=0;ii<NN;ii++)
    {
#ifdef GaussInt
      RP_out[ii] = RP_out[ii] * rex * dphi;
      IP_out[ii] = IP_out[ii] * rex * dphi;
#else	
      RP_out[ii] = RP_out[ii] * rex * dphi * dcostheta;
      IP_out[ii] = IP_out[ii] * rex * dphi * dcostheta;
#endif      
    }
//|------+  Communicate and sum the results from each processor.

    MPI_Allreduce(RP_out,RP,NN,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(IP_out,IP,NN,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

//|------= Free memory.

    delete[] pox[0];
    delete[] pox[1];
    delete[] pox[2];
    delete[] shellf;
    delete[] RP_out; delete[] IP_out;
    DG_List->clearList();
   }
