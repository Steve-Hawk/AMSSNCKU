//$Id: misc.C,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
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

#define PI M_PI

void misc::tillherecheck(int myrank)
{
 int atp=1,tatp;
 MPI_Allreduce(&atp,&tatp,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
 if(myrank==0) cout<<" here now: "<<tatp<<" processors."<<endl;
}
void misc::tillherecheck(const char str[])
{
 int myrank;
 MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
 int atp=1,tatp;
 MPI_Allreduce(&atp,&tatp,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
 if(myrank==0)
 {
   cout<<" here now: "<<tatp<<" processors."<<endl;
   cout<<str<<endl;
 }
}
// pick out value from input string
int misc::parse_parts(string str, string& sgrp, string& skey, string& sval, int& ind)
{
      int pos1, pos2;
      string s0;
      
      ind = 0;
      
      // remove comments
      str = str.substr(0, str.find("#") );
      if ( rTrim(str).empty() ) return 0;   // continue;
      
      // parse {group, key, val}
      pos1 = str.find("::");  pos2 = str.find("=");
      if (pos1 == string::npos || pos2 == string::npos) return -1;
      
      s0 = str.substr(     0, pos1       );  sgrp = lTrim( s0 );
      s0 = str.substr(pos1+2, pos2-pos1-2);  skey = rTrim( s0 );
      s0 = str.substr(pos2+1)             ;  sval = Trim( s0 );  

      pos1 = sval.find("\"");  pos2 = sval.rfind("\"");
      if ( pos1 != string::npos ) {
         sval = sval.substr(1, pos2-1); 
      }

      pos1 = skey.find("[");  pos2 = skey.find("]");
      if ( pos1 != string::npos ) {
         s0   = skey.substr(0,pos1);
         ind = atoi( skey.substr(pos1+1 ,pos2-pos1-1).c_str() );
         skey  = s0; 
      }
      
      return 1;
}
int misc::parse_parts(string str, string& sgrp, string& skey, string& sval, int& ind1, int& ind2)
{
      int pos1, pos2;
      string s0,s1;
      
      ind1 = ind2 = 0;
      
      // remove comments
      str = str.substr(0, str.find("#") );
      if ( rTrim(str).empty() ) return 0;   // continue;
      
      // parse {group, key, val}
      pos1 = str.find("::");  pos2 = str.find("=");
      if (pos1 == string::npos || pos2 == string::npos) return -1;
      
      s0 = str.substr(     0, pos1       );  sgrp = lTrim( s0 );
      s0 = str.substr(pos1+2, pos2-pos1-2);  skey = rTrim( s0 );
      s0 = str.substr(pos2+1)             ;  sval = Trim( s0 );  

      pos1 = sval.find("\"");  pos2 = sval.rfind("\"");
      if ( pos1 != string::npos ) {
         sval = sval.substr(1, pos2-1); 
      }
      
      pos1 = skey.find("[");  pos2 = skey.find("]");
      if ( pos1 != string::npos ) {
         s0   = skey.substr(0,pos1);
	 s1   = skey.substr(pos2+1);
         ind1 = atoi( skey.substr(pos1+1 ,pos2-pos1-1).c_str() );
         skey  = s0; 
      }

      pos1 = s1.find("[");  pos2 = s1.find("]");
      if ( pos1 != string::npos ) {
         s0   = s1.substr(pos2+1);
         ind2 = atoi( s1.substr(pos1+1 ,pos2-pos1-1).c_str() );
      }
      
      return 1;
}
int misc::parse_parts(string str, string& sgrp, string& skey, string& sval, int& ind1, int& ind2, int& ind3)
{
      int pos1, pos2;
      string s0,s1;
      
      ind1 = ind2 = ind3 = 0;
      
      // remove comments
      str = str.substr(0, str.find("#") );
      if ( rTrim(str).empty() ) return 0;   // continue;
      
      // parse {group, key, val}
      pos1 = str.find("::");  pos2 = str.find("=");
      if (pos1 == string::npos || pos2 == string::npos) return -1;
      
      s0 = str.substr(     0, pos1       );  sgrp = lTrim( s0 );
      s0 = str.substr(pos1+2, pos2-pos1-2);  skey = rTrim( s0 );
      s0 = str.substr(pos2+1)             ;  sval = Trim( s0 );  

      pos1 = sval.find("\"");  pos2 = sval.rfind("\"");
      if ( pos1 != string::npos ) {
         sval = sval.substr(1, pos2-1); 
      }
      
      pos1 = skey.find("[");  pos2 = skey.find("]");
      if ( pos1 != string::npos ) {
         s0   = skey.substr(0,pos1);
	 s1   = skey.substr(pos2+1);
         ind1 = atoi( skey.substr(pos1+1 ,pos2-pos1-1).c_str() );
         skey  = s0; 
      }

      pos1 = s1.find("[");  pos2 = s1.find("]");
      if ( pos1 != string::npos ) {
         s0   = s1.substr(pos2+1);
         ind2 = atoi( s1.substr(pos1+1 ,pos2-pos1-1).c_str() );
      }

      pos1 = s0.find("[");  pos2 = s0.find("]");
      if ( pos1 != string::npos ) {
         ind3 = atoi( s0.substr(pos1+1 ,pos2-pos1-1).c_str() );
      }
      
      return 1;
}
// sent me from Roman Gold on 2010-10-8
void misc::gaulegf(double x1, double x2, double *x, double *w, int n)
  {
    int i, j, m;
    double eps = 1.2E-16;
    double p1, p2, p3, pp, xl, xm, z, z1;
  
    m = (n+1)/2;
    xm = 0.5*(x2+x1);
    xl = 0.5*(x2-x1);
    for(i=0; i<m; i++)
    {
      z = cos(PI*((double)i+0.75)/((double)n+0.5));
      do
      {
        p1 = 1.0;
        p2 = 0.0;
        for(j=0; j<n; j++)
        {
          p3 = p2;
          p2 = p1;
          p1 = ((2*(double)j+1)*z*p2-(double)j*p3)/((double)j+1);
        }
        pp = n*(z*p1-p2)/(z*z-1.0);
        z1 = z;
        z = z1 - p1/pp;
      }while(fabs(z-z1) > eps);
      x[i] = xm - xl*z;
      x[n-1-i] = xm + xl*z;
      w[i] = 2.0*xl/((1.0-z*z)*pp*pp);
      w[n-1-i] = w[i];
    }
  } /* end gaulegf */
void misc::inversearray(double *aa,int NN)
{
    int i, m;
    m = (NN+1)/2;
    double rr;
    for(i=0; i<m; i++)
    {
       rr = aa[i];
       aa[i] = aa[NN-1-i];
       aa[NN-1-i] = rr;
    }
}
//Eq.(42) of PRD 77, 024027 (2008)
double misc::Wigner_d_function(int l,int m, int s, double costheta)
{
//we consider only theta in [0,pi]
  int C1=max(0,m-s),C2=min(l+m,l-s);
  double vv=0;
  double sinht=sqrt((1-costheta)/2.0),cosht=sqrt((1+costheta)/2.0);
  if(C1%2==0)
  {
    for(int t=C1;t<C2+1;t+=2)
       vv=vv+pow(cosht,2*l+m-s-2*t)*pow(sinht,2*t+s-m)/
	     (fact(l+m-t)*fact(l-s-t)*fact(t)*fact(t+s-m));
    for(int t=C1+1;t<C2+1;t+=2)
       vv=vv-pow(cosht,2*l+m-s-2*t)*pow(sinht,2*t+s-m)/
	     (fact(l+m-t)*fact(l-s-t)*fact(t)*fact(t+s-m));
  }
  else
  {
    for(int t=C1;t<C2+1;t+=2)
       vv=vv-pow(cosht,2*l+m-s-2*t)*pow(sinht,2*t+s-m)/
	     (fact(l+m-t)*fact(l-s-t)*fact(t)*fact(t+s-m));
    for(int t=C1+1;t<C2+1;t+=2)
       vv=vv+pow(cosht,2*l+m-s-2*t)*pow(sinht,2*t+s-m)/
	     (fact(l+m-t)*fact(l-s-t)*fact(t)*fact(t+s-m));
  }
  return vv*sqrt(fact(l+m)*fact(l-m)*fact(l+s)*fact(l-s));
}
double misc::fact(int N)
{
  if(N < 0) cout<<"error input for factorial."<<endl;
  double f;
  if(N==0)f=1;
  else    f=N*fact(N-1);
  return f;
}
