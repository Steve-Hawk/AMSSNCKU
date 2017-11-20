//$Id: ri2aphi.C,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
#ifdef newc
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstring>
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

#define PI M_PI

// usage:
// ri2aphi infile outfile
// transform data in infile from real part and imaginary part form to ampulitude and phase form
// and store in outfile
int main(int argc, char * argv[])
{
    for(int i = 1; i < argc; i++)
      if(!(strncmp("-h",argv[i],2)))
      {
         cout<<" usage:"<<endl<<"ri2aphi infile outfile"<<endl;
	 cout<<"transform data in infile from real part and imaginary part form to ampulitude and phase form"<<endl;
         cout<<"and store in outfile"<<endl;
         exit(0);
      }
    if(argc!=3)
    {
      cout<<"syntax erro, type 'pickdata -h' for help"<<endl;
      exit(0);
    }

    ifstream infile;
    ofstream outfile;
    infile.open(argv[1]);
    outfile.open(argv[2]);
    if(!infile)
    {
        cerr << "Can't open " << argv[1] << " for input." << endl;
	exit(1);
    }
    if(!outfile)
    {
        cerr << "Can't open " << argv[2] << " for output." << endl;
	exit(1);
    }

// deal with the three lines head information
    int headi = 0,NN;

    char c[10000]; infile.getline(c,10000);
    while(int(c[0])<48 || int(c[0])>57)
    {
       headi++; infile.getline(c,10000);
    }
// find out the number of colomns    
    {
      char *c1=c;
      NN=1;
      while(*(c1+1))
      {
//                 0               9                   " "	      
	if(int(*c1)>47 && int(*c1)<58 && int(*(c1+1))==32) NN++;
	c1++;
      }
    }

    infile.close();
    if((NN-1)%2 != 0)
    {
	cout<<NN<<" colomns? For time, real part and imaginary part?"<<endl;
        outfile.close();
	exit(0);
    }
    infile.open(argv[1]);
    for(int i=0;i<headi;i++) {infile.getline(c,10000); outfile<<c<<endl;}

    int N = (NN-1)/2;

    double *phio,*phin,*tphi,*am;
    phio = new double[N];
    phin = new double[N];
    tphi = new double[N];
    am   = new double[N];
    double t,*a;
    a = new double[NN-1];

    if(infile.eof()) {cout<<"empty "<<argv[2]<<endl; exit(0);}

    void convertaphi(double rp,double ip, double &am, double &phi);

    infile>>t;
    outfile<<setprecision(10)<<t<<" ";
    for(int i=0;i<NN-1;i++) infile>>a[i];
    for(int i=0;i<N;i++)
    {
       convertaphi(a[2*i],a[2*i+1],am[i],phin[i]);
       tphi[i]=phin[i];
       if(i<N-1)
       {
            outfile<<setw(10)<<setprecision(10)<<am[i]<<" "
	           <<setw(10)<<setprecision(10)<<tphi[i]<<" ";
       }
       else
       {
            outfile<<setw(10)<<setprecision(10)<<am[i]<<" "
	           <<setw(10)<<setprecision(10)<<tphi[i];
       }
    }
    outfile<<endl;

    while(!infile.eof())
    {
	double tmp;
        infile>>t;
	if(infile.eof()) break;
        outfile<<setprecision(10)<<t<<" ";
        for(int i=0;i<NN-1;i++) infile>>a[i];
	for(int i=0;i<N;i++)
	{
	  phio[i]=phin[i];
	  convertaphi(a[2*i],a[2*i+1],am[i],phin[i]);
	  tmp = phin[i]-phio[i];
	  tphi[i] += tmp;
	  if     (tmp < -PI/2) tphi[i] += PI;
	  else if(tmp >  PI/2) tphi[i] -= PI;
	  if(i<N-1)
	  {
            outfile<<setw(10)<<setprecision(10)<<am[i]<<" "
	           <<setw(10)<<setprecision(10)<<tphi[i]<<" ";
	  }
	  else
	  {
            outfile<<setw(10)<<setprecision(10)<<am[i]<<" "
	           <<setw(10)<<setprecision(10)<<tphi[i];
	  }
        }
        outfile<<endl;
    }

      infile.close();
      outfile.close();

      delete[] phin; delete[] phio; delete[] tphi; 
      delete[] am; delete[] a;

      return 0;
}
 void convertaphi(double rp,double ip, double &am, double &phi)
{
   am  = sqrt(rp*rp+ip*ip);
   phi = atan(ip/rp);
}
