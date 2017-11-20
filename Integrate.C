//$Id: Integrate.C,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
#ifdef newc
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cstdio>
using namespace std;
#else
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#endif

int main(int argc, char * argv[])
{
    for(int i = 1; i < argc; i++)
      if(!(strncmp("-h",argv[i],2)))
      {
         cout<<" usage:"<<endl<<"Integrate file ord"<<endl;
	 cout<<"sum up data for every collomn, but the first collomn which is time, of file with ord-th order scheme"<<endl;
         cout<<"and store in sum_file"<<endl;
         exit(0);
      }
    if(argc!=3)
    {
      cout<<"syntax erro, type 'doingsum -h' for help"<<endl;
      exit(0);
    }
    int ORD;
    ORD=atoi(argv[2]);

      ifstream infile;
      ofstream outfile;
      infile.open(argv[1]);
      if(!infile)
       {
        cerr << "Can't open " << argv[1] << " for input." << endl;
	exit(1);
       }      

      char filename[50];
      sprintf(filename,"sum_%s",argv[1]);
      outfile.open(filename);
      if(!outfile)
       {
        cerr << "Can't open " << filename << " for output." << endl;
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
//                 0               9                   " "	      TAB   
	if(int(*c1)>47 && int(*c1)<58 && (int(*(c1+1))==32 || int(*(c1+1))==9)) NN++;
	c1++;
      }
    }
    cout<<"find "<<NN<<" columns in file "<<argv[1]<<endl;

    infile.close();
    infile.open(argv[1]);
    for(int i=0;i<headi;i++) {infile.getline(c,10000); outfile<<c<<endl;}
 
      double dT,tmp;
      for(int j=0;j<2;j++) 
      for(int i=0;i<NN;i++)
	{
	  infile>>tmp;
          if(infile.eof())
	  {
	    cout<<"too few data to be sum!"<<endl;
	    exit(0);
	  }
	  if(j==0 && i==0) dT=tmp;
	  if(j==1 && i==0) dT=tmp-dT;
	}
      cout<<"get dT = "<<dT<<endl;
      infile.close();

      infile.open(argv[1]);
      for(int i=0;i<headi;i++) {infile.getline(c,10000); outfile<<c<<endl;}

      int i;
      double **a,*sum;
      double *coef;
      if(ORD<1){cout<<"bad order parameter: "<<ORD<<endl; exit(0);}
      a  = new double*[ORD];
      coef = new double[ORD];
      switch(ORD)
      {
	case(1):
	coef[0] = 1;
	break;
	case(2):
	coef[0] = coef[1] = 0.5;
	break;
	case(3):
	coef[0] = coef[2] = 1.0/3;
	coef[1] = 4.0/3;
	break;
	case(4):
	coef[0] = coef[3] = 3.0/8;
	coef[1] = coef[2] = 9.0/8;
	break;
	case(5):
	coef[0] = coef[4] = 14.0/45;
	coef[1] = coef[3] = 64.0/45;
	coef[2] = 24.0/45;
	break;
	case(6):
	coef[0] = coef[5] = 95.0/288;
	coef[1] = coef[4] = 375.0/288;
	coef[2] = coef[3] = 250.0/288;
	break;
	case(7):
	coef[0] = coef[6] = 41.0/140;
	coef[1] = coef[5] = 216.0/140;
	coef[2] = coef[4] = 27.0/140;
	coef[3] = 272.0/140;
	break;
	case(8):
	coef[0] = coef[7] = 5257.0/17280;
	coef[1] = coef[6] = 25039.0/17280;
	coef[2] = coef[5] = 9261.0/17280;
	coef[3] = coef[4] = 20923.0/17280;
	break;

	default:
	cout<<"we do not support this order integration at present! ORD = "<<ORD<<endl;
        exit(0);
      }
      for(i=0;i<ORD;i++) a[i] = new double[NN];
      sum=new double[NN];

      for(i=0;i<NN;i++)sum[i]=0.0;
      
      for(i=0;i<NN;i++)
      {
        if(infile.eof())
	  {
	    cout<<"too few data to be sum!"<<endl;
	    exit(0);
	  }
        infile>>a[ORD-1][i];
      }

      cout<<"integrating "<< argv[1] <<" with "<<ORD<<"-th order scheme ..."<<endl;
      while(!infile.eof())
      {
        if(ORD > 1)for(i=0;i<NN;i++) a[0][i] = a[ORD-1][i];    //starting point for integration, head-tail connecting
	for(int ordi=1;ordi<ORD;ordi++)
	{
         for(i=0;i<NN;i++)
	 {
          if(infile.eof()) 
	  {
            infile.close();
            outfile.close();

            for(i=0;i<ORD;i++) delete[] a[i];
            delete[] a; delete[] coef;
            delete[] sum;

	    exit(0);
	  }
	  infile>>a[ordi][i];
	 }
	}
        for(int ordi=0;ordi<ORD;ordi++) 
	{
	   sum[0] = a[ordi][0];
           for(i=1;i<NN;i++)
	    sum[i]=sum[i]+coef[ordi]*a[ordi][i]*dT;  
	}

	outfile<<setw(16) << setprecision(12)<<a[ORD-1][0]<<" ";
        for(i=1;i<NN;i++)   outfile<<setw(16) << setprecision(12)<<a[ORD-1][i]<<" "<<setw(16) << setprecision(12)<<sum[i];
	outfile<<endl;
//  special treatment for ORD == 1	
	for(int ordi=ORD-1;ordi<1;ordi++)
	{
         for(i=0;i<NN;i++)
	 {
          if(infile.eof()) 
	  {
            infile.close();
            outfile.close();

            for(i=0;i<ORD;i++) delete[] a[i];
            delete[] a; delete[] coef;
            delete[] sum;

	    exit(0);
	  }
	  infile>>a[ordi][i];
	 }
	}
      }

      infile.close();
      outfile.close();

      for(i=0;i<ORD;i++) delete[] a[i];
      delete[] a; delete[] coef;
      delete[] sum;

      return 0;
}
