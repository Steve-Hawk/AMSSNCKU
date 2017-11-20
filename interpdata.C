//$Id: interpdata.C,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
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
// usage:
// interpdata tf dt file1 file2 ord
// interpolation of data from file1 for tf, tf+dt, tf+2dt, ......
// and store in file2
int main(int argc, char * argv[])
{
    for(int i = 1; i < argc; i++)
      if(!(strncmp("-h",argv[i],2)))
      {
         cout<<" usage:"<<endl<<"interpdata tf dt file1 file2 ord"<<endl;
	 cout<<"interpolation with order ord of data from file1 for tf, tf+dt, tf+2dt, ......"<<endl;
         cout<<"and store in file2"<<endl;
         exit(0);
      }
    if(argc!=6)
    {
      cout<<"syntax erro, type 'pickdata -h' for help"<<endl;
      exit(0);
    }

    double tf,dt;
    tf=atof(argv[1]);
    dt=atof(argv[2]);
    int ORD;
    ORD=atoi(argv[5]);

    ifstream infile;
    ofstream outfile;
    infile.open(argv[3]);
    outfile.open(argv[4]);
    if(!infile)
    {
      cerr << "Can't open " << argv[3] << " for input." << endl;
      exit(1);
    }
    if(!outfile)
    {
      cerr << "Can't open " << argv[4] << " for output." << endl;
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
	if(int(*c1)>47 && int(*c1)<58 && int(*(c1+1))==32) NN++;
	c1++;
      }
    }

    infile.close();
    infile.open(argv[3]);
    for(int i=0;i<headi;i++) {infile.getline(c,10000); outfile<<c<<endl;}

    double *ti,**a;
    ti = new double[ORD];
    a = new double*[NN-1];
    for(int i=0;i<NN-1;i++) a[i] = new double[ORD];

    int ncount=0;
    for(int j=0;j<ORD;j++)
    {     
     if(infile.eof())
     {
	cout<<"end of file "<<argv[3]<< " without any output"<<endl;
	exit(0);
     }
     infile>>ti[j];
     for(int i=0;i<NN-1;i++) infile>>a[i][j];
    }
// Let me call this two runners method
// runner 1:
    while(tf+ncount*dt<ti[ORD/2]) ncount++;

    double ar[NN-1];

    while(!infile.eof())
    {
// runner 2:	    
      while(tf+ncount*dt > ti[ORD/2+1])
      {
	for(int i=1;i<ORD;i++)
	{
          ti[i-1]=ti[i];
	  for(int j=0;j<NN-1;j++) a[j][i-1]=a[j][i];
	}
        if(infile.eof()) break;
        infile>>ti[ORD-1];
        for(int i=0;i<NN-1;i++) infile>>a[i][ORD-1];
      }
      
      if(infile.eof()) break;

      void polyinterp(double t,double &rr,double *ti, double *ri,const int ORD);
      for(int j=0;j<NN-1;j++) polyinterp(tf+ncount*dt,ar[j],ti,a[j],ORD);

      outfile<<setw(10) << setprecision(10)<<tf+ncount*dt<<" ";
      for(int j=0;j<NN-1;j++) outfile<<setw(10) << setprecision(10)<<ar[j]<<" ";
      outfile<<endl;

      ncount++;
    }

      infile.close();
      outfile.close();

      delete[] ti; 
      for(int i=0;i<NN-1;i++) delete[] a[i];
      delete[] a;

      return 0;
}
 void polyinterp(double t,double &rr,double *ti, double *ri,const int ORD)
{
//  (x  -x_1)...(x  -x_i-1)(x  -x_i+1)...(x  -x_N)
// ------------------------------------------------f_i
//  (x_i-x_1)...(x_i-x_i-1)(x_i-x_i+1)...(x_i-x_N)

   rr=0;
   for(int i=0;i<ORD;i++)
   {
     double ss=1,xx=1;
     for(int j=0;j<ORD;j++)
     {
       if(j != i)
       {
        ss*=t-ti[j];
	xx*=ti[i]-ti[j];
       }
     }
     rr+=ss/xx*ri[i];
   }

}
