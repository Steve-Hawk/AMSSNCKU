//$Id: monitor.C,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
#ifdef newc
#include <cstdio>
using namespace std;
#else
#include <stdio.h>
#endif

#include "monitor.h"

monitor::monitor(const char fname[],int myrank,string head):outfile(0)
    {
      I_Print = (myrank == 0);

      if (I_Print)
      {
        outfile.open(fname,ios::trunc);
	
	time_t tnow;
	time(&tnow);
	struct tm *loc_time;
	loc_time = localtime(&tnow);

	outfile<< "# File created on " << asctime(loc_time);
	outfile<<"#"<<endl;
	outfile.setf(ios::left);
	outfile<<head<<endl;
      }
    }
monitor:: ~monitor() 
  {
   if (I_Print) outfile.close();
  }
void monitor::writefile(double time,int NN,double *DDAT)
    {
      if (I_Print) 
      {
	outfile << setprecision(8);
	outfile << setw(14) << time;
	for(int countlm=0;countlm<NN;countlm++)
	{
	   outfile<<" "<< setw(15) << DDAT[countlm];
	}
	outfile << endl;
      }
    }
void monitor::writefile(double time,int NN,double *DDAT1,double *DDAT2)
    {
      if (I_Print) 
      {
	outfile << setprecision(8);
	outfile << setw(14) << time;
	for(int countlm=0;countlm<NN;countlm++)
	{
	   outfile<<" "<< setw(15) << DDAT1[countlm]
		  <<" "<< setw(15) << DDAT2[countlm];
	}
	outfile << endl;
      }
    }
