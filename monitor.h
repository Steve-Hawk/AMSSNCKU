//$Id: monitor.h,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
#ifndef MONITOR_H
#define MONITOR_H

#ifdef newc
#include <iostream>
#include <iomanip>
#include <strstream>
#include <fstream>
#include <string>
using namespace std;
#else
#include <iostream.h>
#include <iomanip.h>
#include <strstream>
#include <fstream.h>
#include <string.h>
#endif
#include <time.h>

#include <mpi.h>

class monitor {

public:
  ofstream outfile;

private:
  bool I_Print;

public:

  monitor(const char fname[],int myrank,string head);

  ~monitor();

  void writefile(double time, int NN, double *DDAT);
  void writefile(double time,int NN,double *DDAT1,double *DDAT2);
};

#endif   /* MONITOR */
