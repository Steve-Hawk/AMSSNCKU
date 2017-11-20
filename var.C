//$Id: var.C,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
using namespace std;

#include <time.h>
#include <mpi.h>

#include "var.h"

var::var(char *namei, int sgfni,
	 const double SYM1, const double SYM2, const double SYM3):
	sgfn(sgfni)
   {
     char *p=namei;
     int i=0;
     while(*(p++)) i++;
     if(i>20) cout<<"too long name for var: "<<namei<<endl;
     sprintf(name,namei);
     SoA[0]=SYM1;
     SoA[1]=SYM2;
     SoA[2]=SYM3;
   }

var::~var(){}
