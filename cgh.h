//$Id: cgh.h,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
#ifndef CGH_H
#define CGH_H

#include <mpi.h>
#include "MyList.h"
#include "Patch.h"
#include "microdef.h"
#include "monitor.h"

class cgh{

public:

   int levels,movls;
// information of boxes   
   int *grids; 
   double ***bbox;
   int ***shape;
   double ***handle;

// information of Patch list
   MyList<Patch> **PatL;

protected:
   int ingfs,fngfs;   
   static const double ratio=0.8;
public:
   
   cgh(int ingfsi,int fngfsi,int Symmetry,char *filename,int checkrun,monitor *ErrorMonitor);

   ~cgh();

   void compose_cgh(int nprocs);
   void sethandle(monitor *ErrorMonitor);
   void checkPatchList(MyList<Patch> *PatL,bool buflog);
   void Regrid(int Symmetry,int BH_num,double **Porgbr,double**Porg0,
		 MyList<var> *OldList, MyList<var> *StateList, 
		 MyList<var> *FutureList,MyList<var> *tmList,bool BB,
		 monitor *ErrorMonitor);
   void recompose_cgh(int nprocs,bool *lev_flag,
		        MyList<var> *OldList, MyList<var> *StateList, 
			MyList<var> *FutureList,MyList<var> *tmList,
			int Symmetry, bool BB);
   void read_bbox(int Symmetry,char *filename);
   MyList<Patch> *construct_patchlist(int lev,int Symmetry);
};

#endif   /* CGH_H */
