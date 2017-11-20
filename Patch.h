//$Id: Patch.h,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
#ifndef PATCH_H
#define PATCH_H

#include <mpi.h>
#include "MyList.h"
#include "Block.h"
#include "var.h"
#include "microdef.h" //need dim here; Vertex or Cell; ghost_width

class Patch{

public:
   
   int lev;
   int shape[dim];
   double bbox[2*dim];  //this bbox includes buffer points
   MyList<Block> *blb, *ble;
   int lli[dim],uui[dim];  //denote the buffer points on each boundary

public:
  
   Patch(){}; 
   Patch(int DIM, int *shapei, double *bboxi,int levi,bool buflog,int Symmetry);

   ~Patch();

   void checkPatch(bool buflog);
   void checkBlock();
   void Interp_Points(MyList<var> *VarList,
		int NN,double **XX,
		double *Shellf,int Symmetry);
   double getdX(int dir);
};

#endif   /* PATCH_H */
