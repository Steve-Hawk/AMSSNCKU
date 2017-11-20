//$Id: Block.h,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
#ifndef BLOCK_H
#define BLOCK_H

#include <mpi.h>
#include "microdef.h" //need dim here; Vertex or Cell
#include "var.h"
#include "MyList.h"
class Block{

public:
   int shape[dim];
   double bbox[2*dim];
   double *X[dim];
   int rank;
   int lev;
   int ingfs,fngfs;
   int *(*igfs);
   double *(*fgfs);

public:
   
   Block(int DIM, int *shapei, double *bboxi,int ranki, int ingfsi, int fngfs, int levi);

   ~Block();

   void checkBlock();

   double getdX(int dir);
   void swapList(MyList<var> *VarList1,MyList<var> *VarList2,int myrank);
};

#endif   /* BLOCK_H */
