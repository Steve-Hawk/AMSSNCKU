//$Id: ShellPatch.C,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <new>
using namespace std;

#include "ShellPatch.h"
#include "Parallel.h"
#include "fmisc.h"
#include "misc.h"
#include "shellfunctions.h"

#define PI M_PI

// x   x   x   x   x   o   *
//             *   o   x   x   x   x   x
// each side contribute an overlap points
// so we need half of that
#define overghost ((ghost_width+1)/2+ghost_width)

ss_patch::ss_patch(int ingfsi, int fngfsi,int *shapei,double *bboxi,int myranki):
ingfs(ingfsi),fngfs(fngfsi),myrank(myranki),blb(0),ble(0)
{
  for(int i=0;i<dim;i++) 
  {
      shape[i] = shapei[i];
      bbox[i] = bboxi[i];
      bbox[i+dim] = bboxi[i+dim];
  }
}
ss_patch::~ss_patch()
{
    MyList<Block> *bg;
    while(blb)
    {
      if(blb==ble) break;
      bg=(blb->next)?blb->next:0;
      delete blb->data;
      delete blb;
      blb=bg;
    }
    if(ble)
    {
      delete ble->data;
      delete ble;
    }
    blb=ble=0;
}
// bulk part for given Block within given patch, without extension
MyList<Parallel::gridseg> *ss_patch::build_bulk_gsl(Block *bp)
{
       MyList<Parallel::gridseg> *gs=0;

       gs = new MyList<Parallel::gridseg>;
       gs->data = new Parallel::gridseg;

       for(int i=0;i<dim;i++)
       {
         double DH = bp->getdX(i);
         gs->data->uub[i] = (feq(bp->bbox[dim+i],bbox[dim+i],DH/2)) ? bp->bbox[dim+i] : bp->bbox[dim+i]-ghost_width*DH;
         gs->data->llb[i] = (feq(bp->bbox[i]    ,bbox[i]    ,DH/2)) ? bp->bbox[i]     : bp->bbox[i]    +ghost_width*DH;
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
         gs->data->shape[i] = int((gs->data->uub[i] - gs->data->llb[i])/DH+0.4) + 1;
#else
#ifdef Cell
         gs->data->shape[i] = int((gs->data->uub[i] - gs->data->llb[i])/DH+0.4);
#else
#error Not define Vertex nor Cell
#endif  
#endif 
       }
       gs->data->Bg=bp;
       gs->next=0;

       return gs; 
}
// collect all ghost grid segments or blocks for given patch
MyList<Parallel::gridseg> *ss_patch::build_ghost_gsl()
{
       MyList<Parallel::gridseg> *cgsl=0,*gs,*gsb;
       MyList<Block> *BP=blb;
       while(BP)
       {
          gs = new MyList<Parallel::gridseg>;
	  gs->data = new Parallel::gridseg;
	  
	  for(int i=0;i<dim;i++)
	  {
	     gs->data->llb[i]=BP->data->bbox[i];
	     gs->data->uub[i]=BP->data->bbox[dim+i];
	     gs->data->shape[i]=BP->data->shape[i];
	  }
	  gs->data->Bg=BP->data;
	  gs->next=0;

          gsb = build_bulk_gsl(BP->data);

          if(!cgsl) cgsl = Parallel::gs_subtract(gs,gsb);
	  else      cgsl->catList(Parallel::gs_subtract(gs,gsb));

	  gsb->destroyList();
	  gs->destroyList();

	  if(BP == ble) break;
	  BP=BP->next;
       }

       return cgsl; 
}
// collect all grid segments or blocks without ghost for given patch
// special for Sync usage, so we do not need consider missing points
MyList<Parallel::gridseg> *ss_patch::build_owned_gsl0(int rank_in)
{
       MyList<Parallel::gridseg> *cgsl=0,*gs;
       MyList<Block> *BP=blb;
       while(BP)
       {
	  Block *bp=BP->data;
	  if(bp->rank == rank_in)
	  {
           if(!cgsl)
	   {
	     cgsl = gs = new MyList<Parallel::gridseg>;
	     gs->data = new Parallel::gridseg;
	   }
	   else
	   {
	     gs->next = new MyList<Parallel::gridseg>;
	     gs = gs->next;
	     gs->data = new Parallel::gridseg;
	   }
	    
	   for(int i=0;i<dim;i++)
	   {
	     double DH = bp->getdX(i);
             gs->data->uub[i] = (feq(bp->bbox[dim+i],bbox[dim+i],DH/2)) ? bp->bbox[dim+i] : bp->bbox[dim+i]-ghost_width*DH;
             gs->data->llb[i] = (feq(bp->bbox[i]    ,bbox[i]    ,DH/2)) ? bp->bbox[i]     : bp->bbox[i]    +ghost_width*DH;
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
	     gs->data->shape[i] = int((gs->data->uub[i] - gs->data->llb[i])/DH+0.4) + 1;
#else
#ifdef Cell
	     gs->data->shape[i] = int((gs->data->uub[i] - gs->data->llb[i])/DH+0.4);
#else
#error Not define Vertex nor Cell
#endif  
#endif 
	   }
	   gs->data->Bg=BP->data;
	   gs->next=0;
	  }

	  if(BP == ble) break;
	  BP=BP->next;
       }

       return cgsl; 
}
void ss_patch::Sync(MyList<var> *VarList,int Symmetry)
{
  int cpusize;
  MPI_Comm_size(MPI_COMM_WORLD,&cpusize);

  MyList<Parallel::gridseg> *dst;
  MyList<Parallel::gridseg> **src,**transfer_src,**transfer_dst;
  src=new MyList<Parallel::gridseg> *[cpusize];
  transfer_src=new MyList<Parallel::gridseg> *[cpusize];
  transfer_dst=new MyList<Parallel::gridseg> *[cpusize];

  dst=build_ghost_gsl();  // ghost region only
  for(int node=0;node<cpusize;node++)
  {
    src[node]=build_owned_gsl0(node); //for the part without ghost points and do not extend
    Parallel::build_gstl(src[node], dst, &transfer_src[node],&transfer_dst[node]); //for transfer[node], data locate on cpu#node
  }

  Parallel::transfer(transfer_src,transfer_dst,VarList,VarList,Symmetry);
  
  if(dst) dst->destroyList();
  for(int node=0;node<cpusize;node++)
  {
    if(src[node]) src[node]->destroyList();
    if(transfer_src[node]) transfer_src[node]->destroyList();
    if(transfer_dst[node]) transfer_dst[node]->destroyList();
  }

  delete[] src; delete[] transfer_src; delete[] transfer_dst;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void xp_patch::setupcordtrans()
{    
    MyList<Block> *BP=blb;
    while(BP)
    {
      Block *cg=BP->data;
      if(myrank == cg->rank) 
      {
        f_xp_getxyz(cg->shape,cg->X[0],cg->X[1],cg->X[2], 
                     cg->fgfs[fngfs+ShellPatch::gx],cg->fgfs[fngfs+ShellPatch::gy],cg->fgfs[fngfs+ShellPatch::gz]);
        f_xpm_getjacobian(cg->shape,cg->X[0],cg->X[1],cg->X[2], 
                     cg->fgfs[fngfs+ShellPatch::gx],cg->fgfs[fngfs+ShellPatch::gy],cg->fgfs[fngfs+ShellPatch::gz],
                     cg->fgfs[fngfs+ShellPatch::drhodx],cg->fgfs[fngfs+ShellPatch::drhody],cg->fgfs[fngfs+ShellPatch::drhodz],
                     cg->fgfs[fngfs+ShellPatch::dsigmadx],cg->fgfs[fngfs+ShellPatch::dsigmady],cg->fgfs[fngfs+ShellPatch::dsigmadz],
                     cg->fgfs[fngfs+ShellPatch::dRdx],cg->fgfs[fngfs+ShellPatch::dRdy],cg->fgfs[fngfs+ShellPatch::dRdz],
                     cg->fgfs[fngfs+ShellPatch::drhodxx],cg->fgfs[fngfs+ShellPatch::drhodxy],cg->fgfs[fngfs+ShellPatch::drhodxz],
                     cg->fgfs[fngfs+ShellPatch::drhodyy],cg->fgfs[fngfs+ShellPatch::drhodyz],cg->fgfs[fngfs+ShellPatch::drhodzz],
                     cg->fgfs[fngfs+ShellPatch::dsigmadxx],cg->fgfs[fngfs+ShellPatch::dsigmadxy],cg->fgfs[fngfs+ShellPatch::dsigmadxz],
                     cg->fgfs[fngfs+ShellPatch::dsigmadyy],cg->fgfs[fngfs+ShellPatch::dsigmadyz],cg->fgfs[fngfs+ShellPatch::dsigmadzz],
                     cg->fgfs[fngfs+ShellPatch::dRdxx],cg->fgfs[fngfs+ShellPatch::dRdxy],cg->fgfs[fngfs+ShellPatch::dRdxz],
                     cg->fgfs[fngfs+ShellPatch::dRdyy],cg->fgfs[fngfs+ShellPatch::dRdyz],cg->fgfs[fngfs+ShellPatch::dRdzz]);
      }
      if(BP==ble) break;
      BP=BP->next;
    }
}
void xm_patch::setupcordtrans()
{    
    MyList<Block> *BP=blb;
    while(BP)
    {
      Block *cg=BP->data;
      if(myrank == cg->rank) 
      {
        f_xm_getxyz(cg->shape,cg->X[0],cg->X[1],cg->X[2], 
                     cg->fgfs[fngfs+ShellPatch::gx],cg->fgfs[fngfs+ShellPatch::gy],cg->fgfs[fngfs+ShellPatch::gz]);
        f_xpm_getjacobian(cg->shape,cg->X[0],cg->X[1],cg->X[2], 
                     cg->fgfs[fngfs+ShellPatch::gx],cg->fgfs[fngfs+ShellPatch::gy],cg->fgfs[fngfs+ShellPatch::gz],
                     cg->fgfs[fngfs+ShellPatch::drhodx],cg->fgfs[fngfs+ShellPatch::drhody],cg->fgfs[fngfs+ShellPatch::drhodz],
                     cg->fgfs[fngfs+ShellPatch::dsigmadx],cg->fgfs[fngfs+ShellPatch::dsigmady],cg->fgfs[fngfs+ShellPatch::dsigmadz],
                     cg->fgfs[fngfs+ShellPatch::dRdx],cg->fgfs[fngfs+ShellPatch::dRdy],cg->fgfs[fngfs+ShellPatch::dRdz],
                     cg->fgfs[fngfs+ShellPatch::drhodxx],cg->fgfs[fngfs+ShellPatch::drhodxy],cg->fgfs[fngfs+ShellPatch::drhodxz],
                     cg->fgfs[fngfs+ShellPatch::drhodyy],cg->fgfs[fngfs+ShellPatch::drhodyz],cg->fgfs[fngfs+ShellPatch::drhodzz],
                     cg->fgfs[fngfs+ShellPatch::dsigmadxx],cg->fgfs[fngfs+ShellPatch::dsigmadxy],cg->fgfs[fngfs+ShellPatch::dsigmadxz],
                     cg->fgfs[fngfs+ShellPatch::dsigmadyy],cg->fgfs[fngfs+ShellPatch::dsigmadyz],cg->fgfs[fngfs+ShellPatch::dsigmadzz],
                     cg->fgfs[fngfs+ShellPatch::dRdxx],cg->fgfs[fngfs+ShellPatch::dRdxy],cg->fgfs[fngfs+ShellPatch::dRdxz],
                     cg->fgfs[fngfs+ShellPatch::dRdyy],cg->fgfs[fngfs+ShellPatch::dRdyz],cg->fgfs[fngfs+ShellPatch::dRdzz]);
      }
      if(BP==ble) break;
      BP=BP->next;
    }
}
void yp_patch::setupcordtrans()
{    
    MyList<Block> *BP=blb;
    while(BP)
    {
      Block *cg=BP->data;
      if(myrank == cg->rank) 
      {
        f_yp_getxyz(cg->shape,cg->X[0],cg->X[1],cg->X[2], 
                     cg->fgfs[fngfs+ShellPatch::gx],cg->fgfs[fngfs+ShellPatch::gy],cg->fgfs[fngfs+ShellPatch::gz]);
        f_ypm_getjacobian(cg->shape,cg->X[0],cg->X[1],cg->X[2], 
                     cg->fgfs[fngfs+ShellPatch::gx],cg->fgfs[fngfs+ShellPatch::gy],cg->fgfs[fngfs+ShellPatch::gz],
                     cg->fgfs[fngfs+ShellPatch::drhodx],cg->fgfs[fngfs+ShellPatch::drhody],cg->fgfs[fngfs+ShellPatch::drhodz],
                     cg->fgfs[fngfs+ShellPatch::dsigmadx],cg->fgfs[fngfs+ShellPatch::dsigmady],cg->fgfs[fngfs+ShellPatch::dsigmadz],
                     cg->fgfs[fngfs+ShellPatch::dRdx],cg->fgfs[fngfs+ShellPatch::dRdy],cg->fgfs[fngfs+ShellPatch::dRdz],
                     cg->fgfs[fngfs+ShellPatch::drhodxx],cg->fgfs[fngfs+ShellPatch::drhodxy],cg->fgfs[fngfs+ShellPatch::drhodxz],
                     cg->fgfs[fngfs+ShellPatch::drhodyy],cg->fgfs[fngfs+ShellPatch::drhodyz],cg->fgfs[fngfs+ShellPatch::drhodzz],
                     cg->fgfs[fngfs+ShellPatch::dsigmadxx],cg->fgfs[fngfs+ShellPatch::dsigmadxy],cg->fgfs[fngfs+ShellPatch::dsigmadxz],
                     cg->fgfs[fngfs+ShellPatch::dsigmadyy],cg->fgfs[fngfs+ShellPatch::dsigmadyz],cg->fgfs[fngfs+ShellPatch::dsigmadzz],
                     cg->fgfs[fngfs+ShellPatch::dRdxx],cg->fgfs[fngfs+ShellPatch::dRdxy],cg->fgfs[fngfs+ShellPatch::dRdxz],
                     cg->fgfs[fngfs+ShellPatch::dRdyy],cg->fgfs[fngfs+ShellPatch::dRdyz],cg->fgfs[fngfs+ShellPatch::dRdzz]);
      }
      if(BP==ble) break;
      BP=BP->next;
    }
}
void ym_patch::setupcordtrans()
{    
    MyList<Block> *BP=blb;
    while(BP)
    {
      Block *cg=BP->data;
      if(myrank == cg->rank) 
      {
        f_ym_getxyz(cg->shape,cg->X[0],cg->X[1],cg->X[2], 
                     cg->fgfs[fngfs+ShellPatch::gx],cg->fgfs[fngfs+ShellPatch::gy],cg->fgfs[fngfs+ShellPatch::gz]);
        f_ypm_getjacobian(cg->shape,cg->X[0],cg->X[1],cg->X[2], 
                     cg->fgfs[fngfs+ShellPatch::gx],cg->fgfs[fngfs+ShellPatch::gy],cg->fgfs[fngfs+ShellPatch::gz],
                     cg->fgfs[fngfs+ShellPatch::drhodx],cg->fgfs[fngfs+ShellPatch::drhody],cg->fgfs[fngfs+ShellPatch::drhodz],
                     cg->fgfs[fngfs+ShellPatch::dsigmadx],cg->fgfs[fngfs+ShellPatch::dsigmady],cg->fgfs[fngfs+ShellPatch::dsigmadz],
                     cg->fgfs[fngfs+ShellPatch::dRdx],cg->fgfs[fngfs+ShellPatch::dRdy],cg->fgfs[fngfs+ShellPatch::dRdz],
                     cg->fgfs[fngfs+ShellPatch::drhodxx],cg->fgfs[fngfs+ShellPatch::drhodxy],cg->fgfs[fngfs+ShellPatch::drhodxz],
                     cg->fgfs[fngfs+ShellPatch::drhodyy],cg->fgfs[fngfs+ShellPatch::drhodyz],cg->fgfs[fngfs+ShellPatch::drhodzz],
                     cg->fgfs[fngfs+ShellPatch::dsigmadxx],cg->fgfs[fngfs+ShellPatch::dsigmadxy],cg->fgfs[fngfs+ShellPatch::dsigmadxz],
                     cg->fgfs[fngfs+ShellPatch::dsigmadyy],cg->fgfs[fngfs+ShellPatch::dsigmadyz],cg->fgfs[fngfs+ShellPatch::dsigmadzz],
                     cg->fgfs[fngfs+ShellPatch::dRdxx],cg->fgfs[fngfs+ShellPatch::dRdxy],cg->fgfs[fngfs+ShellPatch::dRdxz],
                     cg->fgfs[fngfs+ShellPatch::dRdyy],cg->fgfs[fngfs+ShellPatch::dRdyz],cg->fgfs[fngfs+ShellPatch::dRdzz]);
      }
      if(BP==ble) break;
      BP=BP->next;
    }
}
void zp_patch::setupcordtrans()
{    
    MyList<Block> *BP=blb;
    while(BP)
    {
      Block *cg=BP->data;
      if(myrank == cg->rank) 
      {
        f_zp_getxyz(cg->shape,cg->X[0],cg->X[1],cg->X[2], 
                     cg->fgfs[fngfs+ShellPatch::gx],cg->fgfs[fngfs+ShellPatch::gy],cg->fgfs[fngfs+ShellPatch::gz]);
        f_zpm_getjacobian(cg->shape,cg->X[0],cg->X[1],cg->X[2], 
                     cg->fgfs[fngfs+ShellPatch::gx],cg->fgfs[fngfs+ShellPatch::gy],cg->fgfs[fngfs+ShellPatch::gz],
                     cg->fgfs[fngfs+ShellPatch::drhodx],cg->fgfs[fngfs+ShellPatch::drhody],cg->fgfs[fngfs+ShellPatch::drhodz],
                     cg->fgfs[fngfs+ShellPatch::dsigmadx],cg->fgfs[fngfs+ShellPatch::dsigmady],cg->fgfs[fngfs+ShellPatch::dsigmadz],
                     cg->fgfs[fngfs+ShellPatch::dRdx],cg->fgfs[fngfs+ShellPatch::dRdy],cg->fgfs[fngfs+ShellPatch::dRdz],
                     cg->fgfs[fngfs+ShellPatch::drhodxx],cg->fgfs[fngfs+ShellPatch::drhodxy],cg->fgfs[fngfs+ShellPatch::drhodxz],
                     cg->fgfs[fngfs+ShellPatch::drhodyy],cg->fgfs[fngfs+ShellPatch::drhodyz],cg->fgfs[fngfs+ShellPatch::drhodzz],
                     cg->fgfs[fngfs+ShellPatch::dsigmadxx],cg->fgfs[fngfs+ShellPatch::dsigmadxy],cg->fgfs[fngfs+ShellPatch::dsigmadxz],
                     cg->fgfs[fngfs+ShellPatch::dsigmadyy],cg->fgfs[fngfs+ShellPatch::dsigmadyz],cg->fgfs[fngfs+ShellPatch::dsigmadzz],
                     cg->fgfs[fngfs+ShellPatch::dRdxx],cg->fgfs[fngfs+ShellPatch::dRdxy],cg->fgfs[fngfs+ShellPatch::dRdxz],
                     cg->fgfs[fngfs+ShellPatch::dRdyy],cg->fgfs[fngfs+ShellPatch::dRdyz],cg->fgfs[fngfs+ShellPatch::dRdzz]);
      }
      if(BP==ble) break;
      BP=BP->next;
    }
}
void zm_patch::setupcordtrans()
{    
    MyList<Block> *BP=blb;
    while(BP)
    {
      Block *cg=BP->data;
      if(myrank == cg->rank) 
      {
        f_zm_getxyz(cg->shape,cg->X[0],cg->X[1],cg->X[2], 
                     cg->fgfs[fngfs+ShellPatch::gx],cg->fgfs[fngfs+ShellPatch::gy],cg->fgfs[fngfs+ShellPatch::gz]);
        f_zpm_getjacobian(cg->shape,cg->X[0],cg->X[1],cg->X[2], 
                     cg->fgfs[fngfs+ShellPatch::gx],cg->fgfs[fngfs+ShellPatch::gy],cg->fgfs[fngfs+ShellPatch::gz],
                     cg->fgfs[fngfs+ShellPatch::drhodx],cg->fgfs[fngfs+ShellPatch::drhody],cg->fgfs[fngfs+ShellPatch::drhodz],
                     cg->fgfs[fngfs+ShellPatch::dsigmadx],cg->fgfs[fngfs+ShellPatch::dsigmady],cg->fgfs[fngfs+ShellPatch::dsigmadz],
                     cg->fgfs[fngfs+ShellPatch::dRdx],cg->fgfs[fngfs+ShellPatch::dRdy],cg->fgfs[fngfs+ShellPatch::dRdz],
                     cg->fgfs[fngfs+ShellPatch::drhodxx],cg->fgfs[fngfs+ShellPatch::drhodxy],cg->fgfs[fngfs+ShellPatch::drhodxz],
                     cg->fgfs[fngfs+ShellPatch::drhodyy],cg->fgfs[fngfs+ShellPatch::drhodyz],cg->fgfs[fngfs+ShellPatch::drhodzz],
                     cg->fgfs[fngfs+ShellPatch::dsigmadxx],cg->fgfs[fngfs+ShellPatch::dsigmadxy],cg->fgfs[fngfs+ShellPatch::dsigmadxz],
                     cg->fgfs[fngfs+ShellPatch::dsigmadyy],cg->fgfs[fngfs+ShellPatch::dsigmadyz],cg->fgfs[fngfs+ShellPatch::dsigmadzz],
                     cg->fgfs[fngfs+ShellPatch::dRdxx],cg->fgfs[fngfs+ShellPatch::dRdxy],cg->fgfs[fngfs+ShellPatch::dRdxz],
                     cg->fgfs[fngfs+ShellPatch::dRdyy],cg->fgfs[fngfs+ShellPatch::dRdyz],cg->fgfs[fngfs+ShellPatch::dRdzz]);
      }
      if(BP==ble) break;
      BP=BP->next;
    }
}
ShellPatch::ShellPatch(int ingfsi, int fngfsi,int *shapei, double *Rrangei,int Symmetry,int myranki):
ingfs(ingfsi),fngfs(fngfsi),myrank(myranki),PatL(0)
{
  for(int i=0;i<dim;i++)
  {
     shape[i] = shapei[i];
// we always assume the input parameter is in cell center style     
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
     shape[i] = shape[i]+1;
#endif
  }
// change from cardisian r to local corrdinate r  
  Rrange[0] = getR(Rrangei[0]);
  Rrange[1] = getR(Rrangei[1]);

  if(myrank==0)
  {
    cout<<"shell's range: ["<<Rrange[0]<<":"<<Rrange[1]<<"]"<<endl
	<<"shape: "<<shape[2]<<endl
	<<"resolution: ["<<getdX(0)<<","<<getdX(1)<<","<<getdX(2)<<"]"<<endl;
  }
// extend buffer points for lower boundary
  Rrange[0] -= buffer_width*getdX(2);
  shape[2] += buffer_width;

  double bbox[2*dim];
  int shape_here[dim];
  bbox[2]=Rrange[0];
  bbox[5]=Rrange[1];
  shape_here[2] = shape[2];

  switch(Symmetry)
  {
     case 0:
        for(int i=0;i<2;i++) shape_here[i] = shape[i]+2*overghost;
	bbox[0]=-PI/4-overghost*getdX(0);
	bbox[1]=-PI/4-overghost*getdX(1);
	bbox[3]= PI/4+overghost*getdX(0);
	bbox[4]= PI/4+overghost*getdX(1);
	PatL = new MyList<ss_patch>;
	PatL->data= new xp_patch(ingfs,fngfs,shape_here,bbox,myrank);
	PatL->insert(new xm_patch(ingfs,fngfs,shape_here,bbox,myrank));
	PatL->insert(new yp_patch(ingfs,fngfs,shape_here,bbox,myrank));
	PatL->insert(new ym_patch(ingfs,fngfs,shape_here,bbox,myrank));
	PatL->insert(new zp_patch(ingfs,fngfs,shape_here,bbox,myrank));
	PatL->insert(new zm_patch(ingfs,fngfs,shape_here,bbox,myrank));
        break;
     case 1:
        for(int i=0;i<2;i++) shape_here[i] = shape[i]+2*overghost;
	bbox[0]=-PI/4-overghost*getdX(0);
	bbox[1]=-PI/4-overghost*getdX(1);
	bbox[3]= PI/4+overghost*getdX(0);
	bbox[4]= PI/4+overghost*getdX(1);
	PatL = new MyList<ss_patch>;
	PatL->data= new zp_patch(ingfs,fngfs,shape_here,bbox,myrank);
        shape_here[0] = shape[0]+2*overghost;
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
        shape_here[1] = (shape[1]+1)/2+overghost;
#else
#ifdef Cell
        shape_here[1] = shape[1]/2+overghost;
#else
#error Not define Vertex nor Cell
#endif  
#endif
	bbox[0]=-PI/4-overghost*getdX(0);
	bbox[1]=0;
	bbox[3]= PI/4+overghost*getdX(0);
	bbox[4]= PI/4+overghost*getdX(1);
	PatL->insert(new xp_patch(ingfs,fngfs,shape_here,bbox,myrank));
	PatL->insert(new yp_patch(ingfs,fngfs,shape_here,bbox,myrank));
	bbox[0]=-PI/4-overghost*getdX(0);
	bbox[1]=-PI/4-overghost*getdX(1);
	bbox[3]= PI/4+overghost*getdX(0);
	bbox[4]= 0;
	PatL->insert(new xm_patch(ingfs,fngfs,shape_here,bbox,myrank));
	PatL->insert(new ym_patch(ingfs,fngfs,shape_here,bbox,myrank));
	break;
     case 2:
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
        for(int i=0;i<2;i++) shape_here[i] = (shape[i]+1)/2+overghost;
#else
#ifdef Cell
        for(int i=0;i<2;i++) shape_here[i] = shape[i]/2+overghost;
#else
#error Not define Vertex nor Cell
#endif  
#endif
	bbox[0]=0;
	bbox[1]=0;
	bbox[3]=PI/4+overghost*getdX(0);
	bbox[4]=PI/4+overghost*getdX(1);
	PatL = new MyList<ss_patch>;
	PatL->data= new zp_patch(ingfs,fngfs,shape_here,bbox,myrank);
	PatL->insert(new xp_patch(ingfs,fngfs,shape_here,bbox,myrank));
	PatL->insert(new yp_patch(ingfs,fngfs,shape_here,bbox,myrank));
	break;
     default:
	cout<<"not recognized Symmetry type"<<endl;
        MPI_Abort(MPI_COMM_WORLD,1);
  }
}
ShellPatch::~ShellPatch()
{
    int nprocs = 1;
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    for(int node=0;node<nprocs;node++)
    {
      if(ss_src[node]) destroypsuList(ss_src[node]);
      if(ss_dst[node]) destroypsuList(ss_dst[node]);
      if(cs_src[node]) destroypsuList(cs_src[node]);
      if(cs_dst[node]) destroypsuList(cs_dst[node]);
    }

    delete[] ss_src; delete[] ss_dst;
    delete[] cs_src; delete[] cs_dst;

    while(PatL)
    {
      ss_patch *sPp=PatL->data;
      MyList<Block> *bg;
      while(sPp->blb)
      {
        if(sPp->blb==sPp->ble) break;
        bg=(sPp->blb->next)?sPp->blb->next:0;
        delete sPp->blb->data;
        delete sPp->blb;
        sPp->blb=bg;
      }
      if(sPp->ble)
      {
        delete sPp->ble->data;
        delete sPp->ble;
      }
      sPp->blb=sPp->ble=0;
      PatL=PatL->next;
   }
   PatL->destroyList();
}
void ShellPatch::destroypsuList(MyList<pointstru> *ct)
{
  MyList<pointstru> *n;
  while(ct)
  {
    n = ct->next;
    if(ct->data->coef) {delete[] ct->data->coef; delete[] ct->data->sind;}
    delete ct->data;
    delete ct;
    ct = n;
  }
}
double ShellPatch::getR(double r)
{
  double A=1,B=0,r0=0,eps=1;
  f_shellcordpar(A,B,r0,eps);
  double f = A*(r-r0)+B*sqrt(1+(r-r0)*(r-r0)/eps);
  return f+A*r0-B*sqrt(1+r0*r0/eps);
}
double ShellPatch::getsr(double R)
{
  double A=1,B=0,r0=0,eps=1;
  f_shellcordpar(A,B,r0,eps);
  double f = R+B;
  return r0+(A*f-B*sqrt(A*A+(f*f-B*B)/eps))/(A*A-B*B/eps);
}
MyList<Block> *ShellPatch::compose_sh(int cpusize)
{
  if(dim != 3)
  {
    cout<<"distrivute: now we only support 3-dimension"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

//  checkPatch();

  bool periodic=false;
  MyList<Block> *BlL=0;

  int split_size,min_size,block_size=0;
     
  int min_width=2*max(ghost_width,buffer_width);
  int nxyz[dim],mmin_width[dim],min_shape[dim];  
  
  MyList<ss_patch> *PLi = PatL;
  for(int i=0;i<dim;i++) min_shape[i] = PLi->data->shape[i];
  PLi = PLi->next;
  while(PLi)
  {
    ss_patch *PP = PLi->data;
    for(int i=0;i<dim;i++) min_shape[i] = min(min_shape[i],PP->shape[i]);
    PLi = PLi->next;
  }  

  for(int i=0;i<dim;i++) mmin_width[i] = min(min_width,min_shape[i]);

  min_size = mmin_width[0];
  for(int i=1;i<dim;i++) min_size = min_size*mmin_width[i];
    
  PLi = PatL;
  while(PLi)
  {
    ss_patch *PP = PLi->data;
//    PP->checkPatch(true);
    int bs=PP->shape[0];
    for(int i=1;i<dim;i++) bs = bs*PP->shape[i];
    block_size=block_size + bs;
    PLi = PLi->next;
  }  
  split_size=max(min_size,block_size/cpusize);
  split_size=max(1,split_size);

  int n_rank=0;
  PLi = PatL;
  int reacpu = 0;
  while(PLi)
  {
    ss_patch *PP = PLi->data;

    reacpu += Parallel::partition3(nxyz, split_size,mmin_width, cpusize,PP->shape);

    Block *ng;
    int shape_here[dim],ibbox_here[2*dim];
    double bbox_here[2*dim],dd;

//ibbox : 0,...N-1
      for(int i=0;i<nxyz[0];i++)
	 for(int j=0;j<nxyz[1];j++)
	   for(int k=0;k<nxyz[2];k++)
     	   {
	    ibbox_here[0]=(PP->shape[0]*i)/nxyz[0];
    	    ibbox_here[3]=(PP->shape[0]*(i+1))/nxyz[0]-1;
    	    ibbox_here[1]=(PP->shape[1]*j)/nxyz[1];
    	    ibbox_here[4]=(PP->shape[1]*(j+1))/nxyz[1]-1;
    	    ibbox_here[2]=(PP->shape[2]*k)/nxyz[2];
    	    ibbox_here[5]=(PP->shape[2]*(k+1))/nxyz[2]-1;

	    if(periodic)
	    {
    	     ibbox_here[0]=ibbox_here[0]-ghost_width;
    	     ibbox_here[3]=ibbox_here[3]+ghost_width;
    	     ibbox_here[1]=ibbox_here[1]-ghost_width;
    	     ibbox_here[4]=ibbox_here[4]+ghost_width;	    
    	     ibbox_here[2]=ibbox_here[2]-ghost_width;
    	     ibbox_here[5]=ibbox_here[5]+ghost_width;
	    }
	    else
	    {
    	     ibbox_here[0]=max(0,ibbox_here[0]-ghost_width);
    	     ibbox_here[3]=min(PP->shape[0]-1,ibbox_here[3]+ghost_width);
    	     ibbox_here[1]=max(0,ibbox_here[1]-ghost_width);
    	     ibbox_here[4]=min(PP->shape[1]-1,ibbox_here[4]+ghost_width);	    
    	     ibbox_here[2]=max(0,ibbox_here[2]-ghost_width);
    	     ibbox_here[5]=min(PP->shape[2]-1,ibbox_here[5]+ghost_width);
	    }
	    
    	    shape_here[0]=ibbox_here[3]-ibbox_here[0]+1;
    	    shape_here[1]=ibbox_here[4]-ibbox_here[1]+1;
    	    shape_here[2]=ibbox_here[5]-ibbox_here[2]+1;
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
            dd=(PP->bbox[3]-PP->bbox[0])/(PP->shape[0]-1);
    	    bbox_here[0]=PP->bbox[0]+ibbox_here[0]*dd;
    	    bbox_here[3]=PP->bbox[0]+ibbox_here[3]*dd;
	    
            dd=(PP->bbox[4]-PP->bbox[1])/(PP->shape[1]-1);
    	    bbox_here[1]=PP->bbox[1]+ibbox_here[1]*dd;
    	    bbox_here[4]=PP->bbox[1]+ibbox_here[4]*dd;

            dd=(PP->bbox[5]-PP->bbox[2])/(PP->shape[2]-1);
    	    bbox_here[2]=PP->bbox[2]+ibbox_here[2]*dd;
    	    bbox_here[5]=PP->bbox[2]+ibbox_here[5]*dd;
#else
#ifdef Cell
            dd=(PP->bbox[3]-PP->bbox[0])/PP->shape[0];
    	    bbox_here[0]=PP->bbox[0]+(ibbox_here[0]  )*dd;
    	    bbox_here[3]=PP->bbox[0]+(ibbox_here[3]+1)*dd;
	    
            dd=(PP->bbox[4]-PP->bbox[1])/PP->shape[1];
    	    bbox_here[1]=PP->bbox[1]+(ibbox_here[1]  )*dd;
    	    bbox_here[4]=PP->bbox[1]+(ibbox_here[4]+1)*dd;

            dd=(PP->bbox[5]-PP->bbox[2])/PP->shape[2];
    	    bbox_here[2]=PP->bbox[2]+(ibbox_here[2]  )*dd;
    	    bbox_here[5]=PP->bbox[2]+(ibbox_here[5]+1)*dd;
#else
#error Not define Vertex nor Cell
#endif  
#endif 

	    ng=new Block(dim,shape_here,bbox_here,n_rank++,ingfs,fngfs+dRdzz+1,0);  //delete through KillBlocks
//	    ng->checkBlock();
    	    if(n_rank==cpusize)n_rank=0;
    	    if(BlL) BlL->insert(ng); 
	    else    BlL = new MyList<Block>(ng);  //delete through KillBlocks

// set PP->blb      
	    if(i == 0 && j == 0 && k == 0) 
	    {
	       MyList<Block> *Bp=BlL;
	       while(Bp->data != ng) Bp=Bp->next;
	       PP->blb = Bp;
	    }
     	   }
// set PP->ble      
    {
      MyList<Block> *Bp=BlL;
      while(Bp->data != ng) Bp=Bp->next;
      PP->ble = Bp;
    }
    PLi = PLi->next;
  }
  if(reacpu < cpusize*2/3)
  {
     int myrank;
     MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
     if(myrank == 0) cout<<"ShellPatch::distribute CAUSTION: uses essencially "<<reacpu<<" processors vs "<<
                       cpusize<<" cpus run, your scientific computation scale is not as large as you estimate."<<endl;
  }

  return BlL;
}
void ShellPatch::getlocalpox(double x,double y,double z,int &sst,double &lx,double &ly,double &lz)
{
     double r;
     r = sqrt(x*x+y*y+z*z);
     lz = getR(r);
     if(fabs(x) <= z && fabs(y)<= z)
     {
	sst=0;
        lx=atan(x/z);
	ly=atan(y/z);
     }
     else if(fabs(x) <= -z && fabs(y)<= -z)
     {
	sst=1;
        lx=atan(x/z);
	ly=atan(y/z);
     }
     else if(fabs(y) <=  x && fabs(z)<=  x)
     {
	sst=2;
        lx=atan(y/x);
        ly=atan(z/x);
     }
     else if(fabs(y) <= -x && fabs(z)<= -x)
     {
	sst=3;
        lx=atan(y/x);
        ly=atan(z/x);
     }
     else if(fabs(x) <=  y && fabs(z)<=  y)
     {
	sst=4;
        lx=atan(x/y);
        ly=atan(z/y);
     }
     else if(fabs(x) <= -y && fabs(z)<= -y)
     {
	sst=5;
        lx=atan(x/y);
        ly=atan(z/y);
     }
     else
     {
	cout<<"ShellPatch::getlocalpox should not come here, something wrong"<<endl;
     }
}
void ShellPatch::getlocalpoxsst(double x,double y,double z,int sst,double &lx,double &ly,double &lz)
{
     double r;
     r = sqrt(x*x+y*y+z*z);
     lz = getR(r);
     switch(sst)
     {
	case -1:
            lx = x;
	    ly = y;
	    lz = z;
	    break;
	case 0:
            lx=atan(x/z);
	    ly=atan(y/z);
            break;
	case 1:
            lx=atan(x/z);
	    ly=atan(y/z);
            break;
	case 2:
            lx=atan(y/x);
	    ly=atan(z/x);
	    break;
	case 3:
            lx=atan(y/x);
	    ly=atan(z/x);
            break;
	case 4:
            lx=atan(x/y);
	    ly=atan(z/y);
            break;
	case 5:
            lx=atan(x/y);
	    ly=atan(z/y);
            break;
	default:
  	    cout<<"ShellPatch::getlocalpoxsst should not come here, something wrong"<<endl;
     }
}
void ShellPatch::getglobalpox(double &x,double &y,double &z,int sst,double lx,double ly,double lz)
{
     double r=getsr(lz);
     switch(sst)
     {
	case 0:
	    x=tan(lx);
	    y=tan(ly);
	    z=r/sqrt(1+x*x+y*y);
	    x=z*x;
	    y=z*y;
            break;
	case 1:
	    x=tan(lx);
	    y=tan(ly);
	    z=-r/sqrt(1+x*x+y*y);
	    x=z*x;
	    y=z*y;
            break;
	case 2:
	    y=tan(lx);
	    z=tan(ly);
	    x=r/sqrt(1+z*z+y*y);
	    y=x*y;
	    z=x*z;
            break;
	case 3:
	    y=tan(lx);
	    z=tan(ly);
	    x=-r/sqrt(1+z*z+y*y);
	    y=x*y;
	    z=x*z;
            break;
	case 4:
	    x=tan(lx);
	    z=tan(ly);
	    y=r/sqrt(1+x*x+z*z);
	    x=y*x;
	    z=y*z;
            break;
	case 5:
	    x=tan(lx);
	    z=tan(ly);
	    y=-r/sqrt(1+x*x+z*z);
	    x=y*x;
	    z=y*z;
            break;
     }
}
int ShellPatch::getdumydimension(int acsst,int posst) // -1 means no dumy dimension
{
   int dms;
   if(acsst==-1 || posst==-1) return -1;
   switch(acsst)
   {
      case 0:
      case 1:
	  switch(posst)
	  {
             case 0:
	     case 1:
		 cout<<"error in ShellPatch::getdumydimension: acsst = "<<acsst<<", posst = "<<posst<<endl;
		 return -1;
	     case 2:
	     case 3:
		 return 0;
	     case 4:
	     case 5:
		 return 1;
	     default:
		 cout<<"error in ShellPatch::getdumydimension: posst = "<<posst<<endl;
		 return -1;
	  }
      case 2:
      case 3:
	  switch(posst)
	  {
             case 0:
	     case 1:
		 return 1;
	     case 2:
	     case 3:
		 cout<<"error in ShellPatch::getdumydimension: acsst = "<<acsst<<", posst = "<<posst<<endl;
		 return -1;
	     case 4:
	     case 5:
		 return 0;
	     default:
		 cout<<"error in ShellPatch::getdumydimension: posst = "<<posst<<endl;
		 return -1;
	  }
      case 4:
      case 5:
	  switch(posst)
	  {
             case 0:
	     case 1:
		 return 1;
	     case 2:
	     case 3:
		 return 0;
	     case 4:
	     case 5:
		 cout<<"error in ShellPatch::getdumydimension: acsst = "<<acsst<<", posst = "<<posst<<endl;
		 return -1;
	     default:
		 cout<<"error in ShellPatch::getdumydimension: posst = "<<posst<<endl;
		 return -1;
	  }
      default:
	  cout<<"error in ShellPatch::getdumydimension: acsst = "<<acsst<<endl;
	  return -1;
   }
}
// used by _dst construction, so these x,y,z must coinside with grid point
// we have considered ghost points now
void ShellPatch::prolongpointstru(MyList<pointstru>*&psul,MyList<ss_patch> *sPpi,double DH[dim],
		                  MyList<Patch> *Ppi,double CDH[dim],MyList<pointstru>*pss)
{
   int n_dst=0;
   MyList<ss_patch> *sPp=sPpi;
   MyList<Patch> *Pp=Ppi;
   MyList<Block> *Bgl;
   Block *Bg;
   double llb[dim],uub[dim];
   double lx,ly,lz;

   if(pss->data->tsst >= 0)
   {
     getlocalpoxsst(pss->data->gpox[0],pss->data->gpox[1],pss->data->gpox[2],pss->data->tsst,
                    lx,ly,lz);
     while(sPp)
     {
        if(sPp->data->sst == pss->data->tsst)
        {
          Bgl = sPp->data->blb;
          while(Bgl)
          {
	     Bg = Bgl->data;
	     {
		for(int j=0;j<dim;j++)
		{
		   llb[j] = Bg->bbox[j];
		   uub[j] = Bg->bbox[j+dim];
		}
		
		if(lx>llb[0]-0.1*DH[0] && lx<uub[0]+0.1*DH[0] &&
                   ly>llb[1]-0.1*DH[1] && ly<uub[1]+0.1*DH[1] &&
                   lz>llb[2]-0.1*DH[2] && lz<uub[2]+0.1*DH[2])
		{
		  MyList<pointstru>* ps=new MyList<pointstru>;
	       	  ps->data=new pointstru;
	       	  ps->next=0;
	       	  for(int i=0;i<dim;i++) ps->data->gpox[i]=pss->data->gpox[i];
		  ps->data->lpox[0]=lx;
		  ps->data->lpox[1]=ly;
		  ps->data->lpox[2]=lz;
	       	  ps->data->ssst=pss->data->ssst;
	       	  ps->data->tsst=sPp->data->sst;
		  ps->data->dumyd=getdumydimension(ps->data->tsst,ps->data->ssst);
		  ps->data->Bg=Bg;
		  ps->data->coef=0;
		  ps->data->sind=0;
	       	  if(psul) psul->catList(ps);
	       	  else     psul         =ps;
	          n_dst++;
		}
	     }
             if(Bgl==sPp->data->ble) break;
	     Bgl=Bgl->next;
	  }
	}
   	sPp = sPp->next;
     }
   }
   else
   {
     if(pss->data->tsst != -1) cout<<"somthing is wrong in ShellPatch::prolongpointstru"<<endl;
     lx=pss->data->gpox[0];
     ly=pss->data->gpox[1];
     lz=pss->data->gpox[2];
     while(Pp)
     {
       Bgl = Pp->data->blb;
       while(Bgl)
       {
         Bg = Bgl->data;
         {
	   for(int j=0;j<dim;j++)
	   {
	      llb[j] = Bg->bbox[j];
	      uub[j] = Bg->bbox[j+dim];
	   }
	   
	   if(lx>llb[0]-0.1*CDH[0] && lx<uub[0]+0.1*CDH[0] &&
              ly>llb[1]-0.1*CDH[1] && ly<uub[1]+0.1*CDH[1] &&
              lz>llb[2]-0.1*CDH[2] && lz<uub[2]+0.1*CDH[2])
	   {
	     MyList<pointstru>* ps=new MyList<pointstru>;
	     ps->data=new pointstru;
	     ps->next=0;
	     for(int i=0;i<dim;i++) ps->data->gpox[i]=pss->data->gpox[i];
	     ps->data->lpox[0]=lx;
	     ps->data->lpox[1]=ly;
	     ps->data->lpox[2]=lz;
	     ps->data->ssst=pss->data->ssst;
	     ps->data->tsst=-1;
             ps->data->dumyd=getdumydimension(ps->data->tsst,ps->data->ssst);
	     ps->data->Bg=Bg;
	     ps->data->coef=0;
	     ps->data->sind=0;
	     if(psul) psul->catList(ps);
	     else     psul         =ps;
	     n_dst++;
	   }
	 }
         if(Bgl==Pp->data->ble) break;
	 Bgl=Bgl->next;
       }
       Pp = Pp->next;
     }
   }
//if n_dst > 0, that's because of ghost_points
   if(n_dst == 0)
   {
     int myrank = 0;
     MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
     if(myrank==0) cout<<"ShellPatch::prolongpointstru fail to find target Block for pointstru:"<<endl;
     check_pointstrul(pss,true);
     if(Pp == Ppi)
     {
        getlocalpoxsst(pss->data->gpox[0],pss->data->gpox[1],pss->data->gpox[2],pss->data->tsst,
                       lx,ly,lz);
	if(myrank==0) cout<<"sst = "<<pss->data->tsst<<", lx,ly,lz = "<<lx<<","<<ly<<","<<lz<<endl;
        checkBlock(pss->data->tsst);
     }
     else
     {
	Pp=Ppi;
	while(Pp)
	{
	   Pp->data->checkBlock();
	   Pp=Pp->next;
	}
     }
     if(myrank==0) MPI_Abort(MPI_COMM_WORLD,1);
   }
   else
   {
      MyList<pointstru>* ts=0;
      for(int i=1;i<n_dst;i++)
      {
	 MyList<pointstru>* ps=new MyList<pointstru>;
	 ps->data=new pointstru;
	 ps->next=(i==n_dst-1)?pss->next:0;
	 for(int i=0;i<dim;i++) {ps->data->gpox[i]=pss->data->gpox[i]; ps->data->lpox[i]=pss->data->lpox[i];}
	 ps->data->ssst=pss->data->ssst;
	 ps->data->tsst=pss->data->tsst;
	 ps->data->dumyd=getdumydimension(ps->data->ssst,ps->data->tsst);
	 ps->data->Bg=pss->data->Bg;
         ps->data->coef=0;
         ps->data->sind=0;
         if(ts) ts->catList(ps);
	 else   ts=         ps;
      }
      if(ts) pss->next=ts;
   }
}
// used by _src construction, so these x,y,z do not coinside with grid point	
bool ShellPatch::prolongpointstru(MyList<pointstru>*&psul,bool ssyn,int tsst, MyList<ss_patch> *sPp, double DH[dim],
                                  MyList<Patch> *Pp,double CDH[dim],double x,double y,double z,int Symmetry,int rank_in)
{
   MyList<Block> *Bgl;
   Block *Bg;
   double llb[dim],uub[dim];
   double lx,ly,lz;

   if(ssyn)
   {
     int sst;
     getlocalpox(x,y,z,sst,lx,ly,lz);
     while(sPp)
     {
        if(sPp->data->sst == sst)
        {
          Bgl = sPp->data->blb;
          while(Bgl)
          {
	     Bg = Bgl->data;
	     if(Bg->rank == rank_in)
	     {
		for(int j=0;j<2;j++)
		{
                  if(feq(Bg->bbox[j],-PI/4-overghost*DH[j],DH[j]/2))   llb[j]=-PI/4;
	          else if(feq(Bg->bbox[j],sPp->data->bbox[j],DH[j]/2)) llb[j]=Bg->bbox[j];
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
		  else                                                 llb[j]=Bg->bbox[j]+(ghost_width-1)*DH[j];
#else
#ifdef Cell
		  else                                                 llb[j]=Bg->bbox[j]+ghost_width*DH[j];
#else
#error Not define Vertex nor Cell
#endif  
#endif 	  
                  if(feq(Bg->bbox[dim+j],PI/4+overghost*DH[j],DH[j]/2))        uub[j]=PI/4;
	          else if(feq(Bg->bbox[dim+j],sPp->data->bbox[dim+j],DH[j]/2)) uub[j]=Bg->bbox[dim+j];
		  else                                                         uub[j]=Bg->bbox[dim+j]-ghost_width*DH[j];
		}
                if(feq(Bg->bbox[2],sPp->data->bbox[2],DH[2]/2)) llb[2]=Bg->bbox[2];
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
		else                                            llb[2]=Bg->bbox[2]+(ghost_width-1)*DH[2];
#else
#ifdef Cell
		else                                            llb[2]=Bg->bbox[2]+ghost_width*DH[2];
#else
#error Not define Vertex nor Cell
#endif  
#endif 	 
                if(feq(Bg->bbox[dim+2],sPp->data->bbox[dim+2],DH[2]/2)) uub[2]=Bg->bbox[dim+2];
		else                                                    uub[2]=Bg->bbox[dim+2]-ghost_width*DH[2];
		if(lx>llb[0]-0.0001*DH[0] && lx<uub[0]+0.0001*DH[0] &&
                   ly>llb[1]-0.0001*DH[1] && ly<uub[1]+0.0001*DH[1] &&
                   lz>llb[2]-0.0001*DH[2] && lz<uub[2]+0.0001*DH[2]) // even ghost_width-1 the region is like |----|----|
			                                             //                                            ^
							             // so for ^ point may miss for vertext center, so we use 0.0001
		{
		  MyList<pointstru>* ps=new MyList<pointstru>;
                  ps->data=new pointstru;
		  ps->data->Bg=Bg;
		  ps->data->gpox[0]=x;
		  ps->data->gpox[1]=y;
		  ps->data->gpox[2]=z;
		  ps->data->lpox[0]=lx;
		  ps->data->lpox[1]=ly;
		  ps->data->lpox[2]=lz;
		  ps->data->ssst = sPp->data->sst;
		  ps->data->tsst = tsst;
		  ps->data->dumyd=getdumydimension(ps->data->ssst,ps->data->tsst);
		  ps->data->coef=0;
		  ps->data->sind=0;
		  ps->next=0;
	       	  if(psul) psul->catList(ps);
	       	  else     psul         =ps;
		  return true;
		}
	     }
             if(Bgl==sPp->data->ble) break;
	     Bgl=Bgl->next;
	  }
	}
	sPp = sPp->next;
     }
   }
   else   
   {
     while(Pp)
     {
       Bgl = Pp->data->blb;
       while(Bgl)
       {
         Bg = Bgl->data;
         if(Bg->rank == rank_in)
         {
           for(int j=0;j<dim;j++)
	   {
             if(feq(Bg->bbox[j],Pp->data->bbox[j],CDH[j]/2)) llb[j]=Bg->bbox[j];
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
	     else                                            llb[j]=Bg->bbox[j]+(ghost_width-1)*CDH[j];
#else
#ifdef Cell
	     else                                            llb[j]=Bg->bbox[j]+ghost_width*CDH[j];
#else
#error Not define Vertex nor Cell
#endif  
#endif 	     
             if(feq(Bg->bbox[dim+j],Pp->data->bbox[dim+j],CDH[j]/2)) uub[j]=Bg->bbox[dim+j];
	     else                                                    uub[j]=Bg->bbox[dim+j]-ghost_width*CDH[j];
	   }
	   if(x>llb[0]-0.0001*CDH[0] && x<uub[0]+0.0001*CDH[0] &&
              y>llb[1]-0.0001*CDH[1] && y<uub[1]+0.0001*CDH[1] &&
              z>llb[2]-0.0001*CDH[2] && z<uub[2]+0.0001*CDH[2])
	   {
	      MyList<pointstru>* ps=new MyList<pointstru>;
              ps->data=new pointstru;
	      ps->data->Bg=Bg;
	      ps->data->gpox[0]=x;
	      ps->data->gpox[1]=y;
	      ps->data->gpox[2]=z;
	      ps->data->lpox[0]=x;
	      ps->data->lpox[1]=y;
	      ps->data->lpox[2]=z;
	      ps->data->ssst = -1;
	      ps->data->tsst = tsst;
              ps->data->dumyd=getdumydimension(ps->data->ssst,ps->data->tsst);
              ps->data->coef=0;
	      ps->data->sind=0;
	      ps->next=0;
	      if(psul) psul->catList(ps);
	      else     psul         =ps;
	      return true;
	   }
	 }
         if(Bgl==Pp->data->ble) break;
	 Bgl=Bgl->next;
       }
       Pp = Pp->next;
     }
   }

   return false;
}
#define CS_width buffer_width
//setup interpatch interpolation stuffs
void ShellPatch::setupintintstuff(int cpusize, MyList<Patch> *CPatL,int Symmetry)
{
  int myrank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  if(myrank==0) cout<<"ShellPatch::setupintintstuff begines..."<<endl;

  ss_src=new MyList<pointstru> *[cpusize];
  ss_dst=new MyList<pointstru> *[cpusize];  
  cs_src=new MyList<pointstru> *[cpusize];
  cs_dst=new MyList<pointstru> *[cpusize];

  MyList<pointstru> *ps,*ts;
  MyList<ss_patch> *sPp;
  MyList<Block> *Bgl;
  MyList<Patch> *Pp;
  Block *Bg;
  double CDH[dim],DH[dim],llb[dim],uub[dim];
  double x,y,z;

  for(int i=0;i<dim;i++)  {CDH[i] = CPatL->data->getdX(i); DH[i] = getdX(i);}

  for(int i=0;i<cpusize;i++)
  {
     ss_src[i]=0;
     cs_src[i]=0;
     ss_dst[i]=0;
     cs_dst[i]=0;
  }

     sPp=PatL;
     while(sPp)
     {
        for(int iz=0;iz<sPp->data->shape[2];iz++)
	  for(int is=0;is<sPp->data->shape[1];is++)
	    for(int ir=0;ir<sPp->data->shape[0];ir++)
	    {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
	       x=sPp->data->bbox[0]+ir*DH[0];
	       y=sPp->data->bbox[1]+is*DH[1];
	       z=sPp->data->bbox[2]+iz*DH[2];
#else
#ifdef Cell
	       x=sPp->data->bbox[0]+(ir+0.5)*DH[0];
	       y=sPp->data->bbox[1]+(is+0.5)*DH[1];
	       z=sPp->data->bbox[2]+(iz+0.5)*DH[2];
#else
#error Not define Vertex nor Cell
#endif  
#endif
	       if(z<sPp->data->bbox[2]+(CS_width+0.0001)*DH[2])
	       {
		  double gx,gy,gz;
		  getglobalpox(gx,gy,gz,sPp->data->sst,x,y,z);
		  bool flag=false;
                  for(int i=0;i<cpusize;i++) 
		  {
		     flag=prolongpointstru(cs_src[i],false,sPp->data->sst,PatL,DH,CPatL,CDH,gx,gy,gz,Symmetry,i);
		     if(flag) break;
		  }
		  if(!flag)
                  {
	            CPatL->data->checkBlock();
	       	    if(myrank==0)
	       	    {
		        cout<<"ShellPatch::prolongpointstru fail to find cardisian source point for"<<endl;
		 	cout<<"sst = "<<sPp->data->sst<<" lx,ly,lz = "<<x<<","<<y<<","<<z<<endl;
		 	cout<<"x,y,z = "<<gx<<","<<gy<<","<<gz<<endl;
                        MPI_Abort(MPI_COMM_WORLD,1);
	       	    }
	       	  }
	       }
//	       else if(x<-PI/4-(overghost-ghost_width-0.0001)*DH[0] || x>PI/4+(overghost-ghost_width-0.0001)*DH[0] ||
//	               y<-PI/4-(overghost-ghost_width-0.0001)*DH[1] || y>PI/4+(overghost-ghost_width-0.0001)*DH[1]   ) //0.0001 is for vertex center
	       if(x<-PI/4-(overghost-ghost_width-0.0001)*DH[0] || x>PI/4+(overghost-ghost_width-0.0001)*DH[0] ||
	          y<-PI/4-(overghost-ghost_width-0.0001)*DH[1] || y>PI/4+(overghost-ghost_width-0.0001)*DH[1]   )
	       {
		  double gx,gy,gz;
		  getglobalpox(gx,gy,gz,sPp->data->sst,x,y,z);
		  bool flag=false;
                  for(int i=0;i<cpusize;i++) 
		  {
		    flag=prolongpointstru(ss_src[i],true,sPp->data->sst,PatL,DH,CPatL,CDH,gx,gy,gz,Symmetry,i);
		    if(flag) break;
		  }
		  if(!flag)
                  {
	       	    if(myrank==0)
	       	    {
		        cout<<"ShellPatch::prolongpointstru fail to find shell source point for"<<endl;
		 	cout<<"sst = "<<sPp->data->sst<<" lx,ly,lz = "<<x<<","<<y<<","<<z<<endl;
                        MPI_Abort(MPI_COMM_WORLD,1);
	       	    }
	       	  }
	       }
	    }
	sPp = sPp->next;
     }
  if(myrank==0) cout<<"ShellPatch::setupintintstuff ss_src completes"<<endl;

     Pp=CPatL;
     while(Pp)
     {  
        double llb[dim],uub[dim];	     
        if(Symmetry>0)   llb[2]=Pp->data->bbox[2]-0.0001*CDH[2];
	else             llb[2]=Pp->data->bbox[2]+(CS_width+0.0001)*CDH[2];
        uub[2]=Pp->data->bbox[dim+2]-(CS_width+0.0001)*CDH[2];
        for(int j=0;j<2;j++)
	{
             if(Symmetry>1)   llb[j]=Pp->data->bbox[j]-0.0001*CDH[j];
	     else             llb[j]=Pp->data->bbox[j]+(CS_width+0.0001)*CDH[j];
	     uub[j]=Pp->data->bbox[dim+j]-(CS_width+0.0001)*CDH[j];
	}
	for(int iz=0;iz<Pp->data->shape[2];iz++)
	  for(int iy=0;iy<Pp->data->shape[1];iy++)
	    for(int ix=0;ix<Pp->data->shape[0];ix++)
	    {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
	       x=Pp->data->bbox[0]+ix*CDH[0];
	       y=Pp->data->bbox[1]+iy*CDH[1];
	       z=Pp->data->bbox[2]+iz*CDH[2];
#else
#ifdef Cell
	       x=Pp->data->bbox[0]+(ix+0.5)*CDH[0];
	       y=Pp->data->bbox[1]+(iy+0.5)*CDH[1];
	       z=Pp->data->bbox[2]+(iz+0.5)*CDH[2];
#else
#error Not define Vertex nor Cell
#endif  
#endif	       
	       if(x<llb[0] || x>uub[0] ||
                  y<llb[1] || y>uub[1] ||
                  z<llb[2] || z>uub[2])
	       {
		  int sst;
		  double lx,ly,lz;
		  bool flag=false;
                  getlocalpox(x,y,z,sst,lx,ly,lz);
                  for(int i=0;i<cpusize;i++) 
		  {
		    flag=prolongpointstru(cs_src[i],true,-1,PatL,DH,CPatL,CDH,x,y,z,Symmetry,i);
                    if(flag) break;
		  }
		  if(!flag)
                  {
	       	    if(myrank==0)
	       	    {
		        cout<<"ShellPatch::prolongpointstru fail to find shell source point for"<<endl;
		 	cout<<"sst = -1, x,y,z = "<<x<<","<<y<<","<<z<<endl;
                        MPI_Abort(MPI_COMM_WORLD,1);
	       	    }
	       	  }
	       }
	    }
        Pp=Pp->next;
     }
  if(myrank==0) cout<<"ShellPatch::setupintintstuff cs_src completes"<<endl;

  for(int i=0;i<cpusize;i++)
  {
     ps=ss_src[i];
     while(ps)
     {
	ts=ps->next;
        prolongpointstru(ss_dst[i],PatL,DH,CPatL,CDH,ps); //ps may be insterted more here
        ps=ts;
     }

     ps=cs_src[i];
     while(ps)
     {
	ts=ps->next;
        prolongpointstru(cs_dst[i],PatL,DH,CPatL,CDH,ps); //ps may be insterted more here
        ps=ts;
     }
  }
  if(myrank==0) cout<<"ShellPatch::setupintintstuff ss_dst and cs_dst complete"<<endl;

/*  
  for(int i=0;i<cpusize;i++)
  {
     ps=ss_src[i];
     ts=ss_dst[i];
     while(ps)
     {
      if(myrank==0) cout<<"src:"<<endl;
      check_pointstrul(ps,1);  
      if(myrank==0) cout<<"dst:"<<endl;
      check_pointstrul(ts,1);  
      ps=ps->next;
      ts=ts->next;
     }
  }
  exit(0);
*/  
}
#undef CS_width
void ShellPatch::setupcordtrans()
{
  MyList<ss_patch> *PP=PatL;
  while(PP)
  {
    PP->data->setupcordtrans();
    PP=PP->next;
  }
}
void ShellPatch::checkPatch()
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  if(myrank==0)
  {
     cout<<"belong to Shell Patchs "<<endl;
     MyList<ss_patch> *Pp=PatL;
     while(Pp)
     {
       cout<<" shape: [";
       for(int i=0;i<dim;i++) 
       {
         cout<<Pp->data->shape[i];
         if(i<dim-1) cout<<",";
         else        cout<<"]"<<endl;
       }
       cout<<"range:"<<"(";
       for(int i=0;i<dim;i++) 
       {
          cout<<Pp->data->bbox[i]<<":"<<Pp->data->bbox[dim+i];
          if(i<dim-1) cout<<",";
          else        cout<<")"<<endl;
       }
       Pp=Pp->next;
     }
  }
}
void ShellPatch::checkBlock(int sst)
{
  if(myrank==0)
  {
     cout<<"checking shell patch sst = "<<sst<<endl;
     MyList<ss_patch> *Pp=PatL;
     while(Pp)
     {
       if(Pp->data->sst == sst)
       {
        MyList<Block> *BP=Pp->data->blb;
        while(BP)
        {
    	  BP->data->checkBlock();
	  if(BP==Pp->data->ble) break;
	  BP=BP->next;
        }
       }
       Pp=Pp->next;
     }
  }
}
double ShellPatch::getdX(int dir)
{
    if(dir < 0 || dir >= dim)
    {
        cout<<"ShellPatch::getdX: error input dir = "<<dir<<", this Patch has direction (0,"<<dim-1<<")"<<endl;
        MPI_Abort(MPI_COMM_WORLD,1);
    }	    
    double h;
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif   
    if(shape[dir] == 1)
    {
        cout<<"ShellPatch::getdX: for direction "<<dir<<", this Patch has only one point. Can not determine dX for vertex center grid."<<endl;
        MPI_Abort(MPI_COMM_WORLD,1);
    } 
    if(dir<2) h = PI/2/(shape[dir]-1);
    else      h = (Rrange[1]-Rrange[0])/(shape[dir]-1);
#else
#ifdef Cell
    if(dir<2) h = PI/2/shape[dir];
    else      h = (Rrange[1]-Rrange[0])/shape[dir];
#else
#error Not define Vertex nor Cell
#endif  
#endif
    return h;
}
void ShellPatch::shellname(char* sn,int i)
{
   switch(i)
   {
      case 0:
        sprintf(sn,"zp");
	return;
      case 1:
        sprintf(sn,"zm");
	return;
      case 2:
        sprintf(sn,"xp");
	return;
      case 3:
        sprintf(sn,"xm");
	return;
      case 4:
        sprintf(sn,"yp");
	return;
      case 5:
        sprintf(sn,"ym");
	return;
   }
}
//Now we dump the data including overlap points
void ShellPatch::Dump_xyz(char *tag,double time,double dT)
{  
  MyList<ss_patch> *PP=PatL;
  while(PP)
  {
//   round at 4 and 5  
  int ncount=int(time/dT + 0.5);
 
   MPI_Status sta;
   int DIM=3;
   double llb[3],uub[3];
   double DX,DY,DZ;

   double *databuffer=0;
   if(myrank==0) 
   {
      databuffer=(double*)malloc(sizeof(double)*PP->data->shape[0]*PP->data->shape[1]*PP->data->shape[2]);
      if(!databuffer)
      {
       cout<<"ShellPatch::Dump_xyz: out of memory when dumping data."<<endl;
       MPI_Abort(MPI_COMM_WORLD,1);
      }
   }

   for(int DumpList=fngfs+gx;DumpList<=fngfs+gz;DumpList++)
   {     
     MyList<Block> *Bp=PP->data->blb;
     while(Bp)
     {
	Block *BP=Bp->data;
	if(BP->rank == 0 && myrank == 0)
	{
	 DX = BP->getdX(0);
	 DY = BP->getdX(1);
	 DZ = BP->getdX(2);
         llb[0] = (feq(BP->bbox[0],PP->data->bbox[0],DX/2)) ? BP->bbox[0] : BP->bbox[0]+ghost_width*DX;
         llb[1] = (feq(BP->bbox[1],PP->data->bbox[1],DY/2)) ? BP->bbox[1] : BP->bbox[1]+ghost_width*DY;
         llb[2] = (feq(BP->bbox[2],PP->data->bbox[2],DZ/2)) ? BP->bbox[2] : BP->bbox[2]+ghost_width*DZ;
         uub[0] = (feq(BP->bbox[3],PP->data->bbox[3],DX/2)) ? BP->bbox[3] : BP->bbox[3]-ghost_width*DX;
         uub[1] = (feq(BP->bbox[4],PP->data->bbox[4],DY/2)) ? BP->bbox[4] : BP->bbox[4]-ghost_width*DY;
         uub[2] = (feq(BP->bbox[5],PP->data->bbox[5],DZ/2)) ? BP->bbox[5] : BP->bbox[5]-ghost_width*DZ;
	 f_copy(DIM,PP->data->bbox,PP->data->bbox+DIM,PP->data->shape,databuffer,BP->bbox,BP->bbox+DIM,BP->shape,BP->fgfs[DumpList],llb,uub);
	}
	else
	{
	  int nnn=(BP->shape[0])*(BP->shape[1])*(BP->shape[2]);
	  if(myrank == 0)
	  {
	      double *bufferhere =(double*)malloc(sizeof(double)*nnn);
              if(!bufferhere)
              {
               cout<<"on node#"<<myrank<<", out of memory when dumping data."<<endl;
               MPI_Abort(MPI_COMM_WORLD,1);
              }
              MPI_Recv(bufferhere,nnn,MPI_DOUBLE,BP->rank,0,MPI_COMM_WORLD,&sta);
	      DX = BP->getdX(0);
     	      DY = BP->getdX(1);
     	      DZ = BP->getdX(2);
              llb[0] = (feq(BP->bbox[0],PP->data->bbox[0],DX/2)) ? BP->bbox[0] : BP->bbox[0]+ghost_width*DX;
              llb[1] = (feq(BP->bbox[1],PP->data->bbox[1],DY/2)) ? BP->bbox[1] : BP->bbox[1]+ghost_width*DY;
              llb[2] = (feq(BP->bbox[2],PP->data->bbox[2],DZ/2)) ? BP->bbox[2] : BP->bbox[2]+ghost_width*DZ;
              uub[0] = (feq(BP->bbox[3],PP->data->bbox[3],DX/2)) ? BP->bbox[3] : BP->bbox[3]-ghost_width*DX;
              uub[1] = (feq(BP->bbox[4],PP->data->bbox[4],DY/2)) ? BP->bbox[4] : BP->bbox[4]-ghost_width*DY;
              uub[2] = (feq(BP->bbox[5],PP->data->bbox[5],DZ/2)) ? BP->bbox[5] : BP->bbox[5]-ghost_width*DZ;
              f_copy(DIM,PP->data->bbox,PP->data->bbox+DIM,PP->data->shape,databuffer,BP->bbox,BP->bbox+DIM,BP->shape,bufferhere,llb,uub);
	      free(bufferhere);
	  }
	  else if(myrank == BP->rank)
	  {
            MPI_Send(BP->fgfs[DumpList],nnn,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
	  }
	}
	if(Bp == PP->data->ble) break;
	Bp=Bp->next;
     }
     if(myrank == 0)
     {
       char filename[100];
       char sn[3];
       shellname(sn,PP->data->sst);
       switch(DumpList-fngfs)
       {
	  case gx:
          if(tag)sprintf(filename,"%s_LevSH-%s_x_%05d.bin",tag,sn,ncount);
          else   sprintf(filename,   "LevSH-%s_x_%05d.bin"    ,sn,ncount);
	  break;
	  case gy:
          if(tag)sprintf(filename,"%s_LevSH-%s_y_%05d.bin",tag,sn,ncount);
          else   sprintf(filename,   "LevSH-%s_y_%05d.bin"    ,sn,ncount);
	  break;
	  case gz:
          if(tag)sprintf(filename,"%s_LevSH-%s_z_%05d.bin",tag,sn,ncount);
          else   sprintf(filename,   "LevSH-%s_z_%05d.bin"    ,sn,ncount);
	  break;
       }
   
       Parallel::writefile(time,PP->data->shape[0],PP->data->shape[1],PP->data->shape[2],
		      PP->data->bbox[0],PP->data->bbox[3],PP->data->bbox[1],PP->data->bbox[4],
	 	      PP->data->bbox[2],PP->data->bbox[5],filename,databuffer);
     }
  }

  if(myrank==0)  free(databuffer);

  PP=PP->next;
  }
}
void ShellPatch::Dump_Data(MyList<var> *DumpListi,char *tag,double time,double dT)
{
  MyList<ss_patch> *PP=PatL;
  while(PP)
  {
//   round at 4 and 5  
  int ncount=int(time/dT + 0.5);
 
   MPI_Status sta;
   int DIM=3;
   double llb[3],uub[3];
   double DX,DY,DZ;

   double *databuffer=0;
   if(myrank==0) 
   {
      databuffer=(double*)malloc(sizeof(double)*PP->data->shape[0]*PP->data->shape[1]*PP->data->shape[2]);
      if(!databuffer)
      {
       cout<<"ShellPatch::Dump_Data: out of memory when dumping data."<<endl;
       MPI_Abort(MPI_COMM_WORLD,1);
      }
   }
   
   MyList<var> *DumpList=DumpListi;
   while(DumpList)
   {
     var   *VP=DumpList->data;
     
     MyList<Block> *Bp=PP->data->blb;
     while(Bp)
     {
	Block *BP=Bp->data;
	if(BP->rank == 0 && myrank == 0)
	{
	 DX = BP->getdX(0);
	 DY = BP->getdX(1);
	 DZ = BP->getdX(2);
         llb[0] = (feq(BP->bbox[0],PP->data->bbox[0],DX/2)) ? BP->bbox[0] : BP->bbox[0]+ghost_width*DX;
         llb[1] = (feq(BP->bbox[1],PP->data->bbox[1],DY/2)) ? BP->bbox[1] : BP->bbox[1]+ghost_width*DY;
         llb[2] = (feq(BP->bbox[2],PP->data->bbox[2],DZ/2)) ? BP->bbox[2] : BP->bbox[2]+ghost_width*DZ;
         uub[0] = (feq(BP->bbox[3],PP->data->bbox[3],DX/2)) ? BP->bbox[3] : BP->bbox[3]-ghost_width*DX;
         uub[1] = (feq(BP->bbox[4],PP->data->bbox[4],DY/2)) ? BP->bbox[4] : BP->bbox[4]-ghost_width*DY;
         uub[2] = (feq(BP->bbox[5],PP->data->bbox[5],DZ/2)) ? BP->bbox[5] : BP->bbox[5]-ghost_width*DZ;
	 f_copy(DIM,PP->data->bbox,PP->data->bbox+DIM,PP->data->shape,databuffer,BP->bbox,BP->bbox+DIM,BP->shape,BP->fgfs[VP->sgfn],llb,uub);
	}
	else
	{
	  int nnn=(BP->shape[0])*(BP->shape[1])*(BP->shape[2]);
	  if(myrank == 0)
	  {
	      double *bufferhere =(double*)malloc(sizeof(double)*nnn);
              if(!bufferhere)
              {
               cout<<"on node#"<<myrank<<", out of memory when dumping data."<<endl;
               MPI_Abort(MPI_COMM_WORLD,1);
              }
              MPI_Recv(bufferhere,nnn,MPI_DOUBLE,BP->rank,0,MPI_COMM_WORLD,&sta);
	      DX = BP->getdX(0);
     	      DY = BP->getdX(1);
     	      DZ = BP->getdX(2);
              llb[0] = (feq(BP->bbox[0],PP->data->bbox[0],DX/2)) ? BP->bbox[0] : BP->bbox[0]+ghost_width*DX;
              llb[1] = (feq(BP->bbox[1],PP->data->bbox[1],DY/2)) ? BP->bbox[1] : BP->bbox[1]+ghost_width*DY;
              llb[2] = (feq(BP->bbox[2],PP->data->bbox[2],DZ/2)) ? BP->bbox[2] : BP->bbox[2]+ghost_width*DZ;
              uub[0] = (feq(BP->bbox[3],PP->data->bbox[3],DX/2)) ? BP->bbox[3] : BP->bbox[3]-ghost_width*DX;
              uub[1] = (feq(BP->bbox[4],PP->data->bbox[4],DY/2)) ? BP->bbox[4] : BP->bbox[4]-ghost_width*DY;
              uub[2] = (feq(BP->bbox[5],PP->data->bbox[5],DZ/2)) ? BP->bbox[5] : BP->bbox[5]-ghost_width*DZ;
              f_copy(DIM,PP->data->bbox,PP->data->bbox+DIM,PP->data->shape,databuffer,BP->bbox,BP->bbox+DIM,BP->shape,bufferhere,llb,uub);
	      free(bufferhere);
	  }
	  else if(myrank == BP->rank)
	  {
            MPI_Send(BP->fgfs[VP->sgfn],nnn,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
	  }
	}
	if(Bp == PP->data->ble) break;
	Bp=Bp->next;
     }
     if(myrank == 0)
     {
       char filename[100];
       char sn[3];
       shellname(sn,PP->data->sst);
       if(tag)sprintf(filename,"%s_LevSH-%s_%s_%05d.bin",tag,sn,VP->name,ncount);
       else   sprintf(filename,   "LevSH-%s_%s_%05d.bin"    ,sn,VP->name,ncount);
   
       Parallel::writefile(time,PP->data->shape[0],PP->data->shape[1],PP->data->shape[2],
		      PP->data->bbox[0],PP->data->bbox[3],PP->data->bbox[1],PP->data->bbox[4],
	 	      PP->data->bbox[2],PP->data->bbox[5],filename,databuffer);
     }
     DumpList=DumpList->next;
  }

  if(myrank==0)  free(databuffer);

  PP=PP->next;
  }
}
void ShellPatch::intertransfer(MyList<pointstru> **src,MyList<pointstru> **dst,
		               MyList<var> *VarList1 /* source */,MyList<var> *VarList2 /*target */,
			       int Symmetry)
{
    int myrank,cpusize;
    MPI_Comm_size(MPI_COMM_WORLD,&cpusize);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

    int node;
    
    MPI_Request *reqs;
    MPI_Status  *stats; 
    reqs=new MPI_Request[2*cpusize];
    stats=new MPI_Status[2*cpusize];  
    int req_no=0;

    double **send_data,**rec_data;
    send_data=new double*[cpusize];
     rec_data=new double*[cpusize];
    int length;

    for(node=0; node<cpusize; node++)
    {
      send_data[node]=rec_data[node]=0;
      if(node == myrank)
      {       
	if(length=interdata_packer(0,src[myrank],dst[myrank],node,PACK,VarList1,VarList2,Symmetry))
	{
         rec_data[node]=new double[length];
	 if(!rec_data[node])
	 {
          cout<<"out of memory when new in short transfer, place 1"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	 }
         interdata_packer(rec_data[node],src[myrank],dst[myrank],node,PACK,VarList1,VarList2,Symmetry);
       }
      }
      else
      {
//send from this cpu to cpu#node    
       if(length=interdata_packer(0,src[myrank],dst[myrank],node,PACK,VarList1,VarList2,Symmetry))
        {
         send_data[node]=new double[length];
	 if(!send_data[node])
	 {
          cout<<"out of memory when new in short transfer, place 2"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	 }
         interdata_packer(send_data[node],src[myrank],dst[myrank],node,PACK,VarList1,VarList2,Symmetry);
       	 MPI_Isend((void*)send_data[node],length,MPI_DOUBLE,node,1,MPI_COMM_WORLD,reqs+req_no++);
        }
//receive from cpu#node to this cpu
       if(length=interdata_packer(0,src[node],dst[node],node,UNPACK,VarList1,VarList2,Symmetry))
       {
	 rec_data[node]=new double[length];
	 if(!rec_data[node])
	 {
          cout<<"out of memory when new in short transfer, place 3"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	 }
	 MPI_Irecv((void*)rec_data[node],length,MPI_DOUBLE,node,1,MPI_COMM_WORLD,reqs+req_no++);
       }
      }
    }
//wait for all requests to complete
   MPI_Waitall(req_no,reqs,stats);

   for(node=0; node<cpusize; node++)
       if(rec_data[node])
         interdata_packer(rec_data[node],src[node],dst[node],node,UNPACK,VarList1,VarList2,Symmetry);

   for(node=0; node<cpusize; node++)
   {
     if(send_data[node]) delete[] send_data[node];
     if(rec_data[node])  delete[] rec_data[node];
   }

   delete[] reqs; delete[] stats; 
   delete[] send_data; delete[] rec_data;
}
//   PACK: prepare target data in 'data'
// UNPACK: copy target data from 'data' to corresponding numerical grids
int ShellPatch::interdata_packer(double *data,MyList<pointstru> *src,MyList<pointstru> *dst, int rank_in, int dir,
                                 MyList<var> *VarLists /* source */,MyList<var> *VarListd /* target */,int Symmetry)
{
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

   int DIM = dim;
   int ordn=2*ghost_width;
  
   if(dir!=PACK && dir!=UNPACK)
   { 
    cout << "error dir " <<dir<<" for data_packer "<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
   }

   int size_out=0;

   if(!src || !dst) return size_out;

   MyList<var> *varls,*varld;

   varls=VarLists;
   varld=VarListd;
   while(varls && varld)
   {
     varls=varls->next;
     varld=varld->next;
   }

   if(varls || varld)
   {
     cout<<"error in short data packer, var lists does not match."<<endl;
     MPI_Abort(MPI_COMM_WORLD,1);
   }
   
   while(src && dst)
   {
    if((dir==PACK && dst->data->Bg->rank==rank_in && src->data->Bg->rank==myrank) ||
       (dir==UNPACK && src->data->Bg->rank==rank_in && dst->data->Bg->rank==myrank))
    {
    varls=VarLists;
    varld=VarListd;
    while(varls && varld)
    {
     if(data)
     {
     if(dir==PACK)
     {
/*	     
         f_global_interp(src->data->Bg->shape,src->data->Bg->X[0],src->data->Bg->X[1],src->data->Bg->X[2],
			 src->data->Bg->fgfs[varls->data->sgfn],data[size_out],
			 src->data->lpox[0],src->data->lpox[1],src->data->lpox[2],ordn,varls->data->SoA,Symmetry);
*/
	 int DIMh=(src->data->dumyd==-1)?dim:1;
         if(src->data->coef==0)
	 {
	    src->data->coef=new double[ordn*DIMh];
	    src->data->sind=new int[dim];
	    if(DIMh==3)
	    {
	      for(int i=0;i<DIMh;i++)
	      {
	       double dd=src->data->Bg->getdX(i);
	       //0.001 instead of 0.4 makes the point locate more center
	       src->data->sind[i] = int((src->data->lpox[i]-src->data->Bg->X[i][0])/dd)-ordn/2+1;  
	       double h1,h2;
	       for(int j=0;j<ordn;j++)
	       {
		  h1=src->data->Bg->X[i][0]+(src->data->sind[i]+j)*dd;
 	          src->data->coef[i*ordn+j]=1;
		  for(int k=0;k<j;k++)
		  {
		     h2=src->data->Bg->X[i][0]+(src->data->sind[i]+k)*dd;
		     src->data->coef[i*ordn+j] *= (src->data->lpox[i]-h2)/(h1-h2);
		  }
		  for(int k=j+1;k<ordn;k++) 
		  {
		     h2=src->data->Bg->X[i][0]+(src->data->sind[i]+k)*dd;
		     src->data->coef[i*ordn+j] *= (src->data->lpox[i]-h2)/(h1-h2);
		  }
	       }
	      }
	    }
	    else
	    {
	       int actd=1-src->data->dumyd;
	       double dd=src->data->Bg->getdX(actd);
	       src->data->sind[0] = int((src->data->lpox[actd]-src->data->Bg->X[actd][0])/dd)-ordn/2+1;
	       double h1,h2;
	       for(int j=0;j<ordn;j++)
	       {
		  h1=src->data->Bg->X[actd][0]+(src->data->sind[0]+j)*dd;
 	          src->data->coef[j]=1;
		  for(int k=0;k<j;k++)
		  {
		     h2=src->data->Bg->X[actd][0]+(src->data->sind[0]+k)*dd;
		     src->data->coef[j] *= (src->data->lpox[actd]-h2)/(h1-h2);
		  }
		  for(int k=j+1;k<ordn;k++) 
		  {
		     h2=src->data->Bg->X[actd][0]+(src->data->sind[0]+k)*dd;
		     src->data->coef[j] *= (src->data->lpox[actd]-h2)/(h1-h2);
		  }
	       }
	      src->data->sind[2] = int((src->data->lpox[2]-src->data->Bg->X[2][0])/src->data->Bg->getdX(2)+0.001);
	      if(!feq(src->data->Bg->X[2][src->data->sind[2]],src->data->lpox[2],src->data->Bg->getdX(2)/2000))
		 cout<<"error in ShellPatch::interdata_packer point = "<<src->data->lpox[2]<<" != grid "<<src->data->Bg->X[2][src->data->sind[2]]<<endl;
	      src->data->sind[1] = int((src->data->lpox[src->data->dumyd]-src->data->Bg->X[src->data->dumyd][0])/
			                src->data->Bg->getdX(src->data->dumyd)+0.001);
	      if(!feq(src->data->Bg->X[src->data->dumyd][src->data->sind[1]],src->data->lpox[src->data->dumyd],src->data->Bg->getdX(src->data->dumyd)/2000))
		 cout<<"error in ShellPatch::interdata_packer for dumy dimension point = "
		     <<src->data->lpox[src->data->dumyd]<<" != grid "<<src->data->Bg->X[src->data->dumyd][src->data->sind[1]]<<endl;
	    }
	 }
//interpolate
         switch(DIMh)
	 {
	  case 3:
            f_global_interpind(src->data->Bg->shape,src->data->Bg->X[0],src->data->Bg->X[1],src->data->Bg->X[2],
			 src->data->Bg->fgfs[varls->data->sgfn],data[size_out],
			 src->data->lpox[0],src->data->lpox[1],src->data->lpox[2],ordn,varls->data->SoA,Symmetry,
			 src->data->sind,src->data->coef,src->data->ssst);
	    break;
	  case 2:
            f_global_interpind2d(src->data->Bg->shape,src->data->Bg->X[0],src->data->Bg->X[1],src->data->Bg->X[2],
			 src->data->Bg->fgfs[varls->data->sgfn],data[size_out],
			 src->data->lpox[0],src->data->lpox[1],src->data->lpox[2],ordn,varls->data->SoA,Symmetry,
			 src->data->sind,src->data->coef,src->data->ssst);
	    break;
	  case 1:
            f_global_interpind1d(src->data->Bg->shape,src->data->Bg->X[0],src->data->Bg->X[1],src->data->Bg->X[2],
			 src->data->Bg->fgfs[varls->data->sgfn],data[size_out],
			 src->data->lpox[0],src->data->lpox[1],src->data->lpox[2],ordn,varls->data->SoA,Symmetry,
			 src->data->sind,src->data->coef,src->data->ssst,src->data->dumyd);
	    break;
	  default:
	    cout<<"ShellPatch::interdata_packer: not recognized DIM = "<<DIMh<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	 }
     }
     if(dir==UNPACK) //from target data to corresponding grid
         f_pointcopy(DIM,dst->data->Bg->bbox,dst->data->Bg->bbox+dim,dst->data->Bg->shape,dst->data->Bg->fgfs[varld->data->sgfn],
		     dst->data->lpox[0],dst->data->lpox[1],dst->data->lpox[2],data[size_out]);
     }
     size_out+=1;
     varls=varls->next;
     varld=varld->next;
    }
    }
    dst=dst->next;
    src=src->next;
   }

   return size_out;
}
void ShellPatch::Synch(MyList<var> *VarList,int Symmetry)
{
  MyList<ss_patch> *Pp = PatL;
  while(Pp)
  {
     Pp->data->Sync(VarList,Symmetry);
     Pp=Pp->next;
  }
	
  intertransfer(ss_src,ss_dst,VarList,VarList,Symmetry);	
}
void ShellPatch::CS_Inter(MyList<var> *VarList,int Symmetry)
{
  intertransfer(cs_src,cs_dst,VarList,VarList,Symmetry);	
}
void ShellPatch::check_pointstrul(MyList<pointstru>* pp,bool first_only)
{
  int myrank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  if(myrank==0)
  {
   if(!pp) cout<<"ShellPatch::check_pointstrul meets empty pointstru"<<endl;
   else    cout<<"checking check_pointstrul..."<<endl;
   while(pp)
   {
    if(pp->data->Bg) cout<<"on node#"<<pp->data->Bg->rank<<endl;
    else             cout<<"virtual pointstru"<<endl;
    cout<<"source sst = "<<pp->data->ssst<<endl;
    cout<<"target sst = "<<pp->data->tsst<<endl;
    cout<<"dumy dimension = "<<pp->data->dumyd<<endl;
    cout<<"global coordinates: (";
    for(int i=0;i<dim;i++) 
    {
       if(i<dim-1) cout<<pp->data->gpox[i]<<",";
       else        cout<<pp->data->gpox[i]<<")"<<endl;
    }
    cout<<"local coordinates: (";
    for(int i=0;i<dim;i++) 
    {
       if(i<dim-1) cout<<pp->data->lpox[i]<<",";
       else        cout<<pp->data->lpox[i]<<")"<<endl;
    }
    if(first_only) return;
    pp=pp->next;
   } 
  }
}
void ShellPatch::check_pointstrul2(MyList<pointstru>* pp,int first_last_only)
{
  int myrank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  if(myrank==0)
  {
   if(!pp) cout<<"ShellPatch::check_pointstrul meets empty pointstru"<<endl;
   else    cout<<"checking check_pointstrul..."<<endl;
   while(pp)
   {
   if(first_last_only==2)
   {
     if(pp->next==0)
     {
    if(pp->data->Bg) cout<<"on node#"<<pp->data->Bg->rank<<endl;
    else             cout<<"virtual pointstru"<<endl;
    cout<<"source sst = "<<pp->data->ssst<<endl;
    cout<<"target sst = "<<pp->data->tsst<<endl;
    cout<<"dumy dimension = "<<pp->data->dumyd<<endl;
    cout<<"global coordinates: (";
    for(int i=0;i<dim;i++) 
    {
       if(i<dim-1) cout<<pp->data->gpox[i]<<",";
       else        cout<<pp->data->gpox[i]<<")"<<endl;
    }
    cout<<"local coordinates: (";
    for(int i=0;i<dim;i++) 
    {
       if(i<dim-1) cout<<pp->data->lpox[i]<<",";
       else        cout<<pp->data->lpox[i]<<")"<<endl;
    }
     }
   }
   else
   {
    if(pp->data->Bg) cout<<"on node#"<<pp->data->Bg->rank<<endl;
    else             cout<<"virtual pointstru"<<endl;
    cout<<"source sst = "<<pp->data->ssst<<endl;
    cout<<"target sst = "<<pp->data->tsst<<endl;
    cout<<"dumy dimension = "<<pp->data->dumyd<<endl;
    cout<<"global coordinates: (";
    for(int i=0;i<dim;i++) 
    {
       if(i<dim-1) cout<<pp->data->gpox[i]<<",";
       else        cout<<pp->data->gpox[i]<<")"<<endl;
    }
    cout<<"local coordinates: (";
    for(int i=0;i<dim;i++) 
    {
       if(i<dim-1) cout<<pp->data->lpox[i]<<",";
       else        cout<<pp->data->lpox[i]<<")"<<endl;
    }
    if(first_last_only==1) return;
   }
    pp=pp->next;
   } 
  }
}
void ShellPatch::matchcheck(MyList<Patch> *CPatL)
{
  double cbd=CPatL->data->bbox[dim];
  for(int i=1;i<dim;i++) cbd = min(cbd,CPatL->data->bbox[dim+i]);
  cbd = cbd - getsr(Rrange[0]);
  double dr,dc;
  dc = CPatL->data->getdX(0);
  dr = getdX(2);
  for(int i=1;i<dim;i++)
  {
    dc = max(dc,CPatL->data->getdX(i));
//    dr = max(dr,getdX(i));
  }

  int ir,ic;
  ir=int(cbd/dr);
  ic=int(cbd/dc);
  if(min(ir,ic)<3*ghost_width)
  {
    int myrank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    if(myrank==0)
    {
       cout<<"Shell Patches insert too shallow:"<<endl;
       cout<<"distantance between these two boundaries is "<<cbd<<", spatial step is "<<max(dc,dr)<<endl;
       MPI_Abort(MPI_COMM_WORLD,1);
    }
  }
}
void ShellPatch::Interp_Points(MyList<var> *VarList,
		int NN,double **XX,  /*input global Cartesian coordinate*/
		double *Shellf,int Symmetry)
{
// NOTE: we do not Synchnize variables here, make sure of that before calling this routine
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

  int ordn=2*ghost_width;
  MyList<var> *varl;
  int num_var=0;
  varl=VarList;
  while(varl)
  {
    num_var++;
    varl=varl->next;
  }
  
  double *shellf;
  shellf=new double[NN*num_var];
  memset(shellf,0,sizeof(double)*NN*num_var);

// we use weight to monitor code, later some day we can move it for optimization  
  int *weight;
  weight=new int[NN];
  memset(weight,0,sizeof(int)*NN);

  double *DH,*llb,*uub;
  DH = new double[dim];
	
  for(int i=0;i<dim;i++)
  {
    DH[i] = getdX(i);
  }
  llb = new double[dim];
  uub = new double[dim];

  for(int j=0; j<NN; j++) // run along points
  {
    double pox[dim];
    int sst;
    getlocalpox(XX[0][j],XX[1][j],XX[2][j],sst,pox[0],pox[1],pox[2]);

    MyList<ss_patch> *sPp=PatL;
    while(sPp->data->sst != sst) sPp = sPp->next;

    if(myrank == 0 && ((!sPp) || pox[2] < Rrange[0] || pox[2] > Rrange[1]))
    {
      cout<<"ShellPatch::Interp_Points: point (";
      for(int k=0;k<dim;k++)
      {
        cout<<XX[k][j];
        if(k < dim-1) cout<<",";
        else          cout<<") is out of the ShellPatch."<<endl;
      }
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    if(!sPp) return;

     MyList<Block> *Bp=sPp->data->blb;
     bool notfind=true;
     while(notfind && Bp) // run along Blocks
     {
	Block *BP=Bp->data;
	  
	bool flag=true;
	for(int i=0;i<dim;i++)
	{
// NOTE: our dividing structure is (exclude ghost)
// -1 0
//       1  2
// so (0,1) does not belong to any part for vertex structure
// here we put (0,0.5) to left part and (0.5,1) to right part
// BUT for cell structure the bbox is (-1.5,0.5) and (0.5,2.5), there is no missing region at all
//
// because of getlocalpox, pox will not goes into overghost region of ss_patch
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
          llb[i] = (feq(BP->bbox[i]    ,sPp->data->bbox[i]    ,DH[i]/2)) ? BP->bbox[i]     : BP->bbox[i]    +(ghost_width-0.5)*DH[i];
          uub[i] = (feq(BP->bbox[dim+i],sPp->data->bbox[dim+i],DH[i]/2)) ? BP->bbox[dim+i] : BP->bbox[dim+i]-(ghost_width-0.5)*DH[i];
#else
#ifdef Cell
          llb[i] = (feq(BP->bbox[i]    ,sPp->data->bbox[i]    ,DH[i]/2)) ? BP->bbox[i]     : BP->bbox[i]    +ghost_width*DH[i];
          uub[i] = (feq(BP->bbox[dim+i],sPp->data->bbox[dim+i],DH[i]/2)) ? BP->bbox[dim+i] : BP->bbox[dim+i]-ghost_width*DH[i];
#else
#error Not define Vertex nor Cell
#endif  
#endif 
	 if(pox[i]-llb[i]<-DH[i]/2 || pox[i]-uub[i]>DH[i]/2) {flag = false; break;}
	}

	if(flag)
	{
         notfind = false;
         if(myrank == BP->rank)
         {
//---> interpolation
	   varl=VarList;
	   int k=0;
	   while(varl) // run along variables
	   {
              f_global_interp_ss(BP->shape,BP->X[0],BP->X[1],BP->X[2],BP->fgfs[varl->data->sgfn],shellf[j*num_var+k],
			         pox[0],pox[1],pox[2],ordn,varl->data->SoA,Symmetry,sst);
	      varl=varl->next;
	      k++;
	   }
	   weight[j] = 1;
	 }
	}
	if(Bp == sPp->data->ble) break;
	Bp=Bp->next;
     }
  }
  
  MPI_Allreduce(shellf,Shellf,NN*num_var,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  int *Weight;
  Weight = new int[NN];
  MPI_Allreduce(weight,Weight,NN,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

  for(int i=0;i<NN;i++)
  {
     if(Weight[i] > 1)
     {
       if(myrank == 0) cout<<"WARNING: ShellPatch::Interp_Points meets multiple weight"<<endl;
       for(int j=0;j<num_var;j++)
	 Shellf[j+i*num_var] = Shellf[j+i*num_var]/Weight[i];
     }
     else if(Weight[i] == 0 && myrank == 0)
     {
       cout<<"ERROR: ShellPatch::Interp_Points fails to find point (";
       for(int j=0;j<dim;j++)
       {
	  cout<<XX[j][i];
	  if(j<dim-1) cout<<",";
	  else        cout<<")";
       }
       cout<<" on ShellPatch ("<<Rrange[0]<<","<<Rrange[1]<<endl;

  cout<<"splited domains:"<<endl;
  MyList<ss_patch> *sPp=PatL;
  while(sPp)
  {
     char sn[3];
     shellname(sn,sPp->data->sst);
     cout<<"ss_patch "<<sn<<":"<<endl;
     MyList<Block> *Bp=sPp->data->blb;
     while(Bp)
     {
	Block *BP=Bp->data;

	for(int i=0;i<dim;i++)
	{
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
          llb[i] = (feq(BP->bbox[i]    ,sPp->data->bbox[i]    ,DH[i]/2)) ? BP->bbox[i]     : BP->bbox[i]    +(ghost_width-0.5)*DH[i];
          uub[i] = (feq(BP->bbox[dim+i],sPp->data->bbox[dim+i],DH[i]/2)) ? BP->bbox[dim+i] : BP->bbox[dim+i]-(ghost_width-0.5)*DH[i];
#else
#ifdef Cell
          llb[i] = (feq(BP->bbox[i]    ,sPp->data->bbox[i]    ,DH[i]/2)) ? BP->bbox[i]     : BP->bbox[i]    +ghost_width*DH[i];
          uub[i] = (feq(BP->bbox[dim+i],sPp->data->bbox[dim+i],DH[i]/2)) ? BP->bbox[dim+i] : BP->bbox[dim+i]-ghost_width*DH[i];
#else
#error Not define Vertex nor Cell
#endif  
#endif 
	}       
       cout<<"(";
       for(int j=0;j<dim;j++)
       {
	  cout<<llb[j]<<":"<<uub[j];
	  if(j<dim-1) cout<<",";
	  else        cout<<")"<<endl;
       }
	if(Bp == sPp->data->ble) break;
	Bp=Bp->next;
     }
     sPp=sPp->next;
  }    
       MPI_Abort(MPI_COMM_WORLD,1);
     }
  }

  delete[] shellf;
  delete[] weight;
  delete[] Weight;
  delete[] DH;
  delete[] llb;
  delete[] uub;
}
