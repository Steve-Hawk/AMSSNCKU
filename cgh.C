//$Id: cgh.C,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
#ifdef newc
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
using namespace std;
#else
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#endif

#include <mpi.h>

#include "microdef.h"
#include "misc.h"
#include "cgh.h"
#include "Parallel.h"

cgh::cgh(int ingfsi,int fngfsi,int Symmetry,char *filename,int checkrun,
         monitor *ErrorMonitor):ingfs(ingfsi),fngfs(fngfsi)
{
   read_bbox(Symmetry,filename);
   if(!checkrun) sethandle(ErrorMonitor);
   for(int lev=0;lev<levels;lev++)  PatL[lev] = construct_patchlist(lev,Symmetry);
}
  cgh::~cgh()
{
  for(int lev=0;lev<levels;lev++)
  {
    for(int grd=0;grd<grids[lev];grd++)
    {
       delete[] bbox[lev][grd];
       delete[] shape[lev][grd];
       delete[] handle[lev][grd];
    }
    delete[] bbox[lev];
    delete[] shape[lev];
    delete[] handle[lev];
    Parallel::KillBlocks(PatL[lev]);
    PatL[lev]->destroyList();
  }
  delete[] grids;
  delete[] bbox;
  delete[] shape;
  delete[] handle;
  delete[] PatL;
}
void cgh::compose_cgh(int nprocs)
{
  for(int lev=0;lev<levels;lev++)
  {
    checkPatchList(PatL[lev],false);
    Parallel::distribute(PatL[lev],nprocs,ingfs,fngfs,false);
  }
}
void cgh::sethandle(monitor *ErrorMonitor)
{
    int BH_num;
    double **Porg;
    char filename[100];
    strcpy(filename,"input.par");
// read parameter from file
    {
      const int LEN = 256;   
      char pline[LEN];
      string str, sgrp, skey, sval;
      int sind;
      ifstream inf(filename,  ifstream::in);
      if ( !inf.good() && ErrorMonitor->outfile) 
      {
        ErrorMonitor->outfile<<"Can not open parameter file "<<filename<<" for inputing information of black holes"<<endl;
       	MPI_Abort(MPI_COMM_WORLD,1);
      }

      for( int i=1; inf.good(); i++)  
      {
        inf.getline(pline, LEN);
        str = pline;
     
        int status = misc::parse_parts(str, sgrp, skey, sval, sind);
        if ( status == -1 ) 
	{
	  if(ErrorMonitor->outfile) ErrorMonitor->outfile<<"error reading parameter file "<<filename<<" in line "<<i<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	}
        else if( status == 0 ) continue;
   
	if(sgrp == "BSSN" && skey == "BH_num") BH_num = atoi(sval.c_str());
	else if(sgrp == "cgh" && skey == "moving levels start from") 
	{
	   movls = atoi(sval.c_str()); 
	   movls = min(movls,levels);
	   movls = max(0,movls);
	}
      }
      inf.close();
    }
    Porg = new double*[BH_num];
    for(int i=0;i<BH_num;i++) Porg[i]=new double[dim];
// read parameter from file
    {
      const int LEN = 256;   
      char pline[LEN];
      string str, sgrp, skey, sval;
      int sind;
      ifstream inf(filename,  ifstream::in);
      if ( !inf.good() && ErrorMonitor->outfile) 
      {
        ErrorMonitor->outfile<<"Can not open parameter file "<<filename
		                      <<" for inputing information of black holes"<<endl;
       	MPI_Abort(MPI_COMM_WORLD,1);
      }

      for( int i=1; inf.good(); i++)  
      {
        inf.getline(pline, LEN);
        str = pline;
     
        int status = misc::parse_parts(str, sgrp, skey, sval, sind);
        if ( status == -1 ) 
	{
	  if(ErrorMonitor->outfile) ErrorMonitor->outfile<<"error reading parameter file "<<filename<<" in line "<<i<<endl; 
	  MPI_Abort(MPI_COMM_WORLD,1);
	}
        else if( status == 0 ) continue;
   
	if(sgrp == "BSSN" && sind < BH_num)
	{
	 if ( skey == "Porgx")  Porg[sind][0] = atof(sval.c_str());
	 else if ( skey == "Porgy")  Porg[sind][1] = atof(sval.c_str());
	 else if ( skey == "Porgz")  Porg[sind][2]= atof(sval.c_str());
	}
      }
      inf.close();
    }

   for(int lev=0;lev<movls;lev++)
     for(int grd=0;grd<grids[lev];grd++)
       for(int i=0;i<dim;i++) handle[lev][grd][i]=0;

  if(movls<levels)
  {
   if(ErrorMonitor->outfile) cout<<"moving levels are lev#"<<movls<<"--"<<levels-1<<endl;
   for(int lev=movls;lev<levels;lev++)
     for(int grd=0;grd<grids[lev];grd++)
     {
	 int bht=0;
	 for(int bhi=0;bhi<BH_num;bhi++)
	 {
	    bool flag=false;
	    for(int i=0;i<dim;i++)
               if(Porg[bhi][i] < bbox[lev][grd][i] || Porg[bhi][i] > bbox[lev][grd][i+dim]) {flag=true; break;}
	    if(flag) continue;
	    bht++;
	    if(bht==1)  for(int i=0;i<dim;i++) handle[lev][grd][i]=Porg[bhi][i];
	    else if(ErrorMonitor->outfile) 
	    {
               ErrorMonitor->outfile<<"cgh::sethandle: lev#"<<lev<<" grd#"<<grd<<" has too many black holes"<<endl;
               MPI_Abort(MPI_COMM_WORLD,1);
	    }
	 }
     }     
  }
  else if(ErrorMonitor->outfile) 
  {
   if(levels>1) cout<<"fixed mesh refinement!"<<endl;
   else         cout<<"unigrid simulation!"<<endl;
  }
    
  for(int i=0;i<BH_num;i++) delete[] Porg[i];
  delete[] Porg;
}
void cgh::checkPatchList(MyList<Patch> *PatL,bool buflog)
{
   while(PatL)
   {
     PatL->data->checkPatch(buflog);
     PatL = PatL->next;
   }
}
void cgh::Regrid(int Symmetry,int BH_num,double **Porgbr,double**Porg0,
		 MyList<var> *OldList, MyList<var> *StateList, 
		 MyList<var> *FutureList,MyList<var> *tmList,bool BB,
		 monitor *ErrorMonitor)
{
// for moving part   
  if(movls<levels)
  {
   bool tot_flag=false;
   bool *lev_flag;
   double **tmpPorg;
   tmpPorg=new double*[BH_num];
   for(int bhi=0;bhi<BH_num;bhi++)
   {
      tmpPorg[bhi]=new double[dim];
      for(int i=0;i<dim;i++) tmpPorg[bhi][i]=Porgbr[bhi][i];
   }
   lev_flag = new bool[levels-movls];
   for(int lev=movls;lev<levels;lev++)
   {
   lev_flag[lev-movls]=false;
   for(int grd=0;grd<grids[lev];grd++)
   {
    int flag;
    int do_every=2;
    double dX = PatL[lev]->data->blb->data->getdX(0);
    double dY = PatL[lev]->data->blb->data->getdX(1);
    double dZ = PatL[lev]->data->blb->data->getdX(2);
    double rr;
// make sure that the grid corresponds to the black hole
    int bhi=0;
    for(bhi=0;bhi<BH_num;bhi++)
    {
// because finner level may also change Porgbr, so we need factor 2	    
      if(feq(Porgbr[bhi][0],handle[lev][grd][0],2*do_every*dX) && 
	 feq(Porgbr[bhi][1],handle[lev][grd][1],2*do_every*dY) &&
	 feq(Porgbr[bhi][2],handle[lev][grd][2],2*do_every*dZ)) break;
    }
    if(bhi == BH_num)
    {
      if(ErrorMonitor->outfile)
      {
	 ErrorMonitor->outfile<<"cgh::Regrid: no black hole matches with grid lev#"<<lev<<" grd#"<<grd
	     <<" with handle ("<<handle[lev][grd][0]<<","<<handle[lev][grd][1]<<","<<handle[lev][grd][2]<<")"<<endl;
	 ErrorMonitor->outfile<<"black holes' old positions:"<<endl;
	 for(bhi=0;bhi<BH_num;bhi++) ErrorMonitor->outfile<<"#"<<bhi<<": ("<<Porgbr[bhi][0]<<","<<Porgbr[bhi][1]<<","<<Porgbr[bhi][2]<<")"<<endl;
         MPI_Abort(MPI_COMM_WORLD,1);
      }

      delete[] lev_flag;
      for(bhi=0;bhi<BH_num;bhi++) delete[] tmpPorg[bhi];
      delete[] tmpPorg;
      return;
    }
// x direction
    rr = (Porg0[bhi][0]-handle[lev][grd][0])/dX;
    if(rr>0)flag = int(rr+0.5)/do_every;
    else    flag = int(rr-0.5)/do_every;
    flag = flag*do_every;
    rr = bbox[lev][grd][0]+flag*dX;
// pay attention to the symmetric case
    if(Symmetry == 2 && rr < 0) rr=-bbox[lev][grd][0];
    else                        rr=flag*dX;

    if(fabs(rr)>dX/2)
    {   
     lev_flag[lev-movls]=tot_flag=true;
     bbox[lev][grd][0]=bbox[lev][grd][0]+rr;
     bbox[lev][grd][3]=bbox[lev][grd][3]+rr;
     handle[lev][grd][0]+= rr;
     tmpPorg[bhi][0] = Porg0[bhi][0];
    }
    
// y direction       
    rr = (Porg0[bhi][1]-handle[lev][grd][1])/dY;
    if(rr>0)flag = int(rr+0.5)/do_every;
    else    flag = int(rr-0.5)/do_every;
    flag = flag*do_every;
    rr = bbox[lev][grd][1]+flag*dY;
// pay attention to the symmetric case
    if(Symmetry == 2 && rr < 0) rr=-bbox[lev][grd][1];
    else                        rr=flag*dY;

    if(fabs(rr)>dY/2)
    {
     lev_flag[lev-movls]=tot_flag=true;
     bbox[lev][grd][1]=bbox[lev][grd][1]+rr;
     bbox[lev][grd][4]=bbox[lev][grd][4]+rr;
     handle[lev][grd][1]+= rr;
     tmpPorg[bhi][1] = Porg0[bhi][1];
    }

// z direction    
    rr = (Porg0[bhi][2]-handle[lev][grd][2])/dZ;
    if(rr>0)flag = int(rr+0.5)/do_every;
    else    flag = int(rr-0.5)/do_every;
    flag = flag*do_every;
    rr = bbox[lev][grd][2]+flag*dZ;
// pay attention to the symmetric case
    if(Symmetry > 0 && rr < 0) rr=-bbox[lev][grd][1];
    else                       rr=flag*dZ;

    if(fabs(rr)>dZ/2)
    {
     lev_flag[lev-movls]=tot_flag=true;
     bbox[lev][grd][2]=bbox[lev][grd][2]+rr;
     bbox[lev][grd][5]=bbox[lev][grd][5]+rr;
     handle[lev][grd][2]+= rr;
     tmpPorg[bhi][2] = Porg0[bhi][2];
    }
   }
//   if(ErrorMonitor->outfile && lev_flag[lev-movls]) cout<<"lev#"<<lev<<"'s boxes moved"<<endl;
   }

   if(tot_flag) 
   {
     int nprocs; 
     MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
     recompose_cgh(nprocs,lev_flag,OldList,StateList,FutureList,tmList,Symmetry,BB);
     for(int bhi=0;bhi<BH_num;bhi++)
     {
        for(int i=0;i<dim;i++) Porgbr[bhi][i]=tmpPorg[bhi][i];
     }
   }

   delete[] lev_flag;
   for(int bhi=0;bhi<BH_num;bhi++) delete[] tmpPorg[bhi];
   delete[] tmpPorg;
  }
}
void cgh::recompose_cgh(int nprocs,bool *lev_flag,
		        MyList<var> *OldList, MyList<var> *StateList, 
			MyList<var> *FutureList,MyList<var> *tmList,
			int Symmetry, bool BB)
{
  for(int lev=movls;lev<levels;lev++)
    if(lev_flag[lev-movls])
    {
     MyList<Patch> *tmPat=0;
     tmPat = construct_patchlist(lev,Symmetry);
// tmPat construction completes   
     Parallel::distribute(tmPat,nprocs,ingfs,fngfs,false);
//    checkPatchList(tmPat,true);
     Parallel::fill_level_data(tmPat,PatL[lev],PatL[lev-1],OldList, StateList, FutureList,tmList,Symmetry,BB);

     Parallel::KillBlocks(PatL[lev]);
     PatL[lev]->destroyList();
     PatL[lev] = tmPat;
    }
}
void cgh::read_bbox(int Symmetry,char *filename)
{
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
// read parameter from file
    {
      const int LEN = 256;   
      char pline[LEN];
      string str, sgrp, skey, sval;
      int sind1,sind2,sind3;
      ifstream inf(filename,  ifstream::in);
      if ( !inf.good() && myrank==0 ) 
      {
        cout<<"cgh::cgh: Can not open parameter file "<<filename<<" for inputing information of black holes"<<endl;
       	MPI_Abort(MPI_COMM_WORLD,1);
      }

      for( int i=1; inf.good(); i++)  
      {
        inf.getline(pline, LEN);
        str = pline;
     
        int status = misc::parse_parts(str, sgrp, skey, sval, sind1);
        if ( status == -1 ) { cout<<"error reading parameter file "<<filename<<" in line "<<i<<endl; MPI_Abort(MPI_COMM_WORLD,1);}
        else if( status == 0 ) continue;
   
	if(sgrp == "cgh" && skey == "levels") { levels = atoi(sval.c_str()); break;}
      }
      inf.close();
    }

    grids = new int[levels];
    shape = new int**[levels];
    handle = new double**[levels];
    bbox = new double**[levels];
    PatL = new MyList<Patch>*[levels];
// read parameter from file
    {
      const int LEN = 256;   
      char pline[LEN];
      string str, sgrp, skey, sval;
      int sind1,sind2,sind3;
      ifstream inf(filename,  ifstream::in);
      if ( !inf.good() && myrank==0 ) 
      {
        cout<<"cgh::cgh: Can not open parameter file "<<filename<<" for inputing information of black holes"<<endl;
       	MPI_Abort(MPI_COMM_WORLD,1);
      }

      for( int i=1; inf.good(); i++)  
      {
        inf.getline(pline, LEN);
        str = pline;
     
        int status = misc::parse_parts(str, sgrp, skey, sval, sind1,sind2,sind3);
        if ( status == -1 ) { cout<<"error reading parameter file "<<filename<<" in line "<<i<<endl; MPI_Abort(MPI_COMM_WORLD,1);}
        else if( status == 0 ) continue;
   
	if(sgrp == "cgh" && skey == "grids" && sind1 < levels) grids[sind1] = atoi(sval.c_str());
      }
      inf.close();
    }

    for(int sind1=0;sind1<levels;sind1++)
    {
      shape[sind1] = new int*[grids[sind1]];
      handle[sind1]  = new double*[grids[sind1]];
      bbox[sind1]  = new double*[grids[sind1]];
      for(int sind2=0;sind2<grids[sind1];sind2++)
      {
        shape[sind1][sind2] = new int[dim];
        handle[sind1][sind2] = new double[dim];
        bbox[sind1][sind2]  = new double[2*dim];
      }
    }
// read parameter from file
    {
      const int LEN = 256;   
      char pline[LEN];
      string str, sgrp, skey, sval;
      int sind1,sind2,sind3;
      ifstream inf(filename,  ifstream::in);
      if ( !inf.good() && myrank==0 ) 
      {
        cout<<"cgh::cgh: Can not open parameter file "<<filename<<" for inputing information of black holes"<<endl;
       	MPI_Abort(MPI_COMM_WORLD,1);
      }

      for( int i=1; inf.good(); i++)  
      {
        inf.getline(pline, LEN);
        str = pline;
     
        int status = misc::parse_parts(str, sgrp, skey, sval, sind1,sind2,sind3);

        if ( status == -1 ) { cout<<"error reading parameter file "<<filename<<" in line "<<i<<endl; MPI_Abort(MPI_COMM_WORLD,1);}
        else if( status == 0 ) continue;
   
	if(sgrp == "cgh" && sind1<levels && sind2<grids[sind1])
	{
              if ( skey == "bbox")  bbox[sind1][sind2][sind3] = atof(sval.c_str());
	 else if ( skey == "shape") shape[sind1][sind2][sind3] = atoi(sval.c_str());
	}
      }
      inf.close();
    }
// we always assume the input parameter is in cell center style    
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif  
    for(int lev=0;lev<levels;lev++)
      for(int grd=0;grd<grids[lev];grd++)
      {
	for(int i=0;i<dim;i++)
	{
  
	  shape[lev][grd][i] = shape[lev][grd][i]+1;
	}
      }
#endif
// boxes align check
    {
        double DH0[dim];
	for(int i=0;i<dim;i++)
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
           DH0[i]=(bbox[0][0][i+dim]-bbox[0][0][i])/(shape[0][0][i]-1);
#else
#ifdef Cell
	   DH0[i]=(bbox[0][0][i+dim]-bbox[0][0][i])/shape[0][0][i];
#else
#error Not define Vertex nor Cell
#endif  
#endif
    for(int lev=0;lev<levels;lev++)
      for(int grd=0;grd<grids[lev];grd++)
	 Parallel::aligncheck(bbox[0][0],bbox[lev][grd],lev,DH0,shape[lev][grd]);
    }
// print information of cgh   
    if(myrank==0)
    { 
     cout<<"cgh has levels: "<<levels<<endl;
     for(int lev=0;lev<levels;lev++)
     {
      cout<<"level #"<<lev<<" has boxes: "<<grids[lev]<<endl;
      for(int grd=0;grd<grids[lev];grd++)
      {
        cout<<"#"<<grd<<"box is"   <<"  ("<<bbox[lev][grd][0]<<":"<<bbox[lev][grd][3]
	 			   <<","<<bbox[lev][grd][1]<<":"<<bbox[lev][grd][4]
	 			   <<","<<bbox[lev][grd][2]<<":"<<bbox[lev][grd][5]
	 			   <<")."<<endl;
      }
     }
    }
}
MyList<Patch> *cgh::construct_patchlist(int lev,int Symmetry)
{
//Construct Patches
      MyList<Patch> *tmPat=0;
//construct box list
      MyList<Parallel::gridseg> *boxes=0,*gs;
      for(int grd=0;grd<grids[lev];grd++)
      {
        if(boxes)
      	{
           gs->next = new MyList<Parallel::gridseg>;
           gs = gs->next;
           gs->data = new Parallel::gridseg;
        }
        else
        {
           boxes = gs = new MyList<Parallel::gridseg>;
           gs->data = new Parallel::gridseg;
        }
        for(int i=0;i<dim;i++)
        {
           gs->data->llb[i]=bbox[lev][grd][i];
           gs->data->uub[i]=bbox[lev][grd][dim+i];
           gs->data->shape[i]=shape[lev][grd][i];
        }
        gs->data->Bg=0;
        gs->next=0;
     }
     Parallel::merge_gsl(boxes,ratio);
     Parallel::cut_gsl(boxes);
     Parallel::add_ghost_touch(boxes);

     MyList<Patch> *gp;
     gs=boxes;
     while(gs)
     {
        double tbb[2*dim];
        if(tmPat)
        {
          gp->next = new MyList<Patch>;
	  gp = gp->next;
	  for(int i=0;i<dim;i++)
	  {
	       tbb[i]=gs->data->llb[i];
	       tbb[dim+i]=gs->data->uub[i];
	  }
#ifdef WithShell   
	  gp->data = new Patch(3,gs->data->shape,tbb,lev,true,Symmetry);
#else	  
	  gp->data = new Patch(3,gs->data->shape,tbb,lev,(lev>0),Symmetry);
#endif	  
	}
	else
	{
          tmPat = gp = new MyList<Patch>;
	  for(int i=0;i<dim;i++)
	  {
	     tbb[i]=gs->data->llb[i];
	     tbb[dim+i]=gs->data->uub[i];
	  }
#ifdef WithShell   
	  gp->data = new Patch(3,gs->data->shape,tbb,lev,true,Symmetry);
#else	  
	  gp->data = new Patch(3,gs->data->shape,tbb,lev,(lev>0),Symmetry);
#endif
	}
	gp->next = 0;

	gs=gs->next;
     }

     boxes->destroyList();

     return tmPat;
}
