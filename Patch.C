//$Id: Patch.C,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <new>
using namespace std;

#include "Patch.h"
#include "Parallel.h"
#include "fmisc.h"

Patch::Patch(int DIM, int *shapei, double *bboxi,int levi,bool buflog,int Symmetry):
lev(levi)
{
  if(DIM != dim)
  {
     cout<<"dimension is not consistent in Patch construction"<<endl;
     MPI_Abort(MPI_COMM_WORLD,1);
  }
  for(int i=0;i<dim;i++)
  {
    shape[i]    = shapei[i];
    bbox[i]     = bboxi[i];
    bbox[dim+i] = bboxi[dim+i];
    lli[i] = uui[i] = 0;
    if(buflog)
    {
      double DH;
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
      DH = (bbox[dim+i]-bbox[i])/(shape[i]-1);
#else
#ifdef Cell
      DH = (bbox[dim+i]-bbox[i])/shape[i];
#else
#error Not define Vertex nor Cell
#endif  
#endif
      uui[i] = buffer_width;
      bbox[dim+i] = bbox[dim+i] + uui[i]*DH;
      shape[i] = shape[i]+uui[i];
    }
  }
      
  if(buflog)
  {
    if(DIM != 3)
    {
     cout<<"Symmetry in Patch construction only support 3 yet but dim = "<<DIM<<endl;
     MPI_Abort(MPI_COMM_WORLD,1);
    }
    double tmpb,DH;
    if(Symmetry >0 )
    {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
       DH = (bbox[5]-bbox[2])/(shape[2]-1); 
#else
#ifdef Cell
       DH = (bbox[5]-bbox[2])/shape[2]; 
#else
#error Not define Vertex nor Cell
#endif  
#endif 
       tmpb = max(0.0,bbox[2]-buffer_width*DH);    
       lli[2] = int((bbox[2]-tmpb)/DH+0.4);
       bbox[2] = bbox[2] - lli[2]*DH;
       shape[2] = shape[2]+lli[2];
       if(lli[2] < buffer_width)
       {
	  if(feq(bbox[2],0,DH/2)) lli[2] = 0;
	  else
	  {
              cout<<"Code mistake for lli[2] = "<<lli[2]<<", bbox[2] = "<<bbox[2]<<endl;
              MPI_Abort(MPI_COMM_WORLD,1);
	  }
       }
       if(Symmetry > 1)
       {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
	 DH = (bbox[3]-bbox[0])/(shape[0]-1); 
#else
#ifdef Cell
	 DH = (bbox[3]-bbox[0])/shape[0]; 
#else
#error Not define Vertex nor Cell
#endif  
#endif 
         tmpb = max(0.0,bbox[0]-buffer_width*DH);       
         lli[0] = int((bbox[0]-tmpb)/DH+0.4);
         bbox[0] = bbox[0] - lli[0]*DH;
         shape[0] = shape[0]+lli[0];
         if(lli[0] < buffer_width)
         {
	    if(feq(bbox[0],0,DH/2)) lli[0] = 0;
	    else
	    {
               cout<<"Code mistake for lli[0] = "<<lli[0]<<", bbox[0] = "<<bbox[0]<<endl;
               MPI_Abort(MPI_COMM_WORLD,1);
	    }
         }
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
	 DH = (bbox[4]-bbox[1])/(shape[1]-1); 
#else
#ifdef Cell
	 DH = (bbox[4]-bbox[1])/shape[1]; 
#else
#error Not define Vertex nor Cell
#endif  
#endif 
         tmpb = max(0.0,bbox[1]-buffer_width*DH);          
         lli[1] = int((bbox[1]-tmpb)/DH+0.4);
         bbox[1] = bbox[1] - lli[1]*DH;
         shape[1] = shape[1]+lli[1];
         if(lli[1] < buffer_width)
         {
	   if(feq(bbox[1],0,DH/2)) lli[1] = 0;
	   else
	   {
               cout<<"Code mistake for lli[1] = "<<lli[1]<<", bbox[1] = "<<bbox[1]<<endl;
               MPI_Abort(MPI_COMM_WORLD,1);
	   }
         }
       }
       else
       {
         for(int i=0;i<2;i++)
         {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
            DH = (bbox[dim+i]-bbox[i])/(shape[i]-1);
#else
#ifdef Cell
            DH = (bbox[dim+i]-bbox[i])/shape[i];
#else
#error Not define Vertex nor Cell
#endif  
#endif
	    lli[i] = buffer_width;
            bbox[i] = bbox[i] - lli[i]*DH;
            shape[i] = shape[i]+lli[i];
	 }
       }
    }
    else
    {
       for(int i=0;i<dim;i++)
       {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
          DH = (bbox[dim+i]-bbox[i])/(shape[i]-1);
#else
#ifdef Cell
          DH = (bbox[dim+i]-bbox[i])/shape[i];
#else
#error Not define Vertex nor Cell
#endif  
#endif
	  lli[i] = buffer_width;
          bbox[i] = bbox[i] - lli[i]*DH;
          shape[i] = shape[i]+lli[i];
       }
    }
  }

  blb = ble = 0;
}
Patch::~Patch()
{

}
// buflog 1: with buffer points; 0 without
void Patch::checkPatch(bool buflog)
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  if(myrank==0)
  {
    if(buflog)
    {
     cout<<"belong to level "<<lev<<endl;
     cout<<" shape: [";
     for(int i=0;i<dim;i++) 
     {
       cout<<shape[i];
       if(i<dim-1) cout<<",";
       else        cout<<"]";
     }
     cout<<" resolution: [";
     for(int i=0;i<dim;i++) 
     {
       cout<<getdX(i);
       if(i < dim-1) cout<<",";
       else          cout<<"]"<<endl;
     }
     cout<<"range:"<<"(";
     for(int i=0;i<dim;i++) 
     {
       cout<<bbox[i]<<":"<<bbox[dim+i];
       if(i<dim-1) cout<<",";
       else        cout<<")"<<endl;
     }
    }
    else
    {
     cout<<"belong to level "<<lev<<endl;
     cout<<" shape: [";
     for(int i=0;i<dim;i++) 
     {
       cout<<shape[i]-lli[i]-uui[i];
       if(i<dim-1) cout<<",";
       else        cout<<"]";
     }
     cout<<" resolution: [";
     for(int i=0;i<dim;i++) 
     {
       cout<<getdX(i);
       if(i < dim-1) cout<<",";
       else          cout<<"]"<<endl;
     }
     cout<<"range:"<<"(";
     for(int i=0;i<dim;i++) 
     {
       cout<<bbox[i]+lli[i]*getdX(i)<<":"<<bbox[dim+i]-uui[i]*getdX(i);
       if(i<dim-1) cout<<",";
       else        cout<<")"<<endl;
     }
    }
  }
}
void Patch::Interp_Points(MyList<var> *VarList,
		int NN,double **XX,
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
    for(int i=0;i<dim;i++)
    {
      pox[i] = XX[i][j];
      if(myrank == 0 && (XX[i][j] < bbox[i]+lli[i]*DH[i] || XX[i][j] > bbox[dim+i]-uui[i]*DH[i]))
      {
        cout<<"Patch::Interp_Points: point (";
	for(int k=0;k<dim;k++)
	{
	  cout<<XX[k][j];
	  if(k < dim-1) cout<<",";
	  else          cout<<") is out of current Patch."<<endl;
	}
        MPI_Abort(MPI_COMM_WORLD,1);
      }
    }

     MyList<Block> *Bp=blb;
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
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
          llb[i] = (feq(BP->bbox[i]    ,bbox[i]    ,DH[i]/2)) ? BP->bbox[i]+lli[i]*DH[i]     : BP->bbox[i]    +(ghost_width-0.5)*DH[i];
          uub[i] = (feq(BP->bbox[dim+i],bbox[dim+i],DH[i]/2)) ? BP->bbox[dim+i]-uui[i]*DH[i] : BP->bbox[dim+i]-(ghost_width-0.5)*DH[i];
#else
#ifdef Cell
          llb[i] = (feq(BP->bbox[i]    ,bbox[i]    ,DH[i]/2)) ? BP->bbox[i]+lli[i]*DH[i]     : BP->bbox[i]    +ghost_width*DH[i];
          uub[i] = (feq(BP->bbox[dim+i],bbox[dim+i],DH[i]/2)) ? BP->bbox[dim+i]-uui[i]*DH[i] : BP->bbox[dim+i]-ghost_width*DH[i];
#else
#error Not define Vertex nor Cell
#endif  
#endif 
	  if(XX[i][j]-llb[i]<-DH[i]/2 || XX[i][j]-uub[i]>DH[i]/2) {flag = false; break;}
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
//              shellf[j*num_var+k] = Parallel::global_interp(dim,BP->shape,BP->X,BP->fgfs[varl->data->sgfn],
//	  		                                    pox,ordn,varl->data->SoA,Symmetry);
              f_global_interp(BP->shape,BP->X[0],BP->X[1],BP->X[2],BP->fgfs[varl->data->sgfn],shellf[j*num_var+k],
			      pox[0],pox[1],pox[2],ordn,varl->data->SoA,Symmetry);
	      varl=varl->next;
	      k++;
	   }
	   weight[j] = 1;
	 }
	}
	if(Bp == ble) break;
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
       if(myrank == 0) cout<<"WARNING: Patch::Interp_Points meets multiple weight"<<endl;
       for(int j=0;j<num_var;j++)
	 Shellf[j+i*num_var] = Shellf[j+i*num_var]/Weight[i];
     }
     else if(Weight[i] == 0 && myrank == 0)
     {
       cout<<"ERROR: Patch::Interp_Points fails to find point (";
       for(int j=0;j<dim;j++)
       {
	  cout<<XX[j][i];
	  if(j<dim-1) cout<<",";
	  else        cout<<")";
       }
       cout<<" on Patch (";
       for(int j=0;j<dim;j++)
       {
	  cout<<bbox[j]<<"+"<<lli[j]*getdX(j);
	  if(j<dim-1) cout<<",";
	  else        cout<<")--";
       }
       cout<<"(";
       for(int j=0;j<dim;j++)
       {
	  cout<<bbox[dim+j]<<"-"<<uui[j]*getdX(j);
	  if(j<dim-1) cout<<",";
	  else        cout<<")"<<endl;
       }
#if 0
       checkBlock();
#else
  cout<<"splited domains:"<<endl;
  {
     MyList<Block> *Bp=blb;
     while(Bp)
     {
	Block *BP=Bp->data;

	for(int i=0;i<dim;i++)
	{
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
          llb[i] = (feq(BP->bbox[i]    ,bbox[i]    ,DH[i]/2)) ? BP->bbox[i]+lli[i]*DH[i]     : BP->bbox[i]    +(ghost_width-0.5)*DH[i];
          uub[i] = (feq(BP->bbox[dim+i],bbox[dim+i],DH[i]/2)) ? BP->bbox[dim+i]-uui[i]*DH[i] : BP->bbox[dim+i]-(ghost_width-0.5)*DH[i];
#else
#ifdef Cell
          llb[i] = (feq(BP->bbox[i]    ,bbox[i]    ,DH[i]/2)) ? BP->bbox[i]+lli[i]*DH[i]     : BP->bbox[i]    +ghost_width*DH[i];
          uub[i] = (feq(BP->bbox[dim+i],bbox[dim+i],DH[i]/2)) ? BP->bbox[dim+i]-uui[i]*DH[i] : BP->bbox[dim+i]-ghost_width*DH[i];
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
	if(Bp == ble) break;
	Bp=Bp->next;
     }
  }
#endif       
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
void Patch::checkBlock()
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  if(myrank==0)
  {
     MyList<Block> *BP=blb;
     while(BP)
     {
	BP->data->checkBlock();
	if(BP==ble) break;
	BP=BP->next;
     }
  }
}
double Patch::getdX(int dir)
{
    if(dir < 0 || dir >= dim)
    {
        cout<<"Patch::getdX: error input dir = "<<dir<<", this Patch has direction (0,"<<dim-1<<")"<<endl;
        MPI_Abort(MPI_COMM_WORLD,1);
    }	    
    double h;
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif   
    if(shape[dir] == 1)
    {
        cout<<"Patch::getdX: for direction "<<dir<<", this Patch has only one point. Can not determine dX for vertex center grid."<<endl;
        MPI_Abort(MPI_COMM_WORLD,1);
    } 
    h = (bbox[dim+dir]-bbox[dir])/(shape[dir]-1);
#else
#ifdef Cell
    h = (bbox[dim+dir]-bbox[dir])/shape[dir];
#else
#error Not define Vertex nor Cell
#endif  
#endif
    return h;
}
