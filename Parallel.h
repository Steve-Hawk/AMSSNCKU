//$Id: Parallel.h,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
#ifndef PARALLEL_H
#define PARALLEL_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <new>
using namespace std;


#include "var.h"
#include "Patch.h"
#include "Block.h"
#include "MyList.h"
#include "microdef.h" //need dim; ghost_width; CONTRACT
namespace Parallel
{
  struct gridseg
  {
    double llb[dim];
    double uub[dim];
    int shape[dim];
    Block *Bg;
  };

  int partition3(int *nxyz,int split_size,int *min_width,int cpusize,int *shape);
  MyList<Block> *distribute(MyList<Patch> *PatchLIST, int cpusize, int ingfsi, int fngfs,bool periodic); //produce corresponding Blocks
  void KillBlocks(MyList<Patch> *PatchLIST);

  void setfunction(MyList<Block> *BlL, var *vn, double func(double x, double y, double z));
  void setfunction(int rank,MyList<Block> *BlL, var *vn, double func(double x, double y, double z));
  void writefile(double time,int nx,int ny,int nz,double xmin,double xmax,double ymin,double ymax,
	       double zmin,double zmax,char *filename,double *data_out);
  void getarrayindex(int DIM,int *shape,int *index,int n);
  int getarraylocation(int DIM,int *shape,int *index);
  void copy(int DIM,double *llbout,double *uubout,int *Dshape,double *DD,double *llbin,double *uubin,
	    int *shape,double *datain,double *llb,double *uub);
  void Dump_CPU_Data(MyList<Block> *BlL,MyList<var> *DumpList,char *tag,double time,double dT);
  void Dump_Data(MyList<Patch> *PL,MyList<var> *DumpList,char *tag,double time,double dT);
  void Dump_Data(Patch *PP,MyList<var> *DumpList,char *tag,double time,double dT,int grd);
  void Dump_Data0(Patch *PP,MyList<var> *DumpList,char *tag,double time,double dT);
  double global_interp(int DIM,int *ext,double **CoX,double *datain,
		       double *poX, int ordn, double *SoA, int Symmetry);
  double global_interp(int DIM,int *ext,double **CoX,double *datain,
		       double *poX, int ordn);
  double Lagrangian_Int(double x,int npts,double *xpts,double *funcvals);
  double LagrangePoly(double x,int pt,int npts, double *xpts);
  MyList<gridseg> * build_complete_gsl(Patch *Pat);
  MyList<gridseg> * build_complete_gsl(MyList<Patch> *PatL);
  MyList<gridseg> * build_complete_gsl_virtual(MyList<Patch> *PatL);
  MyList<gridseg> * build_complete_gsl_virtual2(MyList<Patch> *PatL); // - buffer
  MyList<gridseg> * build_owned_gsl0(Patch *Pat,int rank_in);  // - ghost without extension, special for Sync usage
  MyList<gridseg> * build_owned_gsl1(Patch *Pat,int rank_in);  // - ghost, similar to build_owned_gsl0 but extend one point on left side for vertex grid
  MyList<gridseg> * build_owned_gsl2(Patch *Pat,int rank_in);  // - buffer - ghost
  MyList<gridseg> * build_owned_gsl3(Patch *Pat,int rank_in,int Symmetry);  // - ghost - BD ghost
  MyList<gridseg> * build_owned_gsl4(Patch *Pat,int rank_in,int Symmetry);  // - buffer - ghost - BD ghost
  MyList<gridseg> * build_owned_gsl5(Patch *Pat,int rank_in);  // similar to build_owned_gsl2 but no extension
  MyList<gridseg> * build_owned_gsl(MyList<Patch> *PatL,int rank_in,int type,int Symmetry);
  void build_gstl(MyList<gridseg> *srci,MyList<gridseg> *dsti,MyList<gridseg> **out_src,MyList<gridseg> **out_dst);
  int data_packer(double *data,MyList<gridseg> *src,MyList<gridseg> *dst, int rank_in, int dir,
                      MyList<var> *VarLists,MyList<var> *VarListd,int Symmetry);
  void transfer(MyList<gridseg> **src,MyList<gridseg> **dst,
		MyList<var> *VarList1 /* source */,MyList<var> *VarList2 /*target */,
		int Symmetry);
  void Sync(Patch *Pat,MyList<var> *VarList,int Symmetry);
  void Sync(MyList<Patch> *PatL,MyList<var> *VarList,int Symmetry);
  void OutBdLow2Hi(Patch *Patc,Patch *Patf,
	           MyList<var> *VarList1 /* source */,MyList<var> *VarList2 /* target */,
		   int Symmetry);
  void OutBdLow2Hi(MyList<Patch> *PatcL,MyList<Patch> *PatfL,
	           MyList<var> *VarList1 /* source */,MyList<var> *VarList2 /* target */,
		   int Symmetry);
  void  Prolong(Patch *Patc,Patch *Patf,
	        MyList<var> *VarList1 /* source */,MyList<var> *VarList2 /* target */,
		int Symmetry);
  void  Prolongint(Patch *Patc,Patch *Patf,
	           MyList<var> *VarList1 /* source */,MyList<var> *VarList2 /* target */,
		   int Symmetry);
  void Restrict(MyList<Patch> *PatcL,MyList<Patch> *PatfL,
	        MyList<var> *VarList1 /* source */,MyList<var> *VarList2 /* target */,
		int Symmetry);
  MyList<Parallel::gridseg> *build_PhysBD_gsl(Patch *Pat);
  MyList<Parallel::gridseg> *build_ghost_gsl(MyList<Patch> *PatL);
  MyList<Parallel::gridseg> *build_ghost_gsl(Patch *Pat);
  MyList<Parallel::gridseg> *build_buffer_gsl(Patch *Pat);
  MyList<Parallel::gridseg> *build_buffer_gsl(MyList<Patch> *PatL);
  MyList<Parallel::gridseg> *gsl_subtract(MyList<Parallel::gridseg> *A,MyList<Parallel::gridseg> *B);
  MyList<Parallel::gridseg> *gs_subtract(MyList<Parallel::gridseg> *A,MyList<Parallel::gridseg> *B);
  MyList<Parallel::gridseg> *gsl_and(MyList<Parallel::gridseg> *A,MyList<Parallel::gridseg> *B);
  MyList<Parallel::gridseg> *gs_and(MyList<Parallel::gridseg> *A,MyList<Parallel::gridseg> *B);
  MyList<Parallel::gridseg> *clone_gsl(MyList<Parallel::gridseg> *p,bool first_only);
  MyList<Parallel::gridseg> *build_bulk_gsl(Patch *Pat);  // similar to build_owned_gsl0 but does not care rank issue
  MyList<Parallel::gridseg> *build_bulk_gsl(Block *bp, Patch *Pat);
  void build_PhysBD_gstl(Patch *Pat,MyList<Parallel::gridseg> *srci,MyList<Parallel::gridseg> *dsti,
		         MyList<Parallel::gridseg> **out_src,MyList<Parallel::gridseg> **out_dst);
  void PeriodicBD(Patch *Pat,MyList<var> *VarList,int Symmetry);
  double L2Norm(Patch *Pat, var* vf);
  void checkgsl(MyList<Parallel::gridseg>* pp,bool first_only);
  void checkvarl(MyList<var>* pp,bool first_only);
  MyList<Parallel::gridseg> *divide_gsl(MyList<Parallel::gridseg> *p, Patch *Pat);
  MyList<Parallel::gridseg> *divide_gs(MyList<Parallel::gridseg> *p, Patch *Pat);
  void prepare_inter_time_level(Patch *Pat,
		MyList<var> *VarList1 /* source (t+dt) */,MyList<var> *VarList2 /* source (t) */,            		 
		MyList<var> *VarList3 /* target (t+a*dt) */, int tindex);
  void prepare_inter_time_level(Patch *Pat, 
                MyList<var> *VarList1 /* source (t+dt) */,MyList<var> *VarList2 /* source (t) */,
		MyList<var> *VarList3 /* source (t-dt) */,MyList<var> *VarList4 /* target (t+a*dt) */, int tindex);
  void prepare_inter_time_level(MyList<Patch> *PatL,
		MyList<var> *VarList1 /* source (t+dt) */,MyList<var> *VarList2 /* source (t) */,            		 
		MyList<var> *VarList3 /* target (t+a*dt) */, int tindex);
  void prepare_inter_time_level(MyList<Patch> *Pat, 
                MyList<var> *VarList1 /* source (t+dt) */,MyList<var> *VarList2 /* source (t) */,
		MyList<var> *VarList3 /* source (t-dt) */,MyList<var> *VarList4 /* target (t+a*dt) */, int tindex);
  void merge_gsl(MyList<gridseg> *&A,const double ratio);
  bool merge_gs(MyList<gridseg> *D,MyList<gridseg> *B,MyList<gridseg> *&C,const double ratio);
// Add ghost region to tangent plane
// we assume the grids have the same resolution
  void add_ghost_touch(MyList<gridseg> *&A);
  void cut_gsl(MyList<gridseg> *&A);
  bool cut_gs(MyList<gridseg> *D,MyList<gridseg> *B,MyList<gridseg> *&C);
  MyList<Parallel::gridseg> *gs_subtract_virtual(MyList<Parallel::gridseg> *A,MyList<Parallel::gridseg> *B);
  void fill_level_data(MyList<Patch> *PatLd,MyList<Patch> *PatLs,MyList<Patch> *PatcL,
		               MyList<var> *OldList, MyList<var> *StateList, MyList<var> *FutureList,
			       MyList<var> *tmList,int Symmetry,bool BB);
  bool PatList_Interp_Points(MyList<Patch> *PatL,MyList<var> *VarList,
		                     int NN,double **XX,
		                     double *Shellf,int Symmetry);
  void aligncheck(double *bbox0,double *bboxl,int lev,double *DH0,int *shape);
  bool point_locat_gsl(double *pox,MyList<Parallel::gridseg> *gsl);
  void checkpatchlist(MyList<Patch> *PatL,bool buflog);
}
#endif    /*PARALLEL_H */
