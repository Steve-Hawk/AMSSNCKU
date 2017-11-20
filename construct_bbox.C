//$Id: construct_bbox.C,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
//-----------------------------------------------------------------------
// Read binary files and do fancy things with them...
//-----------------------------------------------------------------------
#ifdef newc
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>
using namespace std;
#else
#include <math.h>
#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fstream.h>
#endif

#include "microdef.fh"
#include "microdef.h"
#include "misc.h"

int main(int argc, char *argv[])
{
     int Symmetry=1;
     int levels,movls;	
     char filename[100];
     strcpy(filename,"input.par");
     int *grids; 
     double ***bbox;
     int ***shape;
     int BH_num;
     double **BHPorg;
// read parameter from file
    {
      const int LEN = 256;   
      char pline[LEN];
      string str, sgrp, skey, sval;
      int sind1,sind2,sind3;
      ifstream inf(filename,  ifstream::in);
      if ( !inf.good() ) 
      {
        cout<<"cgh::cgh: Can not open parameter file "<<filename<<" for inputing information of black holes"<<endl;
       	exit(0);
      }

      for( int i=1; inf.good(); i++)  
      {
        inf.getline(pline, LEN);
        str = pline;
     
        int status = misc::parse_parts(str, sgrp, skey, sval, sind1);
        if ( status == -1 ) { cout<<"error reading parameter file "<<filename<<" in line "<<i<<endl; exit(0);}
        else if( status == 0 ) continue;
   
	if(sgrp == "cgh" && skey == "levels") levels = atoi(sval.c_str());
	else if(sgrp == "cgh" && skey == "moving levels start from") movls = atoi(sval.c_str()); 
	else if(sgrp == "BSSN" && skey == "BH_num") BH_num = atoi(sval.c_str());
      }
      inf.close();
    }

    movls = min(movls,levels);
    movls = max(0,movls);

    BHPorg = new double*[BH_num];
    for(int i=0;i<BH_num;i++) BHPorg[i]=new double[dim];

    grids = new int[levels];
    shape = new int**[levels];
    bbox = new double**[levels];
// read parameter from file
    {
      const int LEN = 256;   
      char pline[LEN];
      string str, sgrp, skey, sval;
      int sind1,sind2,sind3;
      ifstream inf(filename,  ifstream::in);
      if ( !inf.good() ) 
      {
        cout<<"cgh::cgh: Can not open parameter file "<<filename<<" for inputing information of black holes"<<endl;
       	exit(0);
      }

      for( int i=1; inf.good(); i++)  
      {
        inf.getline(pline, LEN);
        str = pline;
     
        int status = misc::parse_parts(str, sgrp, skey, sval, sind1,sind2,sind3);
        if ( status == -1 ) { cout<<"error reading parameter file "<<filename<<" in line "<<i<<endl; exit(0);}
        else if( status == 0 ) continue;
   
	if(sgrp == "cgh" && skey == "grids" && sind1 < levels) grids[sind1] = atoi(sval.c_str());	
	else if(sgrp == "BSSN" && sind1 < BH_num)
	{
	 if ( skey == "Porgx")  BHPorg[sind1][0] = atof(sval.c_str());
	 else if ( skey == "Porgy")  BHPorg[sind1][1] = atof(sval.c_str());
	 else if ( skey == "Porgz")  BHPorg[sind1][2]= atof(sval.c_str());
	}
      }
      inf.close();
    }

    for(int sind1=0;sind1<levels;sind1++)
    {
      shape[sind1] = new int*[grids[sind1]];
      bbox[sind1]  = new double*[grids[sind1]];
      for(int sind2=0;sind2<grids[sind1];sind2++)
      {
        shape[sind1][sind2] = new int[dim];
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
      if ( !inf.good() ) 
      {
        cout<<"cgh::cgh: Can not open parameter file "<<filename<<" for inputing information of black holes"<<endl;
       	exit(0);
      }

      for( int i=1; inf.good(); i++)  
      {
        inf.getline(pline, LEN);
        str = pline;
     
        int status = misc::parse_parts(str, sgrp, skey, sval, sind1,sind2,sind3);

        if ( status == -1 ) { cout<<"error reading parameter file "<<filename<<" in line "<<i<<endl; exit(0);}
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
// we assume the number of boexes for moving levels equals to the number of black holes
    
    double DH[dim];
    for(int i=0;i<dim;i++)
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
      DH[i] = (bbox[0][0][i+dim]-bbox[0][0][i])/(shape[0][0][i]-1);
#else
#ifdef Cell
      DH[i] = (bbox[0][0][i+dim]-bbox[0][0][i])/shape[0][0][i];
#else
#error Not define Vertex nor Cell
#endif  
#endif
    for(int lev=movls;lev<levels;lev++)
    {
        if(BH_num == grids[lev])
	{
          for(int grd=0;grd<grids[lev];grd++)
	  {
             for(int i=0;i<dim;i++)
	     {
		 double th=DH[i]*pow(0.5,lev);
		 double rr;
		 rr=BHPorg[grd][i]-bbox[0][0][i];
		 if(rr<0) cout<<"WARNING: black hole is out of the base level"<<endl;
		 rr=bbox[0][0][i]+int(rr/th+0.4)*th;
                 bbox[lev][grd][i] = rr-shape[lev][grd][i]/2*th;
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
                 bbox[lev][grd][i+dim] = bbox[lev][grd][i]+(shape[lev][grd][i]-1)*th;
#else
#ifdef Cell
                 bbox[lev][grd][i+dim] = bbox[lev][grd][i]+shape[lev][grd][i]*th;
#else
#error Not define Vertex nor Cell
#endif  
#endif 
	     }
	     if(Symmetry>0)
	     {
		 double th=DH[2]*pow(0.5,lev);
		 double rr;
		 rr=BHPorg[grd][2]-bbox[0][0][2];
		 if(rr<0) cout<<"WARNING: black hole is out of the base level"<<endl;
		 rr=bbox[0][0][2]+int(rr/th+0.4)*th;
                 bbox[lev][grd][2] = max(0,rr-shape[lev][grd][2]/2*th);
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
                 bbox[lev][grd][2+dim] = bbox[lev][grd][2]+(shape[lev][grd][2]-1)*th;
#else
#ifdef Cell
                 bbox[lev][grd][2+dim] = bbox[lev][grd][2]+shape[lev][grd][2]*th;
#else
#error Not define Vertex nor Cell
#endif  
#endif 
	     }
	  }
	}
	else 
	{
	   cout<<"cgh::construct_bbox does not support BH_num = "<<BH_num<<" grids = "<<grids[lev]<<" for level#"<<lev<<endl;
           exit(0);
	}
    }

// print information of cgh   
     cout<<"cgh has levels: "<<levels<<endl;
     cout.setf(ios::scientific,ios::floatfield);
     cout.precision( 16 );
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

  for(int i=0;i<BH_num;i++) delete[] BHPorg[i];
  delete BHPorg;

  delete[] grids;
  for(int lev=0;lev<levels;lev++)
  {
    for(int grd=0;grd<grids[lev];grd++)
    {
       delete[] bbox[lev][grd];
       delete[] shape[lev][grd];
    }
    delete[] bbox[lev];
    delete[] shape[lev];
  }
  delete[] bbox;
  delete[] shape;
}
