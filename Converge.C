//$Id: Converge.C,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
//for vertex center grid
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

#define overghost ((ghost_width+1)/2+ghost_width)
const int ghost=ghost_width;
const double basict=0.1;
const int buffer=buffer_width;

double* readdata(int &nx,int &ny,int &nz,
		 double &xmin,double &xmax,
		 double &ymin,double &ymax,
		 double &zmin,double &zmax,
		 char *fname)
{
  ifstream infile1;
  infile1.open(fname);
  if(!infile1){
    cerr << "\a Can't open " << fname << " for input." << endl;
    exit(1);
  }

  double time;
  infile1.seekg(0, ios::beg);
  infile1.read((char *) &time, sizeof(double));
  infile1.read((char *) &nx, sizeof(int));
  infile1.read((char *) &ny, sizeof(int));
  infile1.read((char *) &nz, sizeof(int));
  infile1.read((char *) &xmin, sizeof(double));
  infile1.read((char *) &xmax, sizeof(double));
  infile1.read((char *) &ymin, sizeof(double));
  infile1.read((char *) &ymax, sizeof(double));
  infile1.read((char *) &zmin, sizeof(double));
  infile1.read((char *) &zmax, sizeof(double));

  double *data;
  data = new double[nx*ny*nz];
  infile1.read((char *) data, nx*ny*nz*sizeof(double));
  infile1.close();

  return data;

}
void writefile(double time,int nx,int ny,int nz,double xmin,double xmax,double ymin,double ymax,
	       double zmin,double zmax,char *filename,double *data_out)
{   
    ofstream outfile;
    outfile.open(filename,ios::out | ios::trunc);
    if(!outfile)
    {
      cout << "Can't open " << filename << " for output." << endl;
      exit(0);
    }
    outfile.write((char *) &time, sizeof(double));
    outfile.write((char *) &nx,sizeof(int));
    outfile.write((char *) &ny,sizeof(int));
    outfile.write((char *) &nz,sizeof(int));
    outfile.write((char *) &xmin,sizeof(double));
    outfile.write((char *) &xmax,sizeof(double));
    outfile.write((char *) &ymin,sizeof(double));
    outfile.write((char *) &ymax,sizeof(double));
    outfile.write((char *) &zmin,sizeof(double));
    outfile.write((char *) &zmax,sizeof(double));
    outfile.write((char *) data_out,nx*ny*nz*sizeof(double));
    outfile.close();
}
// here index excludes ghost
void coarsemfine(int cx, int cy, int cz, double* datac,
                 int fx, int fy, int fz, double* dataf,
		 int bx,int by,int bz,
		 int ux,int uy,int uz)
{
  for(int k=bz;k<cz-ux;k++)
   for(int j=by;j<cy-uy;j++)
    for(int i=bx;i<cx-uz;i++)
    {
      datac[i+j*cx+k*cx*cy] = datac[i+j*cx+k*cx*cy] - dataf[(bx+2*(i-bx))+(by+2*(j-by))*fx+(bz+2*(k-bz))*fx*fy];
    }
}
double pow2(double x)
{
  return x*x;
}
double L2norm(int cx,int cy,int cz, double dX, double dY, double dZ, double *gzz,int rho,
	      int bx,int by,int bz,
	      int ux,int uy,int uz)
{
  double normhelper=0;
  for(int k=bz;k<cz-ux;k=k+rho)
   for(int j=by;j<cy-uy;j=j+rho)
    for(int i=bx;i<cx-uz;i=i+rho)
    {
      normhelper += pow2(gzz[i+j*cx+k*cx*cy]);
    }

  return sqrt(normhelper*dX*dY*dZ);
}
double onenormnew(double time,int rhoc,int rhof,char* dirc,char* dirf,int Symmetry)	
{
  static bool fit=true;
  int cx,cy,cz,fx,fy,fz;
  double dTc=basict/rhoc,dTf=basict/rhof;
  double xmin,ymin,zmin,xmax,ymax,zmax;

  double dX, dY, dZ;

  double ret=0;
  double *FFc,*FFf;
  char filename[100];
//00
  dX=0.2;
  dY=0.2;
  dZ=0.2;
  sprintf(filename,"./%s/Lev00-00_Sphi0_%05d.bin",dirc,int(time/dTc+0.5));
  FFc = readdata(cx,cy,cz,xmin,xmax,ymin,ymax,zmin,zmax,filename);
  sprintf(filename,"./%s/Lev00-00_Sphi0_%05d.bin",dirf,int(time/dTf+0.5)); 
  FFf = readdata(fx,fy,fz,xmin,xmax,ymin,ymax,zmin,zmax,filename);
  if(fit)
  {
    cout<<"cx,cy,cz = "<<cx<<","<<cy<<","<<cz<<endl;
    cout<<"fx,fy,fz = "<<fx<<","<<fy<<","<<fz<<endl;
  }
  if(Symmetry == 0)       coarsemfine(cx, cy, cz, FFc, fx, fy, fz, FFf,buffer,buffer,buffer,buffer,buffer,buffer);
  else if(Symmetry == 1)  coarsemfine(cx, cy, cz, FFc, fx, fy, fz, FFf,buffer,buffer,0,buffer,buffer,buffer);
  else if(Symmetry == 2)  coarsemfine(cx, cy, cz, FFc, fx, fy, fz, FFf,0,0,0,buffer,buffer,buffer);
  sprintf(filename,"./plot/Lev00-00_Sphi0_%05d.bin",int(time/dTc+0.5));
  writefile(time,cx, cy, cz,xmin,xmax,ymin,ymax,zmin,zmax,filename,FFc);
  if(Symmetry == 0)       ret += L2norm(cx, cy, cz,dX,dY,dZ,FFc,rhoc,buffer,buffer,buffer,buffer,buffer,buffer);
  else if(Symmetry == 1)  ret += L2norm(cx, cy, cz,dX,dY,dZ,FFc,rhoc,buffer,buffer,0,buffer,buffer,buffer);
  else if(Symmetry == 2)  ret += L2norm(cx, cy, cz,dX,dY,dZ,FFc,rhoc,0,0,0,buffer,buffer,buffer);
  delete[] FFc; delete[] FFf;
/*
//01
  dX=0.1;
  dY=0.1;
  dZ=0.1;
  sprintf(filename,"./%s/Lev01-00_Sphi0_%05d.bin",dirc,int(time/dTc+0.5));
  FFc = readdata(cx,cy,cz,xmin,xmax,ymin,ymax,zmin,zmax,filename);
  sprintf(filename,"./%s/Lev01-00_Sphi0_%05d.bin",dirf,int(time/dTf+0.5)); 
  FFf = readdata(fx,fy,fz,xmin,xmax,ymin,ymax,zmin,zmax,filename);
  if(fit)
  {
    cout<<"cx,cy,cz = "<<cx<<","<<cy<<","<<cz<<endl;
    cout<<"fx,fy,fz = "<<fx<<","<<fy<<","<<fz<<endl;
  }
  if(Symmetry == 0)       coarsemfine(cx, cy, cz, FFc, fx, fy, fz, FFf,buffer,buffer,buffer,buffer,buffer,buffer);
  else if(Symmetry == 1)  coarsemfine(cx, cy, cz, FFc, fx, fy, fz, FFf,buffer,buffer,0,buffer,buffer,buffer);
  else if(Symmetry == 2)  coarsemfine(cx, cy, cz, FFc, fx, fy, fz, FFf,0,0,0,buffer,buffer,buffer);
  sprintf(filename,"./plot/Lev01-00_Sphi0_%05d.bin",int(time/dTc+0.5));
  writefile(time,cx, cy, cz,xmin,xmax,ymin,ymax,zmin,zmax,filename,FFc);
  if(Symmetry == 0)       ret += L2norm(cx, cy, cz,dX,dY,dZ,FFc,rhoc,buffer,buffer,buffer,buffer,buffer,buffer);
  else if(Symmetry == 1)  ret += L2norm(cx, cy, cz,dX,dY,dZ,FFc,rhoc,buffer,buffer,0,buffer,buffer,buffer);
  else if(Symmetry == 2)  ret += L2norm(cx, cy, cz,dX,dY,dZ,FFc,rhoc,0,0,0,buffer,buffer,buffer);
  delete[] FFc; delete[] FFf;
//zp
  dX=0.0374;
  dY=0.0374;
  dZ=0.2;
  sprintf(filename,"./%s/LevSH-zp_Sphi0_%05d.bin",dirc,int(time/dTc+0.5));
  FFc = readdata(cx,cy,cz,xmin,xmax,ymin,ymax,zmin,zmax,filename);
  sprintf(filename,"./%s/LevSH-zp_Sphi0_%05d.bin",dirf,int(time/dTf+0.5)); 
  FFf = readdata(fx,fy,fz,xmin,xmax,ymin,ymax,zmin,zmax,filename);
  if(fit)
  {
    cout<<"cx,cy,cz = "<<cx<<","<<cy<<","<<cz<<endl;
    cout<<"fx,fy,fz = "<<fx<<","<<fy<<","<<fz<<endl;
  }
  coarsemfine(cx, cy, cz, FFc, fx, fy, fz, FFf,overghost,overghost,buffer,overghost,overghost,0);
  sprintf(filename,"./plot/LevSH-zp_Sphi0_%05d.bin",int(time/dTc+0.5));
  writefile(time,cx, cy, cz,xmin,xmax,ymin,ymax,zmin,zmax,filename,FFc);
  ret += L2norm(cx, cy, cz,dX,dY,dZ,FFc,rhoc,overghost,overghost,buffer,overghost,overghost,0);
  delete[] FFc; delete[] FFf;
//xp
  sprintf(filename,"./%s/LevSH-xp_Sphi0_%05d.bin",dirc,int(time/dTc+0.5));
  FFc = readdata(cx,cy,cz,xmin,xmax,ymin,ymax,zmin,zmax,filename);
  sprintf(filename,"./%s/LevSH-xp_Sphi0_%05d.bin",dirf,int(time/dTf+0.5)); 
  FFf = readdata(fx,fy,fz,xmin,xmax,ymin,ymax,zmin,zmax,filename);
  if(fit)
  {
    cout<<"cx,cy,cz = "<<cx<<","<<cy<<","<<cz<<endl;
    cout<<"fx,fy,fz = "<<fx<<","<<fy<<","<<fz<<endl;
  }
  coarsemfine(cx, cy, cz, FFc, fx, fy, fz, FFf,overghost,0,buffer,overghost,overghost,0);
  sprintf(filename,"./plot/LevSH-xp_Sphi0_%05d.bin",int(time/dTc+0.5));
  writefile(time,cx, cy, cz,xmin,xmax,ymin,ymax,zmin,zmax,filename,FFc);
  ret += L2norm(cx, cy, cz,dX,dY,dZ,FFc,rhoc,overghost,0,buffer,overghost,overghost,0);
  delete[] FFc; delete[] FFf;
//xm
  sprintf(filename,"./%s/LevSH-xm_Sphi0_%05d.bin",dirc,int(time/dTc+0.5));
  FFc = readdata(cx,cy,cz,xmin,xmax,ymin,ymax,zmin,zmax,filename);
  sprintf(filename,"./%s/LevSH-xm_Sphi0_%05d.bin",dirf,int(time/dTf+0.5)); 
  FFf = readdata(fx,fy,fz,xmin,xmax,ymin,ymax,zmin,zmax,filename);
  if(fit)
  {
    cout<<"cx,cy,cz = "<<cx<<","<<cy<<","<<cz<<endl;
    cout<<"fx,fy,fz = "<<fx<<","<<fy<<","<<fz<<endl;
  }
  coarsemfine(cx, cy, cz, FFc, fx, fy, fz, FFf,overghost,overghost,buffer,overghost,0,0);
  sprintf(filename,"./plot/LevSH-xm_Sphi0_%05d.bin",int(time/dTc+0.5));
  writefile(time,cx, cy, cz,xmin,xmax,ymin,ymax,zmin,zmax,filename,FFc);
  ret += L2norm(cx, cy, cz,dX,dY,dZ,FFc,rhoc,overghost,overghost,buffer,overghost,0,0);
  delete[] FFc; delete[] FFf;
//yp
  sprintf(filename,"./%s/LevSH-yp_Sphi0_%05d.bin",dirc,int(time/dTc+0.5));
  FFc = readdata(cx,cy,cz,xmin,xmax,ymin,ymax,zmin,zmax,filename);
  sprintf(filename,"./%s/LevSH-yp_Sphi0_%05d.bin",dirf,int(time/dTf+0.5)); 
  FFf = readdata(fx,fy,fz,xmin,xmax,ymin,ymax,zmin,zmax,filename);
  if(fit)
  {
    cout<<"cx,cy,cz = "<<cx<<","<<cy<<","<<cz<<endl;
    cout<<"fx,fy,fz = "<<fx<<","<<fy<<","<<fz<<endl;
  }
  coarsemfine(cx, cy, cz, FFc, fx, fy, fz, FFf,overghost,0,buffer,overghost,overghost,0);
  sprintf(filename,"./plot/LevSH-yp_Sphi0_%05d.bin",int(time/dTc+0.5));
  writefile(time,cx, cy, cz,xmin,xmax,ymin,ymax,zmin,zmax,filename,FFc);
  ret += L2norm(cx, cy, cz,dX,dY,dZ,FFc,rhoc,overghost,0,buffer,overghost,overghost,0);
  delete[] FFc; delete[] FFf;
//ym
  sprintf(filename,"./%s/LevSH-ym_Sphi0_%05d.bin",dirc,int(time/dTc+0.5));
  FFc = readdata(cx,cy,cz,xmin,xmax,ymin,ymax,zmin,zmax,filename);
  sprintf(filename,"./%s/LevSH-ym_Sphi0_%05d.bin",dirf,int(time/dTf+0.5)); 
  FFf = readdata(fx,fy,fz,xmin,xmax,ymin,ymax,zmin,zmax,filename);
  if(fit)
  {
    cout<<"cx,cy,cz = "<<cx<<","<<cy<<","<<cz<<endl;
    cout<<"fx,fy,fz = "<<fx<<","<<fy<<","<<fz<<endl;
  }
  coarsemfine(cx, cy, cz, FFc, fx, fy, fz, FFf,overghost,overghost,buffer,overghost,0,0);
  sprintf(filename,"./plot/LevSH-ym_Sphi0_%05d.bin",int(time/dTc+0.5));
  writefile(time,cx, cy, cz,xmin,xmax,ymin,ymax,zmin,zmax,filename,FFc);
  ret += L2norm(cx, cy, cz,dX,dY,dZ,FFc,rhoc,overghost,overghost,buffer,overghost,0,0);
  delete[] FFc; delete[] FFf;
*/

  if(fit) fit = false;

  return ret;
}
int main(int argc, char *argv[])
{
  if (argc < 5) 
  { 
    cout << "\aUsage: Converge rhoc dirc dirf Symmetry\n " 
	 << endl;
    exit(1);
  } 

    int rhoc, rhof, Symmetry;
    char dirc[50],dirf[50];
    rhoc=atoi(argv[1]);
    rhof=2*rhoc;
    sprintf(dirc,argv[2]);
    sprintf(dirf,argv[3]);
    Symmetry=atoi(argv[4]);
    if(Symmetry != 1)
    {
	 cout<<"at present we only suport equatorial symmetry!"<<endl;
	 exit(0);
    }

    char filename[100];
    sprintf(filename,"l2norm%05d.dat",rhoc);
    ofstream outfile;
    outfile.open(filename);

    for(double time=0.2;time<=100;time=time+0.2)
    {
      cout<<time<<"..."<<endl;
      outfile<<time<<" "<<onenormnew(time, rhoc, rhof,dirc,dirf,Symmetry)<<endl;
    }

    outfile.close();
}
