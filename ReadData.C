//$Id: ReadData.C,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
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
#include "derivatives.h"

int main(int argc, char *argv[])
{
  //
  // USE:  ReadData flag file1 [ file2 ]
  //
  // where: - flag can be X,Y,Z,XY,XZ,YZ,C
  //
  void set_fname(char *fname);
    
  int optype=1;   //1 data itself; 2 first order derivative; 3 second order derivative
  int Symmetry;
  double SYM1,SYM2,SYM3;
  for(int i = 1; i < argc; i++)
  {
    if(!(strncmp("-h",argv[i],2)))
    { 
     cout << "\aUsage: ReadData flag binaryfile1 [ binaryfile2 ] [-d1 SYM1 SYM2 SYM3 Symmetry] [-d2 SYM1 SYM2 SYM3 Symmetry]\n " 
 	  << "   where: - flag can be X,Y,Z,XY,XZ,YZ,C\n"            //C means center
	  <<"    -d means first order derivative\n"
	  <<"    -dd means second order derivative"
	  << endl;
     exit(1);
    }
    else if(!(strncmp("-d1",argv[i],3)))
    {
      optype=2;
      SYM1 = atof(argv[i+1]);
      SYM2 = atof(argv[i+2]);
      SYM3 = atof(argv[i+3]);
      Symmetry = atoi(argv[i+4]);
    }
    else if(!(strncmp("-d2",argv[i],3)))
    {
      optype=3;
      SYM1 = atof(argv[i+1]);
      SYM2 = atof(argv[i+2]);
      SYM3 = atof(argv[i+3]);
      Symmetry = atoi(argv[i+4]);
    }
  } 

  ifstream infile1;
  infile1.open(argv[2]);
  if(!infile1){
    cerr << "\a Can't open " << argv[2] << " for input." << endl;
    exit(1);
  }

  /* read properties of the binary file */
  double time;
  int nx, ny, nz;
  double xmin, xmax, ymin, ymax, zmin, zmax;
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

  /* get rid of any 4 character suffix */
  set_fname(argv[2]);

  /* sanity check */
  if (nx != ny || nx != nz) {
    cout << "\n" << endl;
    cout << " nx, ny and nz do not agree! Using a symmetry?... "; 
    cout << "\n" << endl;
  }

  cout << "\n Reading file : " << argv[2] << endl;
  cout << "\n  Time        : " << time << endl;
  cout << "  Dimensions  : " << setw(16) << nx   << setw(16) << ny   << setw(16) << nz << endl;
  cout << "   xmin, xmax : " << setw(16) << xmin << setw(16) << xmax << endl;
  cout << "   ymin, ymax : " << setw(16) << ymin << setw(16) << ymax << endl;
  cout << "   zmin, zmax : " << setw(16) << zmin << setw(16) << zmax << endl;
  cout << "\n";

  double *data;
  data = new double[nx*ny*nz];
  int i = 0, j = 0, k = 0; 
  infile1.read((char *) data, nx*ny*nz*sizeof(double));
  infile1.close();
  //
  //
  // if second file given, open second file and subtract from first one!
  //
  //
  if (argc == 4) {
    infile1.open(argv[3]);
    if(!infile1){
      cerr << "\a Can't open " << argv[3] << " for input." << endl;
      exit(1);
    }
    double *indata;
    indata = new double[nx*ny*nz];
    // read in header 
    infile1.seekg(0, ios::beg);
    int nxin, nyin, nzin;
    infile1.read((char *) &time, sizeof(double));
    infile1.read((char *) &nxin, sizeof(int));
    infile1.read((char *) &nyin, sizeof(int));
    infile1.read((char *) &nzin, sizeof(int));
    infile1.read((char *) &xmin, sizeof(double));
    infile1.read((char *) &xmax, sizeof(double));
    infile1.read((char *) &ymin, sizeof(double));
    infile1.read((char *) &ymax, sizeof(double));
    infile1.read((char *) &zmin, sizeof(double));
    infile1.read((char *) &zmax, sizeof(double));
    if (nxin != nx || nyin != ny || nzin != nz) {
      cerr << "\a Number of indices do not agree! " << endl;
      exit(1);
    }
    cout << " Comparing with data at time " << time << "\n" << endl;
    infile1.read((char *) indata, nx*ny*nz*sizeof(double));
    infile1.close();
    for (i = 0; i < nx*ny*nz; i++) data[i] -= indata[i];
  }
  
  double *X, *Y, *Z;
  X = new double[nx];
  Y = new double[ny];
  Z = new double[nz];
  double dd;
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
  dd = (xmax - xmin)/(nx-1);
  for (i = 0; i < nx; i++) X[i] = xmin + i*dd;
  dd = (ymax - ymin)/(ny-1);
  for (j = 0; j < ny; j++) Y[j] = ymin + j*dd;
  dd = (zmax - zmin)/(nz-1);
  for (k = 0; k < nz; k++) Z[k] = zmin + k*dd;
#else
#ifdef Cell
  dd = (xmax - xmin)/nx;
  for (i = 0; i < nx; i++) X[i] = xmin + (i+0.5)*dd;
  dd = (ymax - ymin)/ny;
  for (j = 0; j < ny; j++) Y[j] = ymin + (j+0.5)*dd;
  dd = (zmax - zmin)/nz;
  for (k = 0; k < nz; k++) Z[k] = zmin + (k+0.5)*dd;
#else
#error Not define Vertex nor Cell
#endif  
#endif

  int ext[3];
  ext[0]=nx;
  ext[1]=ny;
  ext[2]=nz;
  void writefile(int *ext, double *XX, double *YY, double *ZZ,double *datain,
             		                      char *filename, const char * flag);

  char outname[50];
  int lev=0;
  switch(optype)
  {
     case 1:
       sprintf(outname,"%s.dat",argv[2]);
       writefile(ext, X, Y, Z, data, outname, argv[1]);
       break;
     case 2:
       double *fx,*fy,*fz;
       fx = new double[nx*ny*nz];
       fy = new double[nx*ny*nz];
       fz = new double[nx*ny*nz];
       f_fderivs(ext,data,fx,fy,fz,X,Y,Z,SYM1,SYM2,SYM3,Symmetry,lev);
       sprintf(outname,"%s_x.dat",argv[2]);
       writefile(ext, X, Y, Z, fx, outname, argv[1]);
       sprintf(outname,"%s_y.dat",argv[2]);
       writefile(ext, X, Y, Z, fy, outname, argv[1]);
       sprintf(outname,"%s_z.dat",argv[2]);
       writefile(ext, X, Y, Z, fz, outname, argv[1]);
       delete[] fx; delete[] fy; delete[] fz;
       break;
     case 3:
       double *fxx,*fxy,*fxz,*fyy,*fyz,*fzz;
       fxx = new double[nx*ny*nz];
       fxy = new double[nx*ny*nz];
       fxz = new double[nx*ny*nz];
       fyy = new double[nx*ny*nz];
       fyz = new double[nx*ny*nz];
       fzz = new double[nx*ny*nz];
       f_fdderivs(ext,data,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM1,SYM2,SYM3,Symmetry,lev);
       sprintf(outname,"%s_xx.dat",argv[2]);
       writefile(ext, X, Y, Z, fxx, outname, argv[1]);
       sprintf(outname,"%s_xy.dat",argv[2]);
       writefile(ext, X, Y, Z, fxy, outname, argv[1]);
       sprintf(outname,"%s_xz.dat",argv[2]);
       writefile(ext, X, Y, Z, fxz, outname, argv[1]);
       sprintf(outname,"%s_yy.dat",argv[2]);
       writefile(ext, X, Y, Z, fyy, outname, argv[1]);
       sprintf(outname,"%s_yz.dat",argv[2]);
       writefile(ext, X, Y, Z, fyz, outname, argv[1]);
       sprintf(outname,"%s_zz.dat",argv[2]);
       writefile(ext, X, Y, Z, fzz, outname, argv[1]);
       delete[] fxx; delete[] fxy; delete[] fxz;
       delete[] fyy; delete[] fyz; delete[] fzz;
       break;
  }

  delete[] data;
  delete[] X;
  delete[] Y;
  delete[] Z;

}

/*-----------------------------------*/
/* get rid of any 4 character suffix */
/*-----------------------------------*/
void set_fname(char *fname)
{
  int len = strlen(fname)-4;
  char * n_fname;
  n_fname = new char[len];

  for (int i = 0; i < len; ++i){n_fname[i] = fname[i]; 
  //     cout << n_fname[i] << " " << i << endl;
  }
  n_fname[len] = '\0';
  
  //      cout << "n_fname: " << n_fname << " fname: " << fname << ", " 
  // 	  << len << endl;
  
  strcpy(fname,n_fname);   /* Send back the old pointer */
  delete n_fname;
}
//|----------------------------------------------------------------------------
//  writefile
//|----------------------------------------------------------------------------
  void writefile(int *ext, double *XX, double *YY, double *ZZ,double *datain,
             		                      char *filename, const char * flag)
{
    int nx = ext[0],ny = ext[1],nz = ext[2];
    int i,j,k;
//|--->open out put file    
    ofstream outfile;
    outfile.open(filename);
   
   if(!strcmp(flag, "X"))
   {
    int exty = ny - 1,extz = nz - 1;
    int qrty = exty/4,qrtz = extz/4;
    int medy = exty/2,medz = extz/2;
    outfile <<"# for Y and Z at "<< setw(10) << setprecision(10) << YY[0] <<","
                                 << setw(10) << setprecision(10) << ZZ[0] << " "
                                 << setw(10) << setprecision(10) << YY[qrty] << ","
                                 << setw(10) << setprecision(10) << ZZ[qrtz] << " "
                                 << setw(10) << setprecision(10) << YY[medy] <<","
                                 << setw(10) << setprecision(10) << ZZ[medz] << " "
                                 << setw(10) << setprecision(10) << YY[exty-qrty] << ","
                                 << setw(10) << setprecision(10) << ZZ[extz-qrtz] << " "
                                 << setw(10) << setprecision(10) << YY[exty] <<","
                                 << setw(10) << setprecision(10) << ZZ[extz] <<endl;
    for (i = 0; i < nx; i++) 
    {
    int ind1 = i;
    int ind2 = i + qrty*nx + qrtz*nx*ny; 
    int ind3 = i + medy*nx + medz*nx*ny; 
    int ind4 = i + (exty - qrty)*nx + (extz - qrtz)*nx*ny; 
    int ind5 = i + exty*nx + extz*nx*ny; 
    outfile << setw(10) << setprecision(10) << XX[i] << " "
	    << setw(16) << setprecision(15) << datain[ind1] << " "
	    << setw(16) << setprecision(15) << datain[ind2] << " "
	    << setw(16) << setprecision(15) << datain[ind3] << " "
	    << setw(16) << setprecision(15) << datain[ind4] << " "
	    << setw(16) << setprecision(15) << datain[ind5] << " "
	    << endl;
    }
   }
   else if(!strcmp(flag, "Z0-X"))
   {
    int exty = ny - 1;
    int qrty = exty/4;
    int medy = exty/2;
    outfile <<"# for Y and Z at "<< setw(10) << setprecision(10) << YY[0] <<","
                                 << setw(10) << setprecision(10) << ZZ[0] << " "
                                 << setw(10) << setprecision(10) << YY[qrty] << ","
                                 << setw(10) << setprecision(10) << ZZ[0] << " "
                                 << setw(10) << setprecision(10) << YY[medy] <<","
                                 << setw(10) << setprecision(10) << ZZ[0] << " "
                                 << setw(10) << setprecision(10) << YY[exty-qrty] << ","
                                 << setw(10) << setprecision(10) << ZZ[0] << " "
                                 << setw(10) << setprecision(10) << YY[exty] <<","
                                 << setw(10) << setprecision(10) << ZZ[0] <<endl;
    for (i = 0; i < nx; i++) 
    {
    int ind1 = i;
    int ind2 = i + qrty*nx; 
    int ind3 = i + medy*nx; 
    int ind4 = i + (exty - qrty)*nx; 
    int ind5 = i + exty*nx; 
    outfile << setw(10) << setprecision(10) << XX[i] << " "
	    << setw(16) << setprecision(15) << datain[ind1] << " "
	    << setw(16) << setprecision(15) << datain[ind2] << " "
	    << setw(16) << setprecision(15) << datain[ind3] << " "
	    << setw(16) << setprecision(15) << datain[ind4] << " "
	    << setw(16) << setprecision(15) << datain[ind5] << " "
	    << endl;
    }
   }
   else if(!strcmp(flag, "Y"))
   {
    int extx = nx - 1,extz = nz - 1;
    int qrtx = extx/4,qrtz = extz/4;
    int medx = extx/2,medz = extz/2;
    outfile <<"# for X and Z at "<< setw(10) << setprecision(10) << XX[0] <<","
                                 << setw(10) << setprecision(10) << ZZ[0] << " "
                                 << setw(10) << setprecision(10) << XX[qrtx] << ","
                                 << setw(10) << setprecision(10) << ZZ[qrtz] << " "
                                 << setw(10) << setprecision(10) << XX[medx] <<","
                                 << setw(10) << setprecision(10) << ZZ[medz] << " "
                                 << setw(10) << setprecision(10) << XX[extx-qrtx] << ","
                                 << setw(10) << setprecision(10) << ZZ[extz-qrtz] << " "
                                 << setw(10) << setprecision(10) << XX[extx] <<","
                                 << setw(10) << setprecision(10) << ZZ[extz] <<endl;
    for (j = 0; j < ny; j++) 
    {
    int ind1 =                 j*nx;
    int ind2 = qrtx          + j*nx + qrtz*nx*ny; 
    int ind3 = medx          + j*nx + medz*nx*ny; 
    int ind4 = (extx - qrtx) + j*nx + (extz - qrtz)*nx*ny; 
    int ind5 = extx          + j*nx + extz*nx*ny; 
    outfile << setw(10) << setprecision(10) << YY[j] << " "
	    << setw(16) << setprecision(15) << datain[ind1] << " "
	    << setw(16) << setprecision(15) << datain[ind2] << " "
	    << setw(16) << setprecision(15) << datain[ind3] << " "
	    << setw(16) << setprecision(15) << datain[ind4] << " "
	    << setw(16) << setprecision(15) << datain[ind5] << " "
	    << endl;
    }
   }
   else if(!strcmp(flag, "Z0-Y"))
   {
    int extx = nx - 1;
    int qrtx = extx/4;
    int medx = extx/2;
    outfile <<"# for X and Z at "<< setw(10) << setprecision(10) << XX[0] <<","
                                 << setw(10) << setprecision(10) << ZZ[0] << " "
                                 << setw(10) << setprecision(10) << XX[qrtx] << ","
                                 << setw(10) << setprecision(10) << ZZ[0] << " "
                                 << setw(10) << setprecision(10) << XX[medx] <<","
                                 << setw(10) << setprecision(10) << ZZ[0] << " "
                                 << setw(10) << setprecision(10) << XX[extx-qrtx] << ","
                                 << setw(10) << setprecision(10) << ZZ[0] << " "
                                 << setw(10) << setprecision(10) << XX[extx] <<","
                                 << setw(10) << setprecision(10) << ZZ[0] <<endl;
    for (j = 0; j < ny; j++) 
    {
    int ind1 =                 j*nx;
    int ind2 = qrtx          + j*nx; 
    int ind3 = medx          + j*nx; 
    int ind4 = (extx - qrtx) + j*nx; 
    int ind5 = extx          + j*nx; 
    outfile << setw(10) << setprecision(10) << YY[j] << " "
	    << setw(16) << setprecision(15) << datain[ind1] << " "
	    << setw(16) << setprecision(15) << datain[ind2] << " "
	    << setw(16) << setprecision(15) << datain[ind3] << " "
	    << setw(16) << setprecision(15) << datain[ind4] << " "
	    << setw(16) << setprecision(15) << datain[ind5] << " "
	    << endl;
    }
   }
   else if(!strcmp(flag, "Z"))
   {
    int extx = nx - 1,exty = ny - 1;
    int qrtx = extx/4,qrty = exty/4;
    int medx = extx/2,medy = exty/2;
    outfile <<"# for X and Y at "<< setw(10) << setprecision(10) << XX[0] <<","
                                 << setw(10) << setprecision(10) << YY[0] << " "
                                 << setw(10) << setprecision(10) << XX[qrtx] << ","
                                 << setw(10) << setprecision(10) << YY[qrty] << " "
                                 << setw(10) << setprecision(10) << XX[medx] <<","
                                 << setw(10) << setprecision(10) << YY[medy] << " "
                                 << setw(10) << setprecision(10) << XX[extx-qrtx] << ","
                                 << setw(10) << setprecision(10) << YY[exty-qrty] << " "
                                 << setw(10) << setprecision(10) << XX[extx] <<","
                                 << setw(10) << setprecision(10) << YY[exty] <<endl;
    for (k = 0; k < nz; k++) 
    {
    int ind1 =                                    k*nx*ny;
    int ind2 = qrtx                   + qrty*nx + k*nx*ny; 
    int ind3 = medx                   + medy*nx + k*nx*ny; 
    int ind4 = (extx - qrtx) + (exty - qrty)*nx + k*nx*ny; 
    int ind5 = extx                   + exty*nx + k*nx*ny; 
    outfile << setw(10) << setprecision(10) << ZZ[k] << " "
	    << setw(16) << setprecision(15) << datain[ind1] << " "
	    << setw(16) << setprecision(15) << datain[ind2] << " "
	    << setw(16) << setprecision(15) << datain[ind3] << " "
	    << setw(16) << setprecision(15) << datain[ind4] << " "
	    << setw(16) << setprecision(15) << datain[ind5] << " "
	    << endl;
    }
   }
   else if(!strcmp(flag, "YZ"))
   {
    int ext = nx - 1;
    int qrt = ext/4;
    int med = ext/2;
    outfile <<"# for X at "<< setw(10) << setprecision(10) << XX[0]   << " "
                           << setw(10) << setprecision(10) << XX[qrt] << " "
                           << setw(10) << setprecision(10) << XX[med] << " "
                           << setw(10) << setprecision(10) << XX[ext-qrt] << " "
                           << setw(10) << setprecision(10) << XX[ext] <<endl;
    for (k = 0; k < nz; k++) 
    {
    	    for (j = 0; j < ny; j++) 
	    {
	  	    int ind1 =       j*nx + k*nx*ny; 
	  	    int ind2 = qrt + j*nx + k*nx*ny; 
	  	    int ind3 = med + j*nx + k*nx*ny; 
	  	    int ind4 = ext - qrt + j*nx + k*nx*ny; 
	  	    int ind5 = ext + j*nx + k*nx*ny; 
	  	    outfile << setw(10) << setprecision(10) << YY[j] << " "
	  		    << setw(10) << setprecision(10) << ZZ[k] << " "
	  		    << setw(16) << setprecision(15) << datain[ind1] << " "
	  		    << setw(16) << setprecision(15) << datain[ind2] << " "
	  		    << setw(16) << setprecision(15) << datain[ind3] << " "
	  		    << setw(16) << setprecision(15) << datain[ind4] << " "
	  		    << setw(16) << setprecision(15) << datain[ind5] << " "
	  		    << endl;
    	    }
    	    outfile << "\n" ;  /* blanck line for gnuplot */
    }
   }
   else if(!strcmp(flag, "XZ"))
   {
    int ext = ny - 1;
    int qrt = ext/4;
    int med = ext/2;
    outfile <<"# for Y at "<< setw(10) << setprecision(10) << YY[0]   << " "
                           << setw(10) << setprecision(10) << YY[qrt] << " "
                           << setw(10) << setprecision(10) << YY[med] << " "
                           << setw(10) << setprecision(10) << YY[ext-qrt] << " "
                           << setw(10) << setprecision(10) << YY[ext] <<endl;
    for (k = 0; k < nz; k++) 
    {
    	    for (i = 0; i < nx; i++) 
	    {
	  	    int ind1 = i          + k*nx*ny; 
	  	    int ind2 = i + qrt*nx + k*nx*ny; 
	  	    int ind3 = i + med*nx + k*nx*ny; 
	  	    int ind4 = i + (ext - qrt)*nx + k*nx*ny; 
	  	    int ind5 = i + ext*nx + k*nx*ny; 
	  	    outfile << setw(10) << setprecision(10) << XX[i] << " "
	  		    << setw(10) << setprecision(10) << ZZ[k] << " "
	  		    << setw(16) << setprecision(15) << datain[ind1] << " "
	  		    << setw(16) << setprecision(15) << datain[ind2] << " "
	  		    << setw(16) << setprecision(15) << datain[ind3] << " "
	  		    << setw(16) << setprecision(15) << datain[ind4] << " "
	  		    << setw(16) << setprecision(15) << datain[ind5] << " "
	  		    << endl;
    	    }
    	    outfile << "\n" ;  /* blanck line for gnuplot */
    }
   }
   else if(!strcmp(flag, "XY"))
   {
    int ext = nz - 1;
    int qrt = ext/4;
    int med = ext/2;
    outfile <<"# for Z at "<< setw(10) << setprecision(10) << ZZ[0]   << " "
                           << setw(10) << setprecision(10) << ZZ[qrt] << " "
                           << setw(10) << setprecision(10) << ZZ[med] << " "
                           << setw(10) << setprecision(10) << ZZ[ext-qrt] << " "
                           << setw(10) << setprecision(10) << ZZ[ext] <<endl;
    for (j = 0; j < ny; j++) 
    {
    	    for (i = 0; i < nx; i++) 
	    {
	  	    int ind1 = i + j*nx; 
	  	    int ind2 = i + j*nx + qrt*nx*ny; 
	  	    int ind3 = i + j*nx + med*nx*ny; 
	  	    int ind4 = i + j*nx + (ext - qrt)*nx*ny; 
	  	    int ind5 = i + j*nx + ext*nx*ny; 
	  	    outfile << setw(10) << setprecision(10) << XX[i] << " "
	  		    << setw(10) << setprecision(10) << YY[j] << " "
	  		    << setw(16) << setprecision(15) << datain[ind1] << " "
	  		    << setw(16) << setprecision(15) << datain[ind2] << " "
	  		    << setw(16) << setprecision(15) << datain[ind3] << " "
	  		    << setw(16) << setprecision(15) << datain[ind4] << " "
	  		    << setw(16) << setprecision(15) << datain[ind5] << " "
	  		    << endl;
    	    }
    	    outfile << "\n" ;  /* blanck line for gnuplot */
    }
   }
   else if(!strcmp(flag, "C"))
   {
    int medx,medy,medz;
    if(fabs(XX[0])<XX[1]-XX[0]) medx=0;   //symmetry case
    else                  medx=nx/2;
    if(fabs(YY[0])<YY[1]-YY[0]) medy=0;
    else                  medy=ny/2;
    if(fabs(ZZ[0])<ZZ[1]-ZZ[0]) medz=0;
    else                  medz=nz/2;

    outfile <<"# for line x at (y,z) = "<< setw(10) << setprecision(10) << YY[medy]   << ","
                                       << setw(10) << setprecision(10) << ZZ[medz] 
	    <<"  for line y at (x,z) = "<< setw(10) << setprecision(10) << XX[medx]   << ","
	                               << setw(10) << setprecision(10) << ZZ[medz]   << ","
	    <<"  for line z at (x,y) = "<< setw(10) << setprecision(10) << XX[medx]   << ","
                                       << setw(10) << setprecision(10) << YY[medy]   <<endl;

    int NN=nx;
    if(NN<ny) NN=ny;
    if(NN<nz) NN=nz;
    for (int ijk = 0; ijk < NN; ijk++) 
    {
      if(ijk<nx)
        outfile << setw(10) << setprecision(10) << XX[ijk] << " "
	        << setw(16) << setprecision(15) << datain[ijk+medy*nx+medz*nx*ny] << " ";
      else
        outfile << setw(10) << setprecision(10) << "NaN" << " "
	        << setw(16) << setprecision(15) << "NaN" << " ";


      if(ijk<ny)
        outfile << setw(10) << setprecision(10) << YY[ijk] << " "
	        << setw(16) << setprecision(15) << datain[medx+ijk*nx+medz*nx*ny] << " ";
      else
        outfile << setw(10) << setprecision(10) << "NaN" << " "
	        << setw(16) << setprecision(15) << "NaN" << " ";

      if(ijk<nz)
        outfile << setw(10) << setprecision(10) << ZZ[ijk] << " "
	        << setw(16) << setprecision(15) << datain[medx+medy*nx+ijk*nx*ny] << " ";
      else
        outfile << setw(10) << setprecision(10) << "NaN" << " "
	        << setw(16) << setprecision(15) << "NaN" << " ";

     outfile << endl;
    }
   }
   else
   {
    cout<<"In output_data: not recognized flag-->"<<flag<<endl;
    exit(0);
   }
    outfile.close();
}
