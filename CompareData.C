//$Id: CompareData.C,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
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

int main(int argc, char *argv[])
{
  //
  // USE:  CompareData file1 file2
  // compare 3D data, output the bigest difference data and the 
  // correspoind position
  //
  // or
  //       CompareData file
  // compare 3D data with 0, output the bigest difference data and the 
  // correspoind position
  //

  if (argc != 2 && argc != 3) 
  { 
    cout << "\aUsage: ReadData binaryfile1 binaryfile2 \n " 
	 << "compare 3D data, output the bigest difference data and the "
	 << "correspoind position\n"
         << "or \n " 
	 << "        ReadData binaryfile\n"
	 << "compare 3D data with 0, output the bigest difference data and the "
	 << "correspoind position"
	 << endl;
    exit(1);
  } 
  ifstream infile1;
  infile1.open(argv[1]);
  if(!infile1){
    cerr << "\a Can't open " << argv[1] << " for input." << endl;
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

  cout << "\n Reading file : " << argv[1] << endl;
  cout << "\n  Time        : " << time << endl;
  cout << "  Dimensions  : " << setw(16) << nx   << setw(16) << ny   << setw(16) << nz << endl;
  cout << "   xmin, xmax : " << setw(16) << xmin << setw(16) << xmax << endl;
  cout << "   ymin, ymax : " << setw(16) << ymin << setw(16) << ymax << endl;
  cout << "   zmin, zmax : " << setw(16) << zmin << setw(16) << zmax << endl;
  cout << "\n";
  
  double *X, *Y, *Z;
  X = new double[nx];
  Y = new double[ny];
  Z = new double[nz];
  double dd;
  int i = 0, j = 0, k = 0; 
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

  double *data;
  data = new double[nx*ny*nz];
  infile1.read((char *) data, nx*ny*nz*sizeof(double));
  infile1.close();

 if(argc == 3)
 {
  infile1.open(argv[2]);
  if(!infile1){
     cerr << "\a Can't open " << argv[2] << " for input." << endl;
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
    cout << " Comparing with data " << argv[2] << "\n" << endl;
    cout << "\n  Time        : " << time << endl;
    cout << "  Dimensions  : " << setw(16) << nx   << setw(16) << ny   << setw(16) << nz << endl;
    cout << "   xmin, xmax : " << setw(16) << xmin << setw(16) << xmax << endl;
    cout << "   ymin, ymax : " << setw(16) << ymin << setw(16) << ymax << endl;
    cout << "   zmin, zmax : " << setw(16) << zmin << setw(16) << zmax << endl;
    cout << "\n";
    infile1.read((char *) indata, nx*ny*nz*sizeof(double));
    infile1.close();
    for (i = 0; i < nx*ny*nz; i++) data[i] -= indata[i];
    delete[] indata;
 }
  
  double a=fabs(data[0]);
  int index=0;
  for (i = 0; i < nx*ny*nz; i++)
    if(!finite(a) || (finite(fabs(data[i])) && a<fabs(data[i]))) {a=fabs(data[i]);index=i;}
  
  if(a==0) 
  {
        if(argc==3)cout<<argv[1]<<" and "<<argv[2]<<" are identical to each other!"<<endl;
   else if(argc==2)cout<<argv[1]<<" are identical to 0!"<<endl;
  }
  else
  {
        if(argc == 3)cout<<"The bigest difference between "<<argv[1]<<" and "<<argv[2]<<" is:"<<endl;
   else if(argc == 2)cout<<"The bigest difference between "<<argv[1]<<" and 0 is:"<<endl;
   k=index/(nx*ny);
   j=(index-k*nx*ny)/nx;
   i=index-j*nx-k*nx*ny;
   cout<<a<<" at point ("<<i+1<<","<<j+1<<","<<k+1<<") or to say (x,y,z) = ("<<X[i]<<","<<
	   Y[j]<<","<<Z[k]<<")"<<endl;
  }

  delete[] data;
}
