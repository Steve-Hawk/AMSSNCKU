//$Id: selectplane.C,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
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

int main(int argc, char *argv[])
{    
  for(int i = 1; i < argc; i++)
  {
    if(!(strncmp("-h",argv[i],2)))
    { 
     cout << "\aUsage: selectplane flag value fname xfname yfname zfname\n " 
 	  << "   where: - flag can be X,Y,Z\n"
	  <<"    value means the plane of flag=value"
	  << endl;
     exit(1);
    }
  } 

  double value=atof(argv[2]);

  ifstream infile1;
  infile1.open(argv[3]);
  if(!infile1){
    cerr << "\a Can't open " << argv[3] << " for input." << endl;
    exit(1);
  }
  ifstream infilex;
  infilex.open(argv[4]);
  if(!infilex){
    cerr << "\a Can't open " << argv[4] << " for input." << endl;
    exit(1);
  }
  ifstream infiley;
  infiley.open(argv[5]);
  if(!infiley){
    cerr << "\a Can't open " << argv[5] << " for input." << endl;
    exit(1);
  }
  ifstream infilez;
  infilez.open(argv[6]);
  if(!infilez){
    cerr << "\a Can't open " << argv[6] << " for input." << endl;
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
  cout << "\n Reading file : " << argv[3] << endl;
  cout << "\n  Time        : " << time << endl;
  cout << "  Dimensions  : " << setw(16) << nx   << setw(16) << ny   << setw(16) << nz << endl;
  cout << "   xmin, xmax : " << setw(16) << xmin << setw(16) << xmax << endl;
  cout << "   ymin, ymax : " << setw(16) << ymin << setw(16) << ymax << endl;
  cout << "   zmin, zmax : " << setw(16) << zmin << setw(16) << zmax << endl;
  cout << "\n";

  const double tor=1e-4;
  double timen;
  int nxn, nyn, nzn;
  double xminn, xmaxn, yminn, ymaxn, zminn, zmaxn;
  infilex.seekg(0, ios::beg);
  infilex.read((char *) &timen, sizeof(double));
  infilex.read((char *) &nxn, sizeof(int));
  infilex.read((char *) &nyn, sizeof(int));
  infilex.read((char *) &nzn, sizeof(int));
  infilex.read((char *) &xminn, sizeof(double));
  infilex.read((char *) &xmaxn, sizeof(double));
  infilex.read((char *) &yminn, sizeof(double));
  infilex.read((char *) &ymaxn, sizeof(double));
  infilex.read((char *) &zminn, sizeof(double));
  infilex.read((char *) &zmaxn, sizeof(double));

  cout << "\n Reading file : " << argv[4] << endl;
  cout << "\n  Time        : " << timen << endl;
  cout << "  Dimensions  : " << setw(16) << nxn   << setw(16) << nyn   << setw(16) << nzn << endl;
  cout << "   xmin, xmax : " << setw(16) << xminn << setw(16) << xmaxn << endl;
  cout << "   ymin, ymax : " << setw(16) << yminn << setw(16) << ymaxn << endl;
  cout << "   zmin, zmax : " << setw(16) << zminn << setw(16) << zmaxn << endl;
  cout << "\n";
  /* sanity check */
  if ( 		  nx != nxn || ny != nyn || nz != nzn ||
       fabs(xmin-xminn)>tor || fabs(ymin-yminn)>tor || fabs(zmin-zminn)>tor ||
       fabs(xmax-xmaxn)>tor || fabs(ymax-ymaxn)>tor || fabs(zmax-zmaxn)>tor) {
    cout << "\n" << endl;
    cout << argv[3]<<" and "<<argv[4]<<" do not agree!"; 
    cout << "\n" << endl;
  }

  infiley.seekg(0, ios::beg);
  infiley.read((char *) &timen, sizeof(double));
  infiley.read((char *) &nxn, sizeof(int));
  infiley.read((char *) &nyn, sizeof(int));
  infiley.read((char *) &nzn, sizeof(int));
  infiley.read((char *) &xminn, sizeof(double));
  infiley.read((char *) &xmaxn, sizeof(double));
  infiley.read((char *) &yminn, sizeof(double));
  infiley.read((char *) &ymaxn, sizeof(double));
  infiley.read((char *) &zminn, sizeof(double));
  infiley.read((char *) &zmaxn, sizeof(double));

  cout << "\n Reading file : " << argv[5] << endl;
  cout << "\n  Time        : " << timen << endl;
  cout << "  Dimensions  : " << setw(16) << nxn   << setw(16) << nyn   << setw(16) << nzn << endl;
  cout << "   xmin, xmax : " << setw(16) << xminn << setw(16) << xmaxn << endl;
  cout << "   ymin, ymax : " << setw(16) << yminn << setw(16) << ymaxn << endl;
  cout << "   zmin, zmax : " << setw(16) << zminn << setw(16) << zmaxn << endl;
  cout << "\n";
  /* sanity check */
  if (		  nx != nxn || ny != nyn || nz != nzn ||
       fabs(xmin-xminn)>tor || fabs(ymin-yminn)>tor || fabs(zmin-zminn)>tor ||
       fabs(xmax-xmaxn)>tor || fabs(ymax-ymaxn)>tor || fabs(zmax-zmaxn)>tor) {
    cout << "\n" << endl;
    cout << argv[3]<<" and "<<argv[5]<<" do not agree!"; 
    cout << "\n" << endl;
  }

  infilez.seekg(0, ios::beg);
  infilez.read((char *) &timen, sizeof(double));
  infilez.read((char *) &nxn, sizeof(int));
  infilez.read((char *) &nyn, sizeof(int));
  infilez.read((char *) &nzn, sizeof(int));
  infilez.read((char *) &xminn, sizeof(double));
  infilez.read((char *) &xmaxn, sizeof(double));
  infilez.read((char *) &yminn, sizeof(double));
  infilez.read((char *) &ymaxn, sizeof(double));
  infilez.read((char *) &zminn, sizeof(double));
  infilez.read((char *) &zmaxn, sizeof(double));

  cout << "\n Reading file : " << argv[6] << endl;
  cout << "\n  Time        : " << timen << endl;
  cout << "  Dimensions  : " << setw(16) << nxn   << setw(16) << nyn   << setw(16) << nzn << endl;
  cout << "   xmin, xmax : " << setw(16) << xminn << setw(16) << xmaxn << endl;
  cout << "   ymin, ymax : " << setw(16) << yminn << setw(16) << ymaxn << endl;
  cout << "   zmin, zmax : " << setw(16) << zminn << setw(16) << zmaxn << endl;
  cout << "\n";
  /* sanity check */
  if (		  nx != nxn || ny != nyn || nz != nzn ||
       fabs(xmin-xminn)>tor || fabs(ymin-yminn)>tor || fabs(zmin-zminn)>tor ||
       fabs(xmax-xmaxn)>tor || fabs(ymax-ymaxn)>tor || fabs(zmax-zmaxn)>tor) {
    cout << "\n" << endl;
    cout << argv[3]<<" and "<<argv[6]<<" do not agree!"; 
    cout << "\n" << endl;
  }

  double *data,*X,*Y,*Z;
  data = new double[nx*ny*nz];
  X = new double[nx*ny*nz];
  Y = new double[nx*ny*nz];
  Z = new double[nx*ny*nz];
  int i; 
  infile1.read((char *) data, nx*ny*nz*sizeof(double));
  infile1.close();
  infilex.read((char *) X, nx*ny*nz*sizeof(double));
  infilex.close();
  infiley.read((char *) Y, nx*ny*nz*sizeof(double));
  infiley.close();
  infilez.read((char *) Z, nx*ny*nz*sizeof(double));
  infilez.close();
   
  /* get rid of any 4 character suffix */
  set_fname(argv[3]);
  char outname[50];
  sprintf(outname,"%s.dat",argv[3]);
  ofstream outfile;
  outfile.open(outname);
  bool first=true;
  double tt;

  if(!strcmp(argv[1], "X"))
  {
    outfile<<"# data for plane X = "<<value<<endl;
    for(int i=0;i<nx*ny*nz;i++)
       if(fabs(X[i]-value)<tor) 
       {
	 outfile<<Y[i]<<" "<<Z[i]<<" "<<data[i]<<endl;
         if(first) {first=false; tt=Y[i];}
	 else      {if(tt>Y[i]) outfile<<endl; tt=Y[i];}
       }
  }
  else if(!strcmp(argv[1], "Y"))
  {
    outfile<<"# data for plane Y = "<<value<<endl;
    for(int i=0;i<nx*ny*nz;i++)
       if(fabs(Y[i]-value)<tor)
       {
	 outfile<<X[i]<<" "<<Z[i]<<" "<<data[i]<<endl;
         if(first) {first=false; tt=X[i];}
	 else      {if(tt>X[i]) outfile<<endl; tt=X[i];}
       } 
  }
  else if(!strcmp(argv[1], "Z"))
  {
    outfile<<"# data for plane Z = "<<value<<endl;
    for(int i=0;i<nx*ny*nz;i++)
       if(fabs(Z[i]-value)<tor)
       {
	 outfile<<X[i]<<" "<<Y[i]<<" "<<data[i]<<endl;
         if(first) {first=false; tt=X[i];}
	 else      {if(tt>X[i]) outfile<<endl; tt=X[i];}
       }
  }
  else
  {
    cout<<"not recognized flag-->"<<argv[1]<<endl;
    exit(0);
  }

  delete[] data;
  delete[] X;
  delete[] Y;
  delete[] Z;
}
