//$Id: microdef.h,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
#ifndef MICRODEF_H
#define MICRODEF_H

#include "microdef.fh"

#define dim 3

//#define Cell or Vertex in "icrodef.fh"

#define buffer_width 6

#if(buffer_width < ghost_width)
#error we always assume buffer_width>ghost_width
#endif

#define PACK 1
#define UNPACK 2

//#define max(a,b) (((a) > (b)) ? (a) : (b))
//#define min(a,b) (((a) < (b)) ? (a) : (b))

#define feq(a,b,d) (fabs(a-b)<d)
#define flt(a,b,d) ((a-b)<d)
#define fgt(a,b,d) ((a-b)>d)

#define TINY 1e-10

//#define GaussInt // for Using Gauss-Legendre quadrature in theta direction

#endif   /* MICRODEF_H */
