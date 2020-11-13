#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>


/*------------------------------------------------------------------*/


static void _die(char *msg, int line, char *file) {
	fprintf(stderr, "Fatal error: %s\n", msg);
	fprintf(stderr, "Death occurred at %s, line %d\n", file, line);
	exit(1);
}
#define die(msg) _die(msg,__LINE__,__FILE__)

static void _die_r(char *msg, int result,int line, char *file) {
	fprintf(stderr, "Fatal error: %s\n", msg);
	fprintf(stderr,"result = %d\n",result);
	fprintf(stderr, "Death occurred at %s, line %d\n", file, line);
	exit(1);
}
#define die_r(msg,res) _die_r(msg,res,__LINE__,__FILE__)

#define MAGIC 1995

//enum { FALSE, TRUE };

#define XX 0
#define YY 1
#define ZZ 2

#ifndef INVSQRT_DONE
#define invsqrt(x) (1.0f/sqrt(x))
#endif

#define Cube(x)  (x)*(x)*(x)
#define Sqr(x)  (x)*(x)


#define AllocMem(a, n, t)  a = (t *) malloc ((n) * sizeof (t))

#define AllocMem2(a, n1, n2, t)                             \
   AllocMem (a, n1, t *);                                   \
   AllocMem (a[0], (n1) * (n2), t);                         \
   for (k = 1; k < n1; k ++) a[k] = a[k - 1] + n2;

 /*
#define VWrap(v,t)                                         \
   if (v[t] >= 0.5 * boxxyz[t])      v[t] -= boxxyz[t];         \
   else if (v[t] < -0.5 * boxxyz[t]) v[t] += boxxyz[t]

#define VWrapAll(v)                                         \
   {VWrap (v, 0);                                           \
   VWrap (v, 1);                                            \
   VWrap (v, 2);}
   */



   /*------------------------------------------------------------------*/