//#include "xdrfile.h"
//#include "xdrfile_xtc.h"
/// #include "define.h"
//#include "xdrfile.c"
//#include "xdrfile_xtc.c"

//#include "stdlib.h"
//#include "stdio.h"
//#include "string.h"
#include "math.h"

#define dump 1
#define run 100000
#define timestep  0.01 //the unit is ps for PE  and tau for KG(0.01 tau)
//#define tau_e 2.32 //the unit is ns for PE and tau for KG(~10000 tau)
 
//check "N_values.dat" from Z1 analysis to make sure that the following 4 parameters are set correctly:
//#define nChain1 50             //the number of long chains. For pure polymer samples, nChain1=0
//#define chainLen1 100
//#define nChain2 500           //the number of short chains. For pure polymer samples, nChain2=nChain
//#define chainLen2 10
//#define nMol 10000   //nMol = (nChain1*chainLen1 + nChain2*chainLen2);
//#define TOT_NFRAME 10000 //= (run / dump): 1~run / dump, without the first frame of the trj file
//#define TOT_NFRAME 100000
//#define BEG_NFRAME 0   //begin from the (BEG_NFRAME+1)th frame actually
#define BEG_time 0
//#define MEAN 580           //# frames to be used to calculate the coordinates of mean path
#define skipped 1
//#define MEMORY 1200000000    //nMol*TOT_NFRAME / skipped
#define MID 80
#define nMID 1 //nChain * MID
//#define OUTER 5
 
//#define TOT_NFRAME 10 //= (run / dump): 1~run / dump, without the first frame of the trj file
//#define BEG_NFRAME 0   //begin from the (BEG_NFRAME+1)th frame actually
//#define MEAN 2           //# frames to be used to calculate the coordinates of mean path



typedef struct{
	double x, y, z;
}VecR;

typedef struct{
	VecR position;
	int id, type;
}Mol;

#define VSet(v,sx,sy,sz)  \
	(v).x = sx, \
	(v).y = sy, \
	(v).z = sz

#define VSub(v1, v2, v3)   \
	(v1).x = (v2).x - (v3).x, \
	(v1).y = (v2).y - (v3).y, \
	(v1).z = (v2).z - (v3).z

#define VCopy(v1, v2)         \
	(v1).x = (v2).x, \
	(v1).y = (v2).y, \
	(v1).z = (v2).z

#define VDot(v1, v2)             \
	((v1).x * (v2).x + (v1).y * (v2).y + (v1).z * (v2).z)

#define VLen(v)  sqrt (VDot (v, v))

#define VLenSq(v)  VDot (v, v)

#define VAdd(v1, v2, v3)            \
	(v1).x = (v2).x + (v3).x, \
	(v1).y = (v2).y + (v3).y, \
	(v1).z = (v2).z + (v3).z

#define VVAdd(v1, v2)  VAdd (v1, v1, v2)

#define VScale(v, s)               \
	(v).x *= s, \
	(v).y *= s, \
	(v).z *= s


#define VSCopy(v2, s1, v1)                                  \
   (v2).x = (s1) * (v1).x,                                  \
   (v2).y = (s1) * (v1).y,                                  \
   (v2).z = (s1) * (v1).z




/*--------------------------------*
void itoa_Linux(int i, char*string)
{
	int power, j;
	j = i;
	for (power = 1; j >= 10; j /= 10)
		power *= 10;
	for (; power>0; power /= 10)
	{
		*string++ = '0' + i / power;
		i %= power;
	}
	*string = '\0';
}

/*--------------------------------*/










