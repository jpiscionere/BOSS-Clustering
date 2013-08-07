#ifndef _JACKKNIFE_IT_H
#define _JACKKNIFE_IT_H

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include "utils.h"
#include "sglib.h"

#define PI (3.141592)

#define SQR(x) ((x)*(x))
#define SQRT(x) (pow(x,0.5))

#define ADD_DIFF_TIME(tstart,tend)          ((tend.tv_sec-tstart.tv_sec) + 1e-6*(tend.tv_usec-tstart.tv_usec))

/* void jackknife_it(int N_jackknife, char *Polygon_File, int *Galaxy_Sector_List, int *Galaxy_Jackknife_ID, int Ngal, double *ra, double *dec, double *area_tot); */
void jackknife_it(int N_Jackknife, char *Polygon_File, int *Galaxy_Sector_Ids, int *Galaxy_Jackknife_Ids, int Ngal, double *ra, double *dec, double *area_tot,
		    double *x,double *y,double *z);


#endif
