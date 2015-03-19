#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI (3.141592)
#include <string.h>
#include <assert.h>
#include <gsl/gsl_integration.h>
#define OMEGA_M (0.25)
#define OMEGA_L (0.75)
#define W_INDEX (-1.0)
#define H_0 (100)
#define SPEED_OF_LIGHT (299792)
#define LITTLE_H (1.0)
#include <gsl/gsl_interp.h>
#define SQR(x) ((x)*(x))
#define SQRT(x) (pow(x,0.5))


double f(double x, void * params)
{
  double alpha = *(double *) params;
  double omega_matter_part=OMEGA_M*(pow(1.0+x,3.0));
  double omega_lambda_part=OMEGA_L*(pow(1.0+x,(3.0*(1.0+W_INDEX))));
  double Ez= sqrt(omega_matter_part + omega_lambda_part);
  double f=(SPEED_OF_LIGHT/H_0) *(1.0/LITTLE_H) *(1.0/Ez);

  return f;
}




int main(int argc, char *argv[])
{



int i=0,j=0;
int size=200000;

FILE *fp1;
char *extinction;

double ra,dec,z,loglumdist,Mag_G;
double petroMag_r,modelMag_r,psfMag_r;
double Mag_Correction;
double distance,distance_modulus;
double *z_extinct,*delta_g,*g_r;
double galactic_extinction;
z_extinct=(double *) calloc(size,sizeof(double)) ;
delta_g=(double *) calloc(size,sizeof(double)) ;
g_r=(double *) calloc(size,sizeof(double)) ;

double mr,z_calibration=0.2;

gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);
  double result, error,redshift_gsl;
  gsl_function F;
  F.function = &f;
  F.params = &redshift_gsl;

	extinction=argv[1];
        fp1=fopen(extinction,"r") ;
        assert( fp1 != NULL );
        
	fprintf(stderr,"HERE\n");	
	i=0;
	
	while(fscanf(fp3,"%lf %lf %lf ",&z_[i],&delta_g[i],&g_r[i])!=EOF){
 		i++;
	}

        fprintf(stderr,"HERE\n");	

	int n=i;
	
  	gsl_interp *interpolation = gsl_interp_alloc (gsl_interp_linear,n);
  	gsl_interp_init(interpolation, z_, delta_g, n);
  	gsl_interp_accel * accelerator =  gsl_interp_accel_alloc();

	i=0;

//	fprintf(stdout,"ra dec z Mag_G distance loglumdist Mag_Correction distance_modulus petroMag_r mr petroMag_G\n");

	while(fscanf(stdin,"%lf %lf %lf %lf %lf ",&ra, &dec, &z, &flux, &galactic_extinction)!=EOF)
        {



		 gsl_integration_qags (&F, 0, z, 0, 1e-7, 1000,
                          w, &result, &error);

		distance=result;
//		Mag_Correction=intercept+coefficient1*z+coefficient2*z*z+coefficient3*z*z*z;
		Mag_Correction= gsl_interp_eval(interpolation, z_extinct, extinct, z, accelerator);
		loglumdist=log10(distance*1e6*(1+z));
    		distance_modulus=5.0*loglumdist - 5.0 + Mag_Correction + z_calibration;

		if(mag==0)	
			mr=flux-galactic_extinction;	
		else
			mr=22.5-2.5*log10(flux) - galactic_extinction;

		fprintf(stdout,"%lf %lf %lf %lf\n",ra, dec, z, Mag_G);

        }









			


return 0;
}
