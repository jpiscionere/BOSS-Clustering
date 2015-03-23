/*
args: galaxy_file other_galaxy_file polygon_file n_jackknife start_bin n_bins
-galaxy file: ra/dec/cz?/Mpc?/sector
sector can be obtains by running polyid.c from mangle or from sdss catalogue
-other galaxy file : same if autocorrelation, different for cross correlation
-polygon_file: where ID's are sector ids, not polygon ids
-start_bin/n_bins
*/

/*
Notes:
-Find better way to allocate space so there isn't too much or too little allocated



*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <string.h>

#include "cellarray.h"

#include "utils.h"
#include "jackknife_it.h"
#include "progressbar.h"
#include <omp.h>
#define CHUNKSIZE   10



#define MAXLEN          (5000)
#define MAXBUFSIZE      (10000)

//Change these values to the correct double-precision ones
#define PI              (3.141592)
#define DEG_TO_RAD      (0.01745328888)
#define RAD_TO_DEG      (57.2957914331)


#define SQR(x) ((x)*(x))
#define SQRT(x) (pow(x,0.5))
#define OMEGA_M (0.3)
#define OMEGA_L (0.7)
#define	W_INDEX (-1.0)
#define H_0 (100)
#define SPEED_OF_LIGHT (299792)
#define LITTLE_H (1)


/* #define ADD_DIFF_TIME(tstart,tend)                        ((tend.tv_sec - tstart.tv_sec) + 1e-6*(tend.tv_usec-tstart.tv_usec)) */
#define MEMORY_INCREASE_FAC                               (1.2)

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

  //////////////////////////***Definitions***////////////////////////////////////////////////////////////////////
   int nthreads=16,chunk=CHUNKSIZE;
  /*Input args/files */

  FILE 	*fp1, /*Spectroscopic Galaxy File */
    *fp2, /*Imaging Galaxy File */
    *fp3,
    *fp4; /*Evolution Correction File*/
  char	*Gxy_Spectro,
    *Gxy_Imaging,
    *Status_File,
    *Polygon_File,
    Evolution_Correction[MAXLEN];
  int	N_Jackknife, /*Number of jackknife samples*/
    N_Bins; /*Number of log bins */
  double	Start_Bin, /*Location of the edge of the smallest bin */
    Max_Separation, /*Maximum rp Separation */
    log_Bin_Size, /*Rp Bin Size in log*/
    Maximum_Redshift=0.0,
    Minimum_Redshift=1000.0; /*Used to calculated maximum serapartion to filter pairs.*/
  int	Normalization_Choice; /*Which normalization should be used for the imaging catalogue 1= Di 2=Ri */
  /* Spectroscopic Galaxy/Randoms Information */

  int	Spectro_Size=1E5; /*This is the assumed length of the galaxy file */
  double	*RA_s, /* Given */
    *Dec_s, /* Given */
    *Redshift_s,
    *Mag_Correction, /*Given */
    *Weight_s, /*The Fiber Collision or Completeness Weight of The Galaxy/Randoms */	
    *Distance_s,
    *Distance_Modulus,
    logLuminosityDistance_s; /* Conversion from Redshift */
  int	*Sector_s, /*Given, use mangle's polyid function to get this if you don't have it. */
    *Jackknife_s; /*Calculated in jackknife_it.c */
  double	*X_s,*Y_s,*Z_s; /*The cartesian elements to calculate cos_Theta*/
  double	area_tot;

  /* Imaging Galaxy/Randoms Information */

  int 	Imaging_Size=4E5; /*This is the assumed length of the imaging file */
  double  *RA_i, /* Given */
    *Dec_i,
    *mag_r,
    extinction; /* Given */

  int	*Sector_i, /*Given, use mangle's polyid function to get this if you don't have it. */
    *Jackknife_i; /*Calculated in jackknife_it.c */ 
  double	*X_i,*Y_i,*Z_i;
  /* Wp calculation information */

  double	**DD,	/*This is not an int because the counts will be weights. It is the shape Nbins X NJackknife */
    *Error,
    *Mean,
    Maximum_Dec_Separation, /*Filter by this dec difference */
    Distance_to_Near_Z=1646.,
    Distance_to_Far_Z=0.0,
    cos_Theta,
    rp; /*distance to inner redshift bin */



  double *jack_x,
	*jack_y,
	*jack_z;

  /*Random Counters and Such */
  int	i=0,j=0,k=0;
  int	Ngal_s=0; /*Number of Galaxies/Randoms in the Spectro Sample */
  int 	Ngal_i=0; /*Number of Galaxies/Randoms in the Imagin Sample */
  /* void gridlink1D(int np,double rmin,double rmax,double rcell,double *z,int *ngrid,int **gridinit,int **gridlist); */
  void gridlink1D_with_struct(int np,double dmin,double dmax,double rcell,double *x1,double *y1,double *z1,double *dec,double *mag,int *ngrid,cellarray **lattice);

  struct timeval t0,t1;

  int nitems,nread;
  char buffer[MAXBUFSIZE];


  /*Read in Args */
  Gxy_Spectro=argv[1];
  Gxy_Imaging=argv[2];
  Polygon_File=argv[3];
  Status_File=argv[4];
  sscanf(argv[5],"%lf",&Start_Bin);
  sscanf(argv[6],"%lf",&Max_Separation);
  sscanf(argv[7],"%d",&N_Bins);
  sscanf(argv[8],"%d",&N_Jackknife);
  sscanf(argv[9],"%d",&Normalization_Choice);
 
 if(Normalization_Choice==1){
	fprintf(stderr,"BOSS Wp > Using Data Imaging Catalogue with Magnitude Cuts.\n");
 }else if(Normalization_Choice==2){
	fprintf(stderr,"BOSS Wp > Using Random Imaging Catalogue with Ri Normalization and No Mag Cuts.\n");
 }else{
	fprintf(stderr,"BOSS Wp > NORMALIZATION CHOICE OUT OF BOUNDS, BAILING OUT\n");
	return -1;
 }

  //log_Bin_Size=(log10(Max_Separation)-log10(Start_Bin))/(N_Bins-1.);
  log_Bin_Size=(log10(Max_Separation)-log10(Start_Bin))/(N_Bins);
  fprintf(stderr,"BOSS Wp > Log Bin size = %lf \n",log_Bin_Size);
  //////////////////////////////*Allocate the Arrays that are going to be used *//////////////////////////////////////////////


/* #ifdef USE_BINLOOKUP */
/*   int *binlookup=NULL; */
/*   const int NBINLOOKUP=5e4; */
/*   binlookup = my_calloc(sizeof(*binlookup),NBINLOOKUP+2); */
/* #ifdef AVOID_SQRT */
/*   setup_squared_bin_lookup(sdss_data_file,&rmin,&rmax,&nbin,NBINLOOKUP,&rupp,binlookup); */
/*   binfac=NBINLOOKUP/(rmax*rmax); */
/* #else */
/*   setup_bin_lookup(sdss_data_file,&rmin,&rmax,&nbin,NBINLOOKUP,&rupp,binlookup); */
/*   binfac=NBINLOOKUP/rmax; */
/* #endif */
/* #endif */



  /*Spectro Arrays*/
  //Variables in the file
  double *z_,*g_r,*delta_g;


  RA_s       = my_calloc(sizeof(*RA_s),Spectro_Size);
  Dec_s      = my_calloc(sizeof(*Dec_s),Spectro_Size);
  Redshift_s = my_calloc(sizeof(*Redshift_s),Spectro_Size);
  Mag_Correction=my_calloc(sizeof(*Mag_Correction),Spectro_Size);
  Weight_s   = my_calloc(sizeof(*Weight_s),Spectro_Size);
  Sector_s   = my_calloc(sizeof(*Sector_s),Spectro_Size);
  jack_x   = my_calloc(sizeof(*jack_x),N_Jackknife);
  jack_y   = my_calloc(sizeof(*jack_y),N_Jackknife);
  jack_z   = my_calloc(sizeof(*jack_z),N_Jackknife);
  delta_g   = my_calloc(sizeof(*delta_g),Spectro_Size);
  z_   = my_calloc(sizeof(*z_),Spectro_Size);
  g_r   = my_calloc(sizeof(*g_r),Spectro_Size);


 snprintf(Evolution_Correction,MAXLEN,"/home/piscioja/Clustering/Boss/Source/Vpac_Codes/correction.txt");
 fp4=fopen(Evolution_Correction,"r") ;
 assert( fp4 != NULL );
 i=0;
  while(fscanf(fp4,"%lf %lf %lf ",&z_[i],&delta_g[i],&g_r[i])!=EOF){
 		i++;
	}

  int n=i;


  gsl_interp *interpolation = gsl_interp_alloc (gsl_interp_linear,n);
  gsl_interp_init(interpolation, z_, delta_g, n);
  gsl_interp_accel * accelerator =  gsl_interp_accel_alloc();


 double mag_1,mag_2,z_calibration=0.2;	

  /////////////////////////////* [ READ IN THE GALAXY FILES AND CONVERT REDSHIFTS TO MPC ] *////////////////////////////////////
  	
  /*Read in Spectro Sample*/
  gettimeofday(&t0,NULL);
  fp1 = my_fopen(Gxy_Spectro,"r") ;
  i=0;
  int flag=0;
  nitems=5;
  /* while(fscanf(fp1,"%lf %lf %lf %lf %d",&RA_s[i],&Dec_s[i],&Redshift_s[i],&Weight_s[i],&Sector_s[i])!=EOF) { */
  while(fgets(buffer,MAXBUFSIZE,fp1)!=NULL) {
    nread=sscanf(buffer,"%lf %lf %lf %lf %d",&RA_s[i],&Dec_s[i],&Redshift_s[i],&Weight_s[i],&Sector_s[i]);
    if (nread == nitems) {
      if(Redshift_s[i] > 10.0) {
	Redshift_s[i]/=SPEED_OF_LIGHT;
	flag=1;
      }
	mag_1= gsl_interp_eval(interpolation, z_, delta_g, Redshift_s[i], accelerator);
        mag_2= gsl_interp_eval(interpolation, z_, g_r, Redshift_s[i], accelerator);
        Mag_Correction[i]=mag_1 - mag_2;  

      
      if(Redshift_s[i] < 0) {
	fprintf(stderr,"BOSS Wp > Warning! Redshift = %lf, NR = %d. Setting to nearly 0.\n",Redshift_s[i],i);
	Redshift_s[i]=0.00001;
      }
      i++;

      if(i==Spectro_Size) {
	fprintf(stderr,"Increasing memory allocation for the spectroscopic sample\n");
	Spectro_Size *= MEMORY_INCREASE_FAC;
	RA_s       = my_realloc(RA_s,sizeof(*RA_s),Spectro_Size,"RA_s");
	Dec_s      = my_realloc(Dec_s,sizeof(*Dec_s),Spectro_Size,"Dec_s");
	Redshift_s = my_realloc(Redshift_s,sizeof(*Redshift_s),Spectro_Size,"Redshift_s");
	Weight_s   = my_realloc(Weight_s,sizeof(*Weight_s),Spectro_Size,"Weight_s");
	Sector_s   = my_realloc(Sector_s,sizeof(*Sector_s),Spectro_Size,"Sector_s");
        Mag_Correction   = my_realloc(Mag_Correction,sizeof(*Mag_Correction),Spectro_Size,"Mag_Correction");
      }
    } else {
      fprintf(stderr,"WARNING: In spectroscopic sample line %d did not contain %d elements...skipping line\n",i,nitems);
    }
  }
  Ngal_s=i;
  fclose(fp1);
  gettimeofday(&t1,NULL);
  
  if(flag!=0)
    fprintf(stderr,"BOSS Wp > Warning! You gave me cz instead of redshift!\n"); 
	
  //Derived variables
  Distance_s = my_calloc(sizeof(*Distance_s),Ngal_s);
  Distance_Modulus = my_calloc(sizeof(*Distance_Modulus),Ngal_s);
  X_s        = my_calloc(sizeof(*X_s),Ngal_s);
  Y_s        = my_calloc(sizeof(*Y_s),Ngal_s);
  Z_s        = my_calloc(sizeof(*Z_s),Ngal_s);
  Jackknife_s= my_calloc(sizeof(*Jackknife_s),Ngal_s);


  if(Ngal_s >= Spectro_Size) {
    fprintf(stderr,"BOSS Wp > Something Terrible Has Happened: SPECTROSCOPIC FILE TOO LONG!!!");
    return EXIT_FAILURE;
    
  }
  fprintf(stderr,"BOSS Wp > There are %d Galaxies in the Spectro Sample. Time taken = %6.2lf sec\n",Ngal_s,ADD_DIFF_TIME(t0,t1));	


  /*Convert Redshift to Comoving Distance in MPC */
  /* Here I am using Simpsons' Numerical Integration Rule To 
   * convert the redshift of the galaxy into Megaparsecs.
   * The details of the integrals I am using is obviously
   * in Hogg's Distance Measures in Cosmology and you can
   * wikipedia Simpsons' Rule.  I am assuming WMAP7 Cosmology
   * throughout.  You can adjust all those parameters in the header.
   * I'm including an extra parameter (the equation of state of dark energy)
   * because I felt like it.
   */


  double mean_distance=0;	


  gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);
  double result, error,redshift_gsl;
  gsl_function F;
  F.function = &f;
  F.params = &redshift_gsl;

  for(i=0;i<Ngal_s;i++) {
    gsl_integration_qags (&F, 0, Redshift_s[i], 0, 1e-7, 1000,
                          w, &result, &error);
    Distance_s[i]=result;
    if(Redshift_s[i] < Minimum_Redshift) {
      Distance_to_Near_Z=Distance_s[i];
      Minimum_Redshift=Redshift_s[i];
    }
    if(Redshift_s[i] > Maximum_Redshift){
        Distance_to_Far_Z=Distance_s[i];
        Maximum_Redshift=Redshift_s[i];
        }
    mean_distance+=Distance_s[i];
    logLuminosityDistance_s=log10(Distance_s[i]*1e6*(1+Redshift_s[i]));
    Distance_Modulus[i]=5.0*logLuminosityDistance_s - 5 + Mag_Correction[i] + z_calibration;

}



  if(Normalization_Choice==2){
      fprintf(stderr,"BOSS Wp > Not using Mag Cuts for Imaging Randoms, Setting Distance Modulus to 0.\n");
      for(i=0;i<Ngal_s;i++)
              Distance_Modulus[i]=0;

   }

  gsl_integration_workspace_free(w);

  fprintf(stderr,"BOSS Wp > Mean Distance = %lf\n",mean_distance/Ngal_s);
  fprintf(stderr,"BOSS Wp > The Distance to the closest redshift is %lf\n",Distance_to_Near_Z);
  fprintf(stderr,"BOSS Wp > The Distance to the furthest redshift %lf is %lf\n",Maximum_Redshift,Distance_to_Far_Z);
  double dist_range=(Distance_to_Far_Z - Distance_to_Near_Z);
  double Volume1=4./3.*PI*pow(Distance_to_Far_Z,3);
  double Volume2=4./3.*PI*pow(Distance_to_Near_Z,3);
  double Volume=(Volume1-Volume2);
  fprintf(stderr,"BOSS Wp > Spherical Volume =%lf\n",Volume);

  Maximum_Dec_Separation=asin(Max_Separation/(2*Distance_to_Near_Z))*2.*RAD_TO_DEG; //The maximum separation that can happen and let's multiply it by 20% more
  fprintf(stderr,"BOSS Wp > Maximum Dec Separation is %lf\n",Maximum_Dec_Separation);
  
  
  /*Read in Imaging File*/
  /*Imaging Arrays */
  RA_i     = my_calloc(sizeof(*RA_i),Imaging_Size);
  Dec_i    = my_calloc(sizeof(*Dec_i),Imaging_Size);
  Sector_i = my_calloc(sizeof(*Sector_i),Imaging_Size);
  mag_r    = my_calloc(sizeof(*mag_r),Imaging_Size);
  extinction=0;
  

  nitems=5;
  gettimeofday(&t0,NULL);
  fp2=my_fopen(Gxy_Imaging,"r") ;
  i=0;


  if(Normalization_Choice==1){
	nitems=5;
	while(fgets(buffer,MAXBUFSIZE,fp2)!=NULL) {
		nread = sscanf(buffer,"%lf %lf %d %lf %lf",&RA_i[i],&Dec_i[i],&Sector_i[i],&mag_r[i],&extinction);
		mag_r[i]=mag_r[i]-extinction;
    		if(nread == nitems) {
      			i++;
      		if(i==Imaging_Size) {
			fprintf(stderr,"Increasing memory allocation for the imaging sample\n");
			Imaging_Size *= MEMORY_INCREASE_FAC;
			RA_i     = my_realloc(RA_i,sizeof(*RA_i),Imaging_Size,"RA_i");
			Dec_i    = my_realloc(Dec_i,sizeof(*Dec_i),Imaging_Size,"Dec_i");
			Sector_i = my_realloc(Sector_i,sizeof(*Sector_i),Imaging_Size,"Sector_i");
			mag_r    = my_realloc(mag_r,sizeof(*mag_r),Imaging_Size,"mag_r");
      		}
    
	} else {
      			fprintf(stderr,"WARNING: line %d did not contain %d elements - skipping\n",i,nitems);
    		}
  }


} else {

	nitems=3;
	while(fgets(buffer,MAXBUFSIZE,fp2)!=NULL) { 
		nread = sscanf(buffer,"%lf %lf %d ",&RA_i[i],&Dec_i[i],&Sector_i[i]);
		mag_r[i]=-22.2; //Setting mag_r by hand so it passes the Mag_G if condition in the main loop
		if(nread == nitems) { 
			i++;
		if(i==Imaging_Size) { 
			fprintf(stderr,"Increasing memory allocation for the imaging sample\n");
			Imaging_Size *= MEMORY_INCREASE_FAC;
			RA_i     = my_realloc(RA_i,sizeof(*RA_i),Imaging_Size,"RA_i");
			Dec_i    = my_realloc(Dec_i,sizeof(*Dec_i),Imaging_Size,"Dec_i");
			mag_r    = my_realloc(mag_r,sizeof(*mag_r),Imaging_Size,"mag_r");
		}
		
		} else {
				fprintf(stderr,"WARNING: line %d did not contain %d elements - skipping\n",i,nitems);
			}
	}
}

 
  fclose(fp2);
  gettimeofday(&t1,NULL);
  Ngal_i=i;
  if(Ngal_i >= Imaging_Size) {
    fprintf(stderr,"BOSS Wp > Something Terrible Has Happened: IMAGING FILE TOO LONG!!!\n");
    return EXIT_FAILURE;
  }

  X_i   = my_calloc(sizeof(*X_i),Ngal_i);
  Y_i   = my_calloc(sizeof(*Y_i),Ngal_i);
  Z_i   = my_calloc(sizeof(*Z_i),Ngal_i);
  Jackknife_i = my_calloc(sizeof(*Jackknife_i),Ngal_i);

  fprintf(stderr,"BOSS Wp > There are %d Galaxies in the Imaging Sample. Time taken = %6.2lf sec\n",Ngal_i,ADD_DIFF_TIME(t0,t1));


  for(i=0;i<Ngal_s;i++) {
    X_s[i]=sin((90-Dec_s[i]) * DEG_TO_RAD)*cos(RA_s[i] * DEG_TO_RAD) ;
    Y_s[i]=sin((90-Dec_s[i]) * DEG_TO_RAD)*sin(RA_s[i] * DEG_TO_RAD) ;
    Z_s[i]=cos((90-Dec_s[i]) * DEG_TO_RAD) ;
  }


  for(i=0;i<Ngal_i;i++){
    X_i[i]=sin((90-Dec_i[i]) * DEG_TO_RAD)*cos(RA_i[i] * DEG_TO_RAD) ;
    Y_i[i]=sin((90-Dec_i[i]) * DEG_TO_RAD)*sin(RA_i[i] * DEG_TO_RAD) ;
    Z_i[i]=cos((90-Dec_i[i]) * DEG_TO_RAD) ;
  }

    /*
   *This is where the jackknife call is going to go.
   *It's going to take the map file,the number of jackknife samples and the observed sectors in the same order as the observed galaxies.
   *It will return the vector of jackknife ID's in the same order the sector list was given to it.
   *The jackknife ID corresponds to the *one* jackknife sample that galaxy doesn't belong in.

   */

  double r,rmin;
  gettimeofday(&t0,NULL);
  jackknife_it(N_Jackknife,Polygon_File,Sector_s,Ngal_s,jack_x,jack_y,jack_z,&area_tot,X_s,Y_s,Z_s);
  
  for(i=0;i<Ngal_s;i++){
                for(j=0;j<N_Jackknife;j++){
                        r=SQRT(SQR(X_s[i]-jack_x[j]) + SQR(Y_s[i]-jack_y[j])+ SQR(Z_s[i]-jack_z[j]));
                                if(r<rmin){
                                Jackknife_s[i]=j;
                                rmin=r;
                        }
                }
                if(Jackknife_s[i]==-1){
                        fprintf(stderr,"Something Terrible Has Gone Wrong with the Jackknife_s\n");
                }
                rmin=10000000;

        }

   
    for(i=0;i<Ngal_i;i++){
                for(j=0;j<N_Jackknife;j++){
                        r=SQRT(SQR(X_i[i]-jack_x[j]) + SQR(Y_i[i]-jack_y[j])+ SQR(Z_i[i]-jack_z[j]));
                                if(r<rmin){
                                Jackknife_i[i]=j;
                                rmin=r;
                        }
                }
                if(Jackknife_i[i]==-1){
                        fprintf(stderr,"Something Terrible Has Gone Wrong with the Jackknife_i\n");
                }
                rmin=10000000;

        }





  gettimeofday(&t1,NULL);
  fprintf(stderr,"Jack-knifes done - total time = %6.2lf sec\n",ADD_DIFF_TIME(t0,t1));
  
  if(area_tot<0.0001) {
    fprintf(stderr,"BOSS Wp> Something has gone wrong in the jackknife function.  EXITING\n");
    return EXIT_FAILURE;
  }

  Volume=Volume*area_tot / (4*PI);
  fprintf(stderr,"BOSS Wp > Volume of Sample = %lf\n",Volume);
  fprintf(stderr,"BOSS Wp > Number Density of Sample in Redshift Range %lf < z < %lf = %8.7f\n",Minimum_Redshift,Maximum_Redshift,Ngal_s/Volume);


  fp3=my_fopen(Status_File,"w");
  for(i=0;i<Ngal_s;i++)
    fprintf(fp3,"%lf %lf %d\n",RA_s[i],Dec_s[i],Jackknife_s[i]);
 
 
  double number_density_of_imaging=Ngal_i/area_tot;
  double Normalization=0.0;
  if(Normalization_Choice==1) {
    for(i=0;i<Ngal_s;i++) {
      Normalization+=Weight_s[i];
    } 
  } else {  
    for(i=0;i<Ngal_s;i++){

      Normalization+=number_density_of_imaging*Weight_s[i]*1./SQR(Distance_s[i]);
    }
	
	
	fprintf(stderr,"BOSS Wp > Imaging Normalization =%lf\n",Normalization); 
  }
  
  
  

  //gridlink the spectroscopic sample
  /*---Gridlink-variables----------------*/
  int ngrid;/* *gridinit1D,*gridlist1D ; */
  double dmin=-90,dmax=90.0;//min/max dec
  double inv_dmax_diff = 1.0/(dmax-dmin);
  cellarray *lattice;
  
  ngrid=0 ;
  /* gridlink1D(Ngal_i,dmin,dmax,Max_Separation,Dec_i,&ngrid,&gridinit1D,&gridlist1D) ; */
  gridlink1D_with_struct(Ngal_i,dmin,dmax,Maximum_Dec_Separation,X_i,Y_i,Z_i,Dec_i,mag_r,&ngrid,&lattice);
  fprintf(stderr,"gridlink1D done. ngrid= %d\n",ngrid) ;

  ////////////////////////////////////****Calculation of Wp****/////////////////////////////////////////////////////////////////////////

  double max_sep_sqr = Max_Separation*Max_Separation;
  double start_bin_sqr = Start_Bin*Start_Bin;
  double inv_start_bin_sqr = 1.0/start_bin_sqr;
  double inv_log_bin_size = 1.0/log_Bin_Size;

  cellarray *cellstruct __attribute__((aligned(ALIGNMENT)));



  int xx=0;
  for(i=0;i<ngrid;i++)
    xx+= lattice[i].nelements;

  if(xx!=Ngal_i) {
    fprintf(stderr,"ERROR: xx=%d is not equal to Ngal_i=%d\n",xx,Ngal_i);
    exit(EXIT_FAILURE);
  }
    
  /*Wp Measurement Arrays */
  Error = my_calloc(sizeof(*Error),N_Bins);
  Mean  = my_calloc(sizeof(*Mean),N_Bins);
  DD    = my_malloc(sizeof(double *),N_Bins);
  for(i=0;i<N_Bins;i++)
    DD[i] = my_calloc(sizeof(double),N_Jackknife + 1);
  double ***DD_threads =my_malloc(sizeof(double *),N_Bins);
  for(i=0;i<N_Bins;i++)
    DD_threads[i] = my_malloc(sizeof(double),N_Jackknife + 1);
  for(i=0;i<N_Bins;i++)
     for(j=0;j<N_Jackknife + 1;j++)
	    DD_threads[i][j] = my_calloc(sizeof(double),nthreads);






  gettimeofday(&t0,NULL);
  omp_set_num_threads(nthreads);
  int counter=0;
  int interrupted=0;
  init_my_progressbar(Ngal_s,&interrupted);

 #pragma omp parallel default(none) shared(Jackknife_s,Jackknife_i,interrupted,stderr,counter,Ngal_s,Dec_s,Weight_s,X_s,Y_s,Z_s,chunk,ngrid,dmin,inv_dmax_diff,Maximum_Dec_Separation,Distance_s,inv_start_bin_sqr,max_sep_sqr,inv_log_bin_size,start_bin_sqr,DD_threads,Distance_Modulus,lattice) 
  { 
    int tid = omp_get_thread_num();
    #pragma omp for schedule(dynamic,chunk)

  for(int ispectro=0;ispectro<Ngal_s;ispectro++) {
  
	#pragma omp atomic
        counter++;
        if(tid==0){
                        my_progressbar(counter,&interrupted);
                }

 
    int icen = (int)(ngrid*(Dec_s[ispectro]-dmin)*inv_dmax_diff); 
    if(icen<0) icen++;
    if(icen>=ngrid) icen = icen--;
    assert(icen >=0 && icen < ngrid && "icen needs to be in [0, ngrid)");
    for(int ii=-BIN_REFINE_FACTOR;ii<=BIN_REFINE_FACTOR;ii++) {
      int icell = icen + ii ;
      /* for(icell=0;icell<ngrid;icell++) { */ // This makes no difference in the output - so the logic is correct
      if(icell>=0 && icell<ngrid)  {
	/*---Loop-over-particles-in-each-cell-----------------*/
	cellarray *cellstruct=&(lattice[icell]);
	double *x1 = cellstruct->x;
	double *y1 = cellstruct->y;
	double *z1 = cellstruct->z;
	double *dec = cellstruct->dec;
        double *mag = cellstruct ->mag;
	int *imaging = cellstruct->index;
	for(int p=0;p<cellstruct->nelements;p++) {
	  if(fabs(Dec_s[ispectro]-dec[p]) <= Maximum_Dec_Separation) {
	    double cos_Theta=X_s[ispectro] * x1[p] + Y_s[ispectro] * y1[p] + Z_s[ispectro] * z1[p];
	    /* rp_sqr=4.0*Distance_s[ispectro]*Distance_s[ispectro]*(1.0 - cos_Theta)*0.5; /\* sin(arccos x) = sqrt(1-x^2) *\/ */
	    double rp_sqr=2.0*Distance_s[ispectro]*Distance_s[ispectro]*(1.0 - cos_Theta); /* sin(arccos x) = sqrt(1-x^2) */
	    if(rp_sqr < max_sep_sqr && rp_sqr >= start_bin_sqr) {
		double Magnitude=mag[p]-Distance_Modulus[ispectro];
		if(Magnitude >= -23.2 && Magnitude < -21.2){
	      		int bin=(int)floor((0.5*log10(rp_sqr*inv_start_bin_sqr))*inv_log_bin_size);
	      /* bin=(int)floor((log10(sqrt(rp_sqr)/Start_Bin))/log_Bin_Size); */
	      		DD_threads[bin][0][tid]+=Weight_s[ispectro]; //Put the Count in the Keeping Track Bin//
	      		DD_threads[bin][Jackknife_s[ispectro]+1][tid]+=Weight_s[ispectro];
	      		if(Jackknife_i[imaging[p]]!=Jackknife_s[ispectro]){
				DD_threads[bin][Jackknife_i[imaging[p]]+1][tid]+=Weight_s[ispectro];
	      		}
		}
	    }
	  }
	}
      }
    }
  }
}

  finish_myprogressbar(&interrupted);
 

  for(i=0;i<N_Bins;i++){
      for(j=0;j<N_Jackknife + 1;j++){
         for(k=0;k<nthreads;k++){
               DD[i][j]+=DD_threads[i][j][k];
		}
	}
  }
	








  gettimeofday(&t1,NULL);
  fprintf(stderr,"Double loop time in main -> %6.2lf sec \n",ADD_DIFF_TIME(t0,t1));
  
  /* #ifndef USE_AVX */
	/* for(p=0;p<cellstruct->nelements;p++) { */
	/*   if(fabs(Dec_s[ispectro]-dec[p]) <= Maximum_Dec_Separation) { */
	/*     cos_Theta=X_s[ispectro] * x1[p] + Y_s[ispectro] * y1[p] + Z_s[ispectro] * z1[p]; */
	/*     rp_sqr=4.0*Distance_s[ispectro]*Distance_s[ispectro]*(1.0 - cos_Theta)*0.5; /\* sin(arccos x) = sqrt(1-x^2) *\/ */
	/*     if(rp_sqr < max_sep_sqr && rp_sqr >= start_bin_sqr) { */
	/*       bin=(int)floor((0.5*log10(rp_sqr*inv_start_bin_sqr))*inv_log_bin_size); */
	/*       /\* bin=(int)floor((log10(sqrt(rp_sqr)/Start_Bin))/log_Bin_Size); *\/ */
	/*       DD[bin][0]+=Weight_s[ispectro]; //Put the Count in the Keeping Track Bin// */
	/*       DD[bin][Jackknife_s[ispectro]+1]+=Weight_s[ispectro]; */
	/*       if(Jackknife_i[imaging[p]]!=Jackknife_s[ispectro]){ */
	/* 	DD[bin][Jackknife_i[imaging[p]]+1]+=Weight_s[ispectro]; */
	/*       } */
	/*     } */
	/*   } */
	/* } */
/* #else */
/* 	double dec_separation[NVECD]; */
/* 	double rp_sqr_array[NVECD],cos_theta_array[NVECD]; */
/* 	for(p=0;(p+NVECD)<cellstruct->nelements;p+=NVECD) { */
/* 	  #pragma vector always */
/* 	  for(int j=0;j<NVECD;j++) { */
/* 	    dec_separation[j]  = fabs(Dec_s[ispectro]-dec[p]); */
/* 	    cos_theta_array[j] = X_s[ispectro] * x1[p+j] + Y_s[ispectro] * y1[p+j] + Z_s[ispectro] * z1[p+j]; */
/* 	    rp_sqr_array[j]    = 4.0*Distance_s[ispectro]*Distance_s[ispectro]*(1.0 - cos_theta_array[j])*0.5; /\* sin(arccos x) = sqrt(1-x^2) *\/ */
/* 	  } */

/* 	  #pragma novector  */
/* 	  for(int j=0;j<NVECD;j++) { */
/* 	    rp_sqr = rp_sqr_array[j]; */
/* 	    if(dec_separation[j] <= Maximum_Dec_Separation) { */
/* 	      if(rp_sqr < max_sep_sqr && rp_sqr >= start_bin_sqr) { */
/* 		bin=(int)floor((0.5*log10(rp_sqr*inv_start_bin_sqr))*inv_log_bin_size); */
/* 		DD[bin][0]+=Weight_s[ispectro]; //Put the Count in the Keeping Track Bin// */
/* 		DD[bin][Jackknife_s[ispectro]+1]+=Weight_s[ispectro]; */
/* 		if(Jackknife_i[imaging[p+j]]!=Jackknife_s[ispectro]){ */
/* 		  DD[bin][Jackknife_i[imaging[p+j]]+1]+=Weight_s[ispectro]; */
/* 		} */
/* 	      } */
/* 	    } */
/* 	  } */
/* 	} */

/* 	//Now serially process the rest */
/* 	p = p > cellstruct->nelements ? p-NVECD:p; */
/* 	for(;p<cellstruct->nelements;p++){ */
/* 	  if(fabs(Dec_s[ispectro]-dec[p]) <= Maximum_Dec_Separation) { */
/* 	    cos_Theta=X_s[ispectro] * x1[p] + Y_s[ispectro] * y1[p] + Z_s[ispectro] * z1[p]; */
/* 	    rp_sqr=4.0*Distance_s[ispectro]*Distance_s[ispectro]*(1.0 - cos_Theta)*0.5; /\* sin(arccos x) = sqrt(1-x^2) *\/ */
/* 	    if(rp_sqr < max_sep_sqr && rp_sqr >= start_bin_sqr) { */
/* 	      bin=(int)floor((0.5*log10(rp_sqr*inv_start_bin_sqr))*inv_log_bin_size); */
/* 	      DD[bin][0]+=Weight_s[ispectro]; //Put the Count in the Keeping Track Bin// */
/* 	      DD[bin][Jackknife_s[ispectro]+1]+=Weight_s[ispectro]; */
/* 	      if(Jackknife_i[imaging[p]]!=Jackknife_s[ispectro]){ */
/* 		DD[bin][Jackknife_i[imaging[p]]+1]+=Weight_s[ispectro]; */
/* 	      } */
/* 	    } */
/* 	  } */
/* 	} */
/* #endif */

  /* for(int ispectro=0;ispectro<Ngal_s;ispectro++){ */
    /* for(int imaging=0;imaging<Ngal_i;imaging++){ */
    /*   if(fabs(Dec_s[ispectro]-Dec_i[imaging]) <= Maximum_Dec_Separation){ */
    /* 	cos_Theta=X_s[ispectro] * X_i[imaging] + Y_s[ispectro] * Y_i[imaging] + Z_s[ispectro] * Z_i[imaging]; */
    /* 	//rp=2.0*Distance_s[ispectro]*SQRT((1.0 - cos_Theta)/2.); /\* sin(arccos x) = sqrt(1-x^2) *\/ */
    /* 	rp_sqr=4.0*Distance_s[ispectro]*Distance_s[ispectro]*(1.0 - cos_Theta)*0.5; /\* sin(arccos x) = sqrt(1-x^2) *\/ */
    /* 	//fprintf(stderr,"distance = %lf,cos_Theta=%lf,rp = %lf\n",Distance_s[ispectro],cos_Theta,rp); */
    /* 	/\* if(rp < Max_Separation && rp>=Start_Bin){ *\/ */
    /* 	if(rp_sqr < max_sep_sqr && rp_sqr >= start_bin_sqr) { */
    /* 	  /\* bin=(int)floor((log10(rp/Start_Bin))/log_Bin_Size); *\/ */
    /* 	  bin=(int)floor((0.5*log10(rp_sqr*inv_start_bin_sqr))*inv_log_bin_size); */
    /* 	  DD[bin][0]+=Weight_s[ispectro]; //Put the Count in the Keeping Track Bin// */
    /* 	  DD[bin][Jackknife_s[ispectro]+1]+=Weight_s[ispectro]; */
    /* 	  if(Jackknife_i[imaging]!=Jackknife_s[ispectro]){ */
    /* 	    //						fprintf(fp3,"%d %lf %d %d %d %d \n",bin, rp,Jackknife_s[ispectro],Sector_s[ispectro],Jackknife_i[imaging],Sector_i[imaging]); */
    /* 	    DD[bin][Jackknife_i[imaging]+1]+=Weight_s[ispectro]; */
    /* 	  } */
    /* 	} */
    /*   } */
    /* } */
  /* } */


  for(i=0;i<N_Bins;i++){
    fprintf(fp3,"%lf %lf \n",pow(10,log_Bin_Size*(i)+log10(Start_Bin)),DD[i][0]);
  }

  for(i=0;i<N_Bins;i++) {
    for(j=0;j<N_Jackknife;j++) {
      DD[i][j+1] =( DD[i][0] - DD[i][j+1])/(Normalization);
      Mean[i]+=DD[i][j+1];
      fprintf(fp3,"%e ", DD[i][j+1]);	
    }
    fprintf(fp3,"\n");
  }

  for(i=0;i<N_Bins;i++)
    Mean[i]/=N_Jackknife;

  for(i=0;i<N_Bins;i++){
    for(k=0;k<N_Jackknife;k++){
      Error[i]+=SQRT(SQR(DD[i][k+1]-Mean[i])/N_Jackknife);
    }
  }



  for(i=0;i<N_Bins;i++) {
    //		fprintf(stderr,"%lf %e %e %e ",pow(10,(log_Bin_Size*(i)+log10(Start_Bin))),DD[i][0]/(Normalization),Mean[i],Error[i]);
    fprintf(stdout,"%lf %e %e %e ",pow(10,(log_Bin_Size*(i)+log10(Start_Bin))),DD[i][0]/(Normalization),Mean[i],Error[i]);
    for(j=0;j<N_Jackknife;j++) {
      //			fprintf(stderr,"%e ",DD[i][j+1]);
      fprintf(stdout,"%e ",DD[i][j+1]);
    }
    fprintf(stdout,"\n");
  }

	


  fclose(fp3);

  /* Free ALL the arrays */	
  free(RA_i);
  free(Dec_i);
  free(Sector_i);
  free(Jackknife_i);
  free(Jackknife_s);

  free(X_s);
  free(Y_s);
  free(Z_s);
  free(X_i);
  free(Y_i);
  free(Z_i);
  free(RA_s);
  free(Dec_s);
  free(Redshift_s);
  free(Distance_s); 
  free(Sector_s);
  free(Weight_s);

	
  for(i=0;i<N_Bins;i++)
    free(DD[i]);
	
  free(DD);	

  free(Error);	
  free(Mean); 

  for(i=0;i<ngrid;i++) {
    free(lattice[i].x);
    free(lattice[i].y);
    free(lattice[i].z);
    free(lattice[i].dec);
    free(lattice[i].index);
  }
  free(lattice);

  return 0;
}
