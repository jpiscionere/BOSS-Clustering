/*
args: galaxy_file other_galaxy_file outfile start_bin end_bin n_bins norm
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
#include <string.h>

#include "cellarray.h"

#include "utils.h"



#define MAXLEN          (5000)
#define MAXBUFSIZE      (10000)

//Change these values to the correct double-precision ones
#define PI              (3.141592)
#define DEG_TO_RAD      (0.01745328888)
#define RAD_TO_DEG      (57.2957914331)


#define SQR(x) ((x)*(x))
#define SQRT(x) (pow(x,0.5))
#define OMEGA_M (0.25)
#define OMEGA_L (0.75)
#define	W_INDEX (-1.0)
#define H_0 (100)
#define SPEED_OF_LIGHT (299792)
#define LITTLE_H (1.0)


#define ADD_DIFF_TIME(tstart,tend)                        ((tend.tv_sec - tstart.tv_sec) + 1e-6*(tend.tv_usec-tstart.tv_usec)) 
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

  /*Input args/files */

  FILE 	*fp1, /*Spectroscopic Galaxy File */
    *fp2; /*Imaging Galaxy File */

  char	*Gxy_Spectro,
    *Gxy_Imaging;


  int	N_Bins; /*Number of log bins */
  double	Start_Bin, /*Location of the edge of the smallest bin */
    Max_Separation, /*Maximum rp Separation */
    log_Bin_Size, /*Rp Bin Size in log*/
  Minimum_Redshift=1000.0, /*Used to calculated maximum serapartion to filter pairs.*/
  Maximum_Redshift=0;
   int	Normalization_Choice; /*Which normalization should be used for the imaging catalogue 1= Di 2=Ri */
  /* Spectroscopic Galaxy/Randoms Information */

  int	Spectro_Size=1E5; /*This is the assumed length of the galaxy file */
  double	*RA_s, /* Given */
    *Dec_s, /* Given */
    *Redshift_s, /*Given */
    *Weight_s, /*The Fiber Collision or Completeness Weight of The Galaxy/Randoms */	
    *Distance_s; /* Conversion from Redshift */

  double	*X_s,*Y_s,*Z_s; /*The cartesian elements to calculate cos_Theta*/
  double	area_tot=4.0*PI;
  fprintf(stderr,"ASSUMMING SPHERE GEOMETRY FOR NORMALIZATION CHOICE!!!!!!!!!!!!!\n");
  /* Imaging Galaxy/Randoms Information */

  int 	Imaging_Size=4E5; /*This is the assumed length of the imaging file */
  double  *RA_i, /* Given */
    *Dec_i; /* Given */



  double	*X_i,*Y_i,*Z_i;
  /* Wp calculation information */

  double	*DD,	/*This is not an int because the counts will be weights. It is the shape Nbins X NJackknife */
    Maximum_Dec_Separation, /*Filter by this dec difference */
    Distance_to_Near_Z=1646., /*distance to inner redshift bin */
    Distance_to_Far_Z=0., /*distance to inner redshift bin */
    cos_Theta,
    rp;




  int	bin;
  /*Random Counters and Such */
  int	i=0,j=0,k=0;
  int	Ngal_s=0; /*Number of Galaxies/Randoms in the Spectro Sample */
  int 	Ngal_i=0; /*Number of Galaxies/Randoms in the Imagin Sample */
  /* void gridlink1D(int np,double rmin,double rmax,double rcell,double *z,int *ngrid,int **gridinit,int **gridlist); */
  void gridlink1D_with_struct(int np,double dmin,double dmax,double rcell,double *x1,double *y1,double *z1,double *dec,int *ngrid,cellarray **lattice);

  struct timeval t0,t1;

  int nitems,nread;
  char buffer[MAXBUFSIZE];


  /*Read in Args */
  Gxy_Spectro=argv[1];
  Gxy_Imaging=argv[2];
  sscanf(argv[3],"%lf",&Start_Bin);
  sscanf(argv[4],"%lf",&Max_Separation);
  sscanf(argv[5],"%d",&N_Bins);
  sscanf(argv[6],"%d",&Normalization_Choice);


  log_Bin_Size=(log10(Max_Separation)-log10(Start_Bin))/(N_Bins);
  //log_Bin_Size=(log10(Max_Separation)-log10(Start_Bin))/(N_Bins-1.);
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
  RA_s       = my_calloc(sizeof(*RA_s),Spectro_Size);
  Dec_s      = my_calloc(sizeof(*Dec_s),Spectro_Size);
  Redshift_s = my_calloc(sizeof(*Redshift_s),Spectro_Size);
  Weight_s   = my_calloc(sizeof(*Weight_s),Spectro_Size);





	

  /////////////////////////////* [ READ IN THE GALAXY FILES AND CONVERT REDSHIFTS TO MPC ] *////////////////////////////////////
  	
  /*Read in Spectro Sample*/
  gettimeofday(&t0,NULL);
  fp1 = my_fopen(Gxy_Spectro,"r") ;
  i=0;
  int flag=0,trash_d;
  nitems=5;
  /* while(fscanf(fp1,"%lf %lf %lf %lf %d",&RA_s[i],&Dec_s[i],&Redshift_s[i],&Weight_s[i],&Sector_s[i])!=EOF) { */
  while(fgets(buffer,MAXBUFSIZE,fp1)!=NULL) {
    nread=sscanf(buffer,"%lf %lf %lf %lf %d",&RA_s[i],&Dec_s[i],&Redshift_s[i],&Weight_s[i],&trash_d);
    if (nread == nitems) {
      if(Redshift_s[i] > 10.0) {
	Redshift_s[i]/=SPEED_OF_LIGHT;
	flag=1;
      }
      
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
  X_s        = my_calloc(sizeof(*X_s),Ngal_s);
  Y_s        = my_calloc(sizeof(*Y_s),Ngal_s);
  Z_s        = my_calloc(sizeof(*Z_s),Ngal_s);



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
  /*GSL Numerical Integration Crap */
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
  }
  gsl_integration_workspace_free(w);
  
  fprintf(stderr,"BOSS Wp > Mean Distance = %lf\n",mean_distance/Ngal_s);	
  fprintf(stderr,"BOSS Wp > The Distance to the closest redshift is %lf\n",Distance_to_Near_Z);
  fprintf(stderr,"BOSS Wp > The Distance to the furthest redshift %lf is %lf\n",Maximum_Redshift,Distance_to_Far_Z);
  double dist_range=(Distance_to_Far_Z - Distance_to_Near_Z);
  double Volume1=4./3.*PI*pow(Distance_to_Far_Z,3);
  double Volume2=4./3.*PI*pow(Distance_to_Near_Z,3);
  double Volume=Volume1-Volume2;
  fprintf(stderr,"BOSS Wp > Spherical Volume =%lf\n",Volume);
  fprintf(stderr,"BOSS Wp > Number Density of Spectro Gal =%17.16f\n",Ngal_s/Volume);
  //	fprintf(stderr,"The Maximum Separation you decided is %lf\n",Max_Separation);	

  Maximum_Dec_Separation=asin(Max_Separation/(2*Distance_to_Near_Z))*2.*RAD_TO_DEG*1.00002; //The maximum separation that can happen and let's multiply it by 20% more
  fprintf(stderr,"BOSS Wp > Maximum Dec Separation is %lf\n",Maximum_Dec_Separation);
  
  
  /*Read in Imaging File*/
  /*Imaging Arrays */
  RA_i     = my_calloc(sizeof(*RA_i),Imaging_Size);
  Dec_i    = my_calloc(sizeof(*Dec_i),Imaging_Size);

  

  nitems=3;
  gettimeofday(&t0,NULL);
  fp2=my_fopen(Gxy_Imaging,"r") ;
  i=0;
  while(fgets(buffer,MAXBUFSIZE,fp2)!=NULL) {
    nread = sscanf(buffer,"%lf %lf %d",&RA_i[i],&Dec_i[i],&trash_d);
    if(nread == nitems) {
      i++;
      if(i==Imaging_Size) {
	fprintf(stderr,"Increasing memory allocation for the imaging sample\n");
	Imaging_Size *= MEMORY_INCREASE_FAC;
	RA_i     = my_realloc(RA_i,sizeof(*RA_i),Imaging_Size,"RA_i");
	Dec_i    = my_realloc(Dec_i,sizeof(*Dec_i),Imaging_Size,"Dec_i");

      }
    } else {
      fprintf(stderr,"WARNING: line %d did not contain %d elements - skipping\n",i,nitems);
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
  
  double number_density_of_imaging=Ngal_i/area_tot;
  double distance_squared=0,Normalization=0.0;
  if(Normalization_Choice==1) {
    for(i=0;i<Ngal_s;i++) {
      Normalization+=Weight_s[i];
    } 
  } else {  
    for(i=0;i<Ngal_s;i++){
      distance_squared+=1./SQR(Distance_s[i]);
      Normalization+=number_density_of_imaging*Weight_s[i]*1./SQR(Distance_s[i]);
    }
//    Normalization=239724.8;
//	Normalization=number_density_of_imaging*1.204988;	 
	fprintf(stderr,"Distance Squared = %lf,Normalization =%lf\n",distance_squared,Normalization); 
  }
 


  

  //gridlink the spectroscopic sample
  /*---Gridlink-variables----------------*/
  int ngrid;/* *gridinit1D,*gridlist1D ; */
  double dmin=-90,dmax=90.0;//min/max dec
  double inv_dmax_diff = 1.0/(dmax-dmin);
  cellarray *lattice;
  
  ngrid=0 ;
  /* gridlink1D(Ngal_i,dmin,dmax,Max_Separation,Dec_i,&ngrid,&gridinit1D,&gridlist1D) ; */
  gridlink1D_with_struct(Ngal_i,dmin,dmax,Maximum_Dec_Separation,X_i,Y_i,Z_i,Dec_i,&ngrid,&lattice);
  fprintf(stderr,"gridlink1D done. ngrid= %d\n",ngrid) ;

  ////////////////////////////////////****Calculation of Wp****/////////////////////////////////////////////////////////////////////////
  double rp_sqr=0.0;
  double max_sep_sqr = Max_Separation*Max_Separation;
  double start_bin_sqr = Start_Bin*Start_Bin;
  double inv_start_bin_sqr = 1.0/start_bin_sqr;
  double inv_log_bin_size = 1.0/log_Bin_Size;
  int icen,icell;
  double *x1,*y1,*z1,*dec;
  int *imaging;
  cellarray *cellstruct __attribute__((aligned(ALIGNMENT)));

  int xx=0;
  for(i=0;i<ngrid;i++)
    xx+= lattice[i].nelements;

  if(xx!=Ngal_i) {
    fprintf(stderr,"ERROR: xx=%d is not equal to Ngal_i=%d\n",xx,Ngal_i);
    exit(EXIT_FAILURE);
  }
    
  /*Wp Measurement Arrays */


  DD    = my_calloc(sizeof(*DD),N_Bins);




  gettimeofday(&t0,NULL);
  for(int ispectro=0;ispectro<Ngal_s;ispectro++) {
    if(fmod(ispectro,10000)==0) fprintf(stderr,"%d of %d\n",ispectro,Ngal_s) ;
    icen = (int)(ngrid*(Dec_s[ispectro]-dmin)*inv_dmax_diff);
    if(icen<0) icen = 0 ;
    if(icen>=ngrid) icen = ngrid-1 ;
    for(int ii=-BIN_REFINE_FACTOR;ii<=BIN_REFINE_FACTOR;ii++) {
      icell = icen + ii ;
      /* for(icell=0;icell<ngrid;icell++) { */ // This makes no difference in the output - so the logic is correct
      if(icell>=0 && icell<ngrid)  {
	/*---Loop-over-particles-in-each-cell-----------------*/
	cellstruct=&(lattice[icell]);
	x1 = cellstruct->x;
	y1 = cellstruct->y;
	z1 = cellstruct->z;
	dec = cellstruct->dec;
	imaging = cellstruct->index;
	for(int p=0;p<cellstruct->nelements;p++) {
	  if(fabs(Dec_s[ispectro]-dec[p]) <= Maximum_Dec_Separation) {
	    cos_Theta=X_s[ispectro] * x1[p] + Y_s[ispectro] * y1[p] + Z_s[ispectro] * z1[p];
	    /* rp_sqr=4.0*Distance_s[ispectro]*Distance_s[ispectro]*(1.0 - cos_Theta)*0.5; /\* sin(arccos x) = sqrt(1-x^2) *\/ */
	    rp_sqr=2.0*Distance_s[ispectro]*Distance_s[ispectro]*(1.0 - cos_Theta); /* sin(arccos x) = sqrt(1-x^2) */
	    if(rp_sqr < max_sep_sqr && rp_sqr >= start_bin_sqr) {
	      bin=(int)floor((0.5*log10(rp_sqr*inv_start_bin_sqr))*inv_log_bin_size);
//	      bin=(int)floor((0.5*log10(rp_sqr*inv_start_bin_sqr))*inv_log_bin_size)-1;
	      /* bin=(int)floor((log10(sqrt(rp_sqr)/Start_Bin))/log_Bin_Size); */
	      DD[bin]+=Weight_s[ispectro]; //Put the Count in the Keeping Track Bin//
	      
	      }
	    }
	  }
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
/* 	for(;p<cellstruct->nelements;p++){ 
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



  for(i=0;i<N_Bins;i++) {
    //		fprintf(stderr,"%lf %e %e %e ",pow(10,(log_Bin_Size*(i)+log10(Start_Bin))),DD[i][0]/(Normalization),Mean[i],Error[i]);
    fprintf(stdout,"%lf %e %lf\n",pow(10,(log_Bin_Size*(i)+log10(Start_Bin))),DD[i]/(Normalization),DD[i]);


  }

	




  /* Free ALL the arrays */	
  free(RA_i);
  free(Dec_i);
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

  free(Weight_s);
  free(DD);	




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
