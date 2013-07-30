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
#include "sglib.h"
#include "jackknife_it.h"
#include "jackknife_it.c"
#include <gsl/gsl_integration.h>
//#include "utils.c"

#define MAXLEN (5000)
#define PI (3.141592)
#include <string.h>
#define SQR(x) ((x)*(x))
#define SQRT(x) (pow(x,0.5))
#define OMEGA_M (0.274)
#define OMEGA_L (0.726)
#define	W_INDEX (-1.0)
#define H_0 (100)
#define SPEED_OF_LIGHT (300000)
#define LITTLE_H (0.7)



double f(double x, void * params){

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
		*fp2, /*Imaging Galaxy File */
		*fp3; /*Polygon File*/
	char	*Gxy_Spectro,
		*Gxy_Imaging,
		*Status_File,
		*Polygon_File;
	int	N_Jackknife, /*Number of jackknife samples*/
		N_Bins; /*Number of log bins */
	double	Start_Bin, /*Location of the edge of the smallest bin */
		Max_Separation, /*Maximum rp Separation */
		log_Bin_Size, /*Rp Bin Size in log*/
		Minimum_Redshift=1000.0; /*Used to calculated maximum serapartion to filter pairs.*/
	int	Normalization_Choice; /*Which normalization should be used for the imaging catalogue 1= Di 2=Ri */
	/* Spectroscopic Galaxy/Randoms Information */

	int	Spectro_Size=5E6; /*This is the assumed length of the galaxy file */
	double	*RA_s, /* Given */
		*Dec_s, /* Given */
		*Redshift_s, /*Given */
		*Weight_s, /*The Fiber Collision or Completeness Weight of The Galaxy/Randoms */	
		*Distance_s; /* Conversion from Redshift */
	int	*Sector_s, /*Given, use mangle's polyid function to get this if you don't have it. */ 
		*Jackknife_s; /*Calculated in jackknife_it.c */
	double	*X_s,*Y_s,*Z_s; /*The cartesian elements to calculate cos_Theta*/


	/* Imaging Galaxy/Randoms Information */

	int 	Imaging_Size=6E6; /*This is the assumed length of the imaging file */
	double  *RA_i, /* Given */
		*Dec_i; /* Given */

	int	*Sector_i, /*Given, use mangle's polyid function to get this if you don't have it. */
		*Jackknife_i; /*Calculated in jackknife_it.c */ 
	double	*X_i,*Y_i,*Z_i;
	/* Wp calculation information */

	double	**DD,	/*This is not an int because the counts will be weights. It is the shape Nbins X NJackknife */
		*Error,
		*Mean,
		Maximum_Dec_Separation, /*Filter by this dec difference */
		Distance_to_Redshift=1646., /*distance to inner redshift bin */
		cos_Theta,
		rp;
	int	bin;
/*Random Counters and Such */
	int	i=0,j=0,k=0;
	int	Ngal_s=0; /*Number of Galaxies/Randoms in the Spectro Sample */
	int 	Ngal_i=0; /*Number of Galaxies/Randoms in the Imagin Sample */



/*GSL Numerical Integration Crap */
	gsl_integration_workspace * w 
         = gsl_integration_workspace_alloc (1000);
	double result, error,redshift_gsl;
	gsl_function F;
       	F.function = &f;
       	F.params = &redshift_gsl;




/*Read in Args */
	fprintf(stderr,"Here1\n");
	Gxy_Spectro=argv[1];
	Gxy_Imaging=argv[2];
	Polygon_File=argv[3];
	Status_File=argv[4];
	sscanf(argv[5],"%lf",&Start_Bin);
	sscanf(argv[6],"%lf",&Max_Separation);
	sscanf(argv[7],"%d",&N_Bins);
	sscanf(argv[8],"%d",&N_Jackknife);
	sscanf(argv[9],"%d",&Normalization_Choice);
	fprintf(stderr,"Here2\n");

	log_Bin_Size=(log10(Max_Separation)-log10(Start_Bin))/(N_Bins-1.);
	fprintf(stderr,"Log Bin size = %lf \n",log_Bin_Size);
//////////////////////////////*Allocate the Arrays that are going to be used *//////////////////////////////////////////////
	fprintf(stderr,"Here3\n");

/*Spectro Arrays*/
	RA_s=(double *)calloc(Spectro_Size,sizeof(double));
	Dec_s=(double *)calloc(Spectro_Size,sizeof(double));
	Redshift_s=(double *)calloc(Spectro_Size,sizeof(double));
	Distance_s=(double *)calloc(Spectro_Size,sizeof(double));
	X_s=(double *)calloc(Spectro_Size,sizeof(double));
	Y_s=(double *)calloc(Spectro_Size,sizeof(double));
	Z_s=(double *)calloc(Spectro_Size,sizeof(double));

	Weight_s=(double *)calloc(Spectro_Size,sizeof(double));
	Sector_s=(int *)calloc(Spectro_Size,sizeof(int));
	Jackknife_s=(int *)calloc(Spectro_Size,sizeof(int));
/*Imaging Arrays */
        RA_i=(double *)calloc(Imaging_Size,sizeof(double));
        Dec_i=(double *)calloc(Imaging_Size,sizeof(double)); 
       	X_i=(double *)calloc(Imaging_Size,sizeof(double));
        Y_i=(double *)calloc(Imaging_Size,sizeof(double));
        Z_i=(double *)calloc(Imaging_Size,sizeof(double));

        Sector_i=(int *)calloc(Imaging_Size,sizeof(int));
        Jackknife_i=(int *)calloc(Imaging_Size,sizeof(int));

/*Wp Measurement Arrays */

	Error=(double *)calloc(N_Bins,sizeof(double));
	Mean=(double *)calloc(N_Bins,sizeof(double));
	DD=malloc(N_Bins * sizeof(double *));
	for(i=0;i<N_Bins;i++)
		DD[i]=(double *)calloc(N_Jackknife + 1,sizeof(double));
	fprintf(stderr,"Here4\n");
/////////////////////////////* [ READ IN THE GALAXY FILES AND CONVERT REDSHIFTS TO MPC ] *////////////////////////////////////
  	
/*Read in Spectro Sample*/

	fp1=fopen(Gxy_Spectro,"r") ;
        assert( fp1 != NULL );
        i=0;
	fprintf(stderr,"Here4.5\n");
	while(fscanf(fp1,"%lf %lf %lf %lf %d",&RA_s[i],&Dec_s[i],&Redshift_s[i],&Weight_s[i],&Sector_s[i])!=EOF)
		{

			i++;
		}


	fclose(fp1);
	fprintf(stderr,"Here5\n");

	Ngal_s=i;
	if(Ngal_s >= Spectro_Size)
		{
			fprintf(stderr,"BOSS Wp > Something Terrible Has Happened: SPECTROSCOPIC FILE TOO LONG!!!");
			return -1;
		
		}
	fprintf(stderr,"There are %d Galaxies in the Spectro Sample.\n",Ngal_s);	


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


	fp3=fopen(Status_File,"w");
	assert(fp3!=NULL);
	fprintf(stderr,"Here6\n");
	for(i=0;i<Ngal_s;i++)
		{


			
			 gsl_integration_qags (&F, 0, Redshift_s[i], 0, 1e-7, 1000,
                             w, &result, &error);
			Distance_s[i]=result;
			if(Redshift_s[i] < Minimum_Redshift){
				Distance_to_Redshift=Distance_s[i];		
				Minimum_Redshift=Redshift_s[i];
			}

			fprintf(fp3,"Redshift %lf Distance %lf\n",Redshift_s[i],Distance_s[i]);
		}

	
	fprintf(stderr,"The Distance to the closest redshift is %lf\n",Distance_to_Redshift);
	fprintf(stderr,"The Maximum Separation you decided is %lf\n",Max_Separation);	

	Maximum_Dec_Separation=asin(Max_Separation/(2*Distance_to_Redshift))*2.*180./PI*1.00002; //The maximum separation that can happen and let's multiply it by 20% more	

	fprintf(stderr,"Maximum Dec Separation is %lf\n",Maximum_Dec_Separation);


/*Read in Imaging File*/
	fprintf(stderr,"Here7\n");

        fp2=fopen(Gxy_Imaging,"r") ;
        assert( fp2 != NULL );
        i=0;

	while(fscanf(fp2,"%lf %lf %lf %lf %d",&RA_s[i],&Dec_s[i],&Redshift_s[i],&Weight_s[i],&Sector_s[i])!=EOF)
                {

                        i++;
                }


        fclose(fp2);

        Ngal_i=i;
        
	if(Ngal_i >= Imaging_Size)
                {
                        fprintf(stderr,"BOSS Wp > Something Terrible Has Happened: IMAGING FILE TOO LONG!!!\n");
                        return -1;

                }

	fprintf(stderr,"There are %d Galaxies in the Imaging Sample.\n",Ngal_i);


/*
 *This is where the jackknife call is going to go.
 *It's going to take the map file,the number of jackknife samples and the observed sectors in the same order as the observed galaxies.
 *It will return the vector of jackknife ID's in the same order the sector list was given to it.
 *The jackknife ID corresponds to the *one* jackknife sample that galaxy doesn't belong in.

*/



	jackknife_it(N_Jackknife,Polygon_File,Sector_s,Jackknife_s,Ngal_s,RA_s,Dec_s);
	jackknife_it(N_Jackknife,Polygon_File,Sector_i,Jackknife_i,Ngal_i,RA_i,Dec_i);

	fprintf(stderr,"Here8\n");
////////////////////////////////////****Calculation of Wp****/////////////////////////////////////////////////////////////////////////
	for(i=0;i<Ngal_s;i++){
		X_s[i]=sin((90-Dec_s[i]) * PI/180)*cos(RA_s[i] * PI/180) ;
		Y_s[i]=sin((90-Dec_s[i]) * PI/180)*sin(RA_s[i] * PI/180) ;
		Z_s[i]=cos((90-Dec_s[i]) * PI/180) ;
		fprintf(fp3,"Spectro> %lf %lf %lf\n",X_s[i],Y_s[i],Z_s[i]);
	}
	fprintf(stderr,"Here9\n");

	for(i=0;i<Ngal_i;i++){
		X_i[i]=sin((90-Dec_i[i]) * PI/180)*cos(RA_i[i] * PI/180) ;
                Y_i[i]=sin((90-Dec_i[i]) * PI/180)*sin(RA_i[i] * PI/180) ;
                Z_i[i]=cos((90-Dec_i[i]) * PI/180) ;
		fprintf(fp3,"Imaging> %lf %lf %lf\n",X_i[i],Y_i[i],Z_i[i]);
	}
	fprintf(stderr,"Here10\n");


	for(i=0;i<Ngal_s;i++){
		for(j=0;j<Ngal_i;j++){
			if(fabs(Dec_s[i]-Dec_i[j]) <= Maximum_Dec_Separation){
				cos_Theta=X_s[i] * X_i[j] + Y_s[i] * Y_i[j] + Z_s[i] * Z_i[j];
				rp=2.0*Distance_s[i]*SQRT((1.0 - cos_Theta)/2.); /* sin(arccos x) = sqrt(1-x^2) */
				//fprintf(stderr,"distance = %lf,cos_Theta=%lf,rp = %lf\n",Distance_s[i],cos_Theta,rp);
				if(rp < Max_Separation && rp>=Start_Bin){
					bin=(int)floor((log10(rp/Start_Bin))/log_Bin_Size);
					DD[bin][0]+=Weight_s[i]; //Put the Count in the Keeping Track Bin//
					DD[bin][Jackknife_s[i]+1]+=Weight_s[i];
					if(Jackknife_i[j]!=Jackknife_s[i]){
						DD[bin][Jackknife_i[j]+1]+=Weight_s[i];
					}
//					fprintf(fp3,"Bin %d, rp %lf, DD[%d][%d]=%lf\n",bin, rp,bin,Jackknife_s[i], DD[bin][Jackknife_s[i]]);
				}
			}
		}

	}

	fprintf(stderr,"Here11\n");

	double Normalization;




//	if(Normalization_Choice==1)
		Normalization=Ngal_s*Ngal_i;;
//	else
//		Normalization=Sum_of_Spectro_Weights*imaging_normalization;	







	for(i=0;i<N_Bins;i++){
		fprintf(fp3,"%lf %lf ",pow(10,log_Bin_Size*(i)+log10(Start_Bin)),DD[i][0]);
		for(j=0;j<N_Jackknife;j++){
			fprintf(fp3,"%e ",DD[i][j+1]);
		}
		fprintf(fp3,"\n");
	}

	for(i=0;i<N_Bins;i++){
		for(j=0;j<N_Jackknife;j++){

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
          		Error[i]+=SQRT(SQR(DD[i][k+1]-Mean[i])*(N_Jackknife-1)/N_Jackknife);
		}

	}



	for(i=0;i<N_Bins;i++){
//		fprintf(stderr,"%lf %e %e %e ",pow(10,(log_Bin_Size*(i)+log10(Start_Bin))),DD[i][0]/(Normalization),Mean[i],Error[i]);
		fprintf(stdout,"%lf %e %e %e ",pow(10,(log_Bin_Size*(i)+log10(Start_Bin))),DD[i][0]/(Normalization),Mean[i],Error[i]);
		for(j=0;j<N_Jackknife;j++){
//			fprintf(stderr,"%e ",DD[i][j+1]);
			fprintf(stdout,"%e ",DD[i][j+1]);
		}
		fprintf(stdout,"\n");
	}
			

	

/* Free ALL the arrays */
	fclose(fp3);

	
	free(RA_i);
	free(Dec_i);
	free(Sector_i);
	free(Jackknife_i);


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
return 0;
}
