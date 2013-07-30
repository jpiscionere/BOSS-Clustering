#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "sglib.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#define MAXLEN (5000)
#define PI (3.141592)
#include <string.h>
#define SQR(x) ((x)*(x))
#define SQRT(x) (pow(x,0.5))
#define OMEGA_M (0.274)
#define OMEGA_L (0.726)
#define W_INDEX (-1.0)
#define H_0 (100)
#define SPEED_OF_LIGHT (300000)
#define LITTLE_H (0.7)



int main(int argc, char *argv[])
{

	FILE 	*fp1,*fp2,*fp3,*fp4,*fp5;
	char	*DsDi,*DsRi,*RsDi,*RsRi,*Covariance_File;
	int	i,j,k,N_Bins,N_Jackknife;
	double	djunk,Number_Density,Imaging_Correction;
	double	**DsDi_jackknife_array,
		**DsRi_jackknife_array,
		**RsDi_jackknife_array,
		**RsRi_jackknife_array,
		**covariance_matrix,
		*Delta,
		*wp_final,
		*wp_mean,
		*rp,
		*DsDi_mean,
		*DsRi_mean,
		*RsDi_mean,
		*RsRi_mean,
		*error;
	void    inv (double**,int),filter(double**,int,int);
	DsDi=argv[1];
	DsRi=argv[2];
	RsDi=argv[3];
	RsRi=argv[4];
	sscanf(argv[5],"%d",&N_Bins);
	sscanf(argv[6],"%d",&N_Jackknife);
	sscanf(argv[7],"%lf",&Number_Density);	
	sscanf(argv[8],"%lf",&Imaging_Correction);
	Covariance_File=argv[9];
	N_Bins=N_Bins-1;

	rp=(double *)calloc(N_Bins, sizeof(double));
	error=(double *)calloc(N_Bins, sizeof(double));
	wp_final=(double *)calloc(N_Bins, sizeof(double));
	Delta=(double *)calloc(N_Bins, sizeof(double));	
	DsDi_mean=(double *)calloc(N_Bins, sizeof(double));
	DsRi_mean=(double *)calloc(N_Bins, sizeof(double));
	RsDi_mean=(double *)calloc(N_Bins, sizeof(double));
	RsRi_mean=(double *)calloc(N_Bins, sizeof(double));
	wp_mean=(double *)calloc(N_Bins, sizeof(double));


	covariance_matrix=malloc(N_Bins * sizeof(double *));
	for(i=0;i<N_Bins;i++)
		covariance_matrix[i]=(double *)calloc(N_Bins, sizeof(double));

	DsDi_jackknife_array=malloc(N_Bins * sizeof(double *));        



	for(i=0;i<N_Bins;i++)
                DsDi_jackknife_array[i]=(double *)calloc(N_Jackknife + 1,sizeof(double));	

	DsRi_jackknife_array=malloc(N_Bins * sizeof(double *));

        for(i=0;i<N_Bins;i++)
                DsRi_jackknife_array[i]=(double *)calloc(N_Jackknife + 1,sizeof(double));
	
	RsDi_jackknife_array=malloc(N_Bins * sizeof(double *));

        for(i=0;i<N_Bins;i++)
                RsDi_jackknife_array[i]=(double *)calloc(N_Jackknife + 1,sizeof(double));
	
	RsRi_jackknife_array=malloc(N_Bins * sizeof(double *));

        for(i=0;i<N_Bins;i++)
                RsRi_jackknife_array[i]=(double *)calloc(N_Jackknife + 1,sizeof(double));


	fp1=fopen(DsDi,"r") ;
        assert( fp1 != NULL );
        i=0;
        fprintf(stderr,"Here4.5\n");
        
	
	for(i=0;i<N_Bins;i++){
                fscanf(fp1,"%lf %lf %lf %lf ",&rp[i],&DsDi_mean[i],&djunk,&djunk);
                for(j=0;j<N_Jackknife;j++){
                        fscanf(fp1,"%lf ",&DsDi_jackknife_array[i][j]);
                }
        }

        fclose(fp1);
	

     	fp2=fopen(DsRi,"r") ;
        assert( fp2 != NULL );



        for(i=0;i<N_Bins;i++){
                fscanf(fp2,"%lf %lf %lf %lf ",&rp[i],&DsRi_mean[i],&djunk,&djunk);
                for(j=0;j<N_Jackknife;j++){
                        fscanf(fp2,"%lf ",&DsRi_jackknife_array[i][j]);
                }
        }

        fclose(fp2);

     	fp3=fopen(RsDi,"r") ;
        assert( fp3 != NULL );

        for(i=0;i<N_Bins;i++){
                fscanf(fp3,"%lf %lf %lf %lf ",&rp[i],&RsDi_mean[i],&djunk,&djunk);
                for(j=0;j<N_Jackknife;j++){
                        fscanf(fp3,"%lf ",&RsDi_jackknife_array[i][j]);
                }
        }

        fclose(fp3);

     	fp4=fopen(RsRi,"r") ;
        assert( fp4 != NULL );

        for(i=0;i<N_Bins;i++){
                fscanf(fp4,"%lf %lf %lf %lf ",&rp[i],&RsRi_mean[i],&djunk,&djunk);
                for(j=0;j<N_Jackknife;j++){
                        fscanf(fp4,"%lf ",&RsRi_jackknife_array[i][j]);
                }
        }

        fclose(fp4);

	for(i=0;i<N_Bins;i++){
		for(j=0;j<N_Jackknife;j++){
			wp_final[i]+=((DsDi_jackknife_array[i][j]/DsRi_jackknife_array[i][j])-(RsDi_jackknife_array[i][j]/RsRi_jackknife_array[i][j]));
			
		}
	}



	for(i=0;i<N_Bins;i++){
		wp_final[i]=wp_final[i]*Imaging_Correction/(N_Jackknife*Number_Density);
		fprintf(stderr,"%e\n",wp_final[i]);
	}
		


	for(i=0;i<N_Bins;i++){
		for(j=0;j<N_Jackknife;j++){

			error[i]+=SQR(((DsDi_jackknife_array[i][j]/DsRi_jackknife_array[i][j])-(RsDi_jackknife_array[i][j]/RsRi_jackknife_array[i][j]))*Imaging_Correction/Number_Density - wp_final[i]);
		}

	}


	for(i=0;i<N_Bins;i++)
		error[i]=SQRT(error[i]*(N_Jackknife-1.)/N_Jackknife);


	for(j=0;j<N_Jackknife;j++)
	{
		for(i=0;i<N_Bins;i++){
			Delta[i]=(((DsDi_jackknife_array[i][j]/DsRi_jackknife_array[i][j])-(RsDi_jackknife_array[i][j]/RsRi_jackknife_array[i][j]))*Imaging_Correction/Number_Density - wp_final[i])/error[i];
		}

		for(i=0;i<N_Bins;i++){
                	for(k=0;k<N_Bins;k++){
                        	covariance_matrix[i][k]+=Delta[i]*Delta[k];

                	}

        	}
	}

	

	double diagonal[N_Bins];
	for(i=0;i<N_Bins;i++){
		for(j=0;j<N_Bins;j++){
			covariance_matrix[i][j]*=((N_Jackknife-1.)/N_Jackknife);	
		}
		diagonal[i]=covariance_matrix[i][i];
	//	fprintf(stderr,"diagonal[%d]=%e\n",i,diagonal[i]);
	}


	

	for(i=0;i<N_Bins;i++){
		for(j=0;j<N_Bins;j++){
	//		fprintf(stderr,"%e %e %e covariance[%d][%d]=%e\n",diagonal[i],diagonal[j], diagonal[i]*diagonal[j],i,j,covariance_matrix[i][j]);		
			covariance_matrix[i][j]/=SQRT(diagonal[i]*diagonal[j]);			
			//fprintf(stderr,"%e %e %e covariance[%d][%d]=%e\n",diagonal[i],diagonal[j], diagonal[i]*diagonal[j],i,j,covariance_matrix[i][j]);		
		}
	}

	inv(covariance_matrix,N_Bins);

	fp5=fopen(Covariance_File,"w");
	assert(fp5!=NULL);

	for(i=0;i<N_Bins;i++){
		for(j=0;j<N_Bins;j++){
			fprintf(fp5,"%e ",covariance_matrix[i][j]);
		}
		fprintf(fp5,"\n");
	}
		
	for(i=0;i<N_Bins;i++){
		wp_mean[i]=(DsDi_mean[i]/DsRi_mean[i] - RsDi_mean[i]/RsRi_mean[i])/Number_Density;
	}


	for(i=0;i<N_Bins;i++)
		fprintf(stdout,"%lf %e %e %e\n",rp[i],wp_final[i],wp_mean[i],error[i]);



        free(rp);
        free(error);
        free(wp_final);
        free(Delta);
        free(DsDi_mean);
        free(DsRi_mean);
        free(RsDi_mean);
        free(RsRi_mean);
        free(wp_mean);


        
        for(i=0;i<N_Bins;i++)
                free(covariance_matrix[i]);
        free(covariance_matrix);

        



        for(i=0;i<N_Bins;i++)
                free(DsDi_jackknife_array[i]);
        
	free(DsDi_jackknife_array);

        

        for(i=0;i<N_Bins;i++)
                free(DsRi_jackknife_array[i]);
        free(DsRi_jackknife_array);

        

        for(i=0;i<N_Bins;i++)
               free(RsDi_jackknife_array[i]);
        free(RsDi_jackknife_array);

        

        for(i=0;i<N_Bins;i++)
                free(RsRi_jackknife_array[i]);
        free(RsRi_jackknife_array);

return 0;
}

void inv(double **a_matrix, int dimension){


        int i,j;



        gsl_matrix * C = gsl_matrix_calloc ( dimension, dimension );
        gsl_matrix * C_i = gsl_matrix_calloc ( dimension, dimension );

        for ( i=0;i<dimension;i++ )
                for ( j=0;j<dimension;j++ ){

                        gsl_matrix_set ( C, i, j, a_matrix[i][j] );
                }


        int s;

        gsl_permutation *p = gsl_permutation_alloc ( dimension );

        gsl_linalg_LU_decomp ( C, p, &s );

        gsl_linalg_LU_invert ( C, p, C_i );

        for ( i=0;i<dimension;i++ )
        {
        for ( j=0;j<dimension;j++ )
                {
                        a_matrix[i][j] = gsl_matrix_get( C_i, i, j ) ;

                }
        }


        gsl_permutation_free( p );
        gsl_matrix_free( C );
        gsl_matrix_free( C_i );


}

