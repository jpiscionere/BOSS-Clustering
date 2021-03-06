//Welcome to the hackiest code that ever hacked



#include "jackknife_it.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "sglib.h"

#define PI (3.141592)
#include <string.h>
#define SQR(x) ((x)*(x))
#define SQRT(x) (pow(x,0.5))



void jackknife_it(int N_Jackknife, char *Polygon_File, int *Galaxy_Sector_Ids, int *Galaxy_Jackknife_Ids, int Ngal, double *ra, double *dec, double *area_tot)
{

	fprintf(stderr,"jackknife_it> In Jackknife Function\n");
	
	int Nmax=1E6;
	int n_masks;
	int i=0,j=0,k=0;
	char *polygons_tag,*polygon_tag,*single_tag,*caps_tag, *weight_tag, *str_tag,*pixel_tag;	
	FILE *fp1;
	int  cap_counter=0;
	int count=0,flag=0;
	double *x,*y,*z;




	x=(double *)calloc(Ngal,sizeof(double));
	y=(double *)calloc(Ngal,sizeof(double));
	z=(double *)calloc(Ngal,sizeof(double));



	for(i=0;i<Ngal;i++){

			x[i]=sin((90.-dec[i]) * PI/180.)*cos(ra[i] * PI/180.) ;
                        y[i]=sin((90.-dec[i]) * PI/180.)*sin(ra[i] * PI/180.) ;
                        z[i]=cos((90.-dec[i]) * PI/180.) ;


	}


	polygons_tag = malloc(9*sizeof(char));

  	
  	polygon_tag = malloc(8*sizeof(char));
  	caps_tag = malloc(6*sizeof(char));
  	weight_tag = malloc(8*sizeof(char));
  	str_tag = malloc(6*sizeof(char));
  	single_tag = malloc(2*sizeof(char));
  	pixel_tag = malloc(7*sizeof(char));

	fp1=fopen(Polygon_File,"r");
        assert(fp1!=NULL);	
	fscanf(fp1,"%d %s\n",&n_masks,polygons_tag);
	fprintf(stderr,"jackknife_it> There are %d masks.\n",n_masks);


	int *mask_id,pixel,n_circ;

	double *weight,*area;

	double x_poly,y_poly,z_poly,dot_poly;

	area=(double *)calloc(n_masks,sizeof(double));
	weight=(double *)calloc(n_masks,sizeof(double));
	mask_id=(int *)calloc(n_masks,sizeof(int));


	char check[100];
	fscanf(fp1,"%s",check);	
	while(strncmp("polygon",check,100)!=0)
		fscanf(fp1,"%s",check);


	i=0;
	
	fscanf(fp1,"%d %s %d %s %lf %s %d %s %lf %s",
                        &mask_id[i],single_tag,&n_circ,caps_tag,&weight[i],weight_tag,&pixel,pixel_tag,&area[i],str_tag);	
	
	if(strncmp("str):",str_tag,100)!=0)
	{
		fprintf(stderr,"jackknife_it> str_tag = %s\n",str_tag);	
		fprintf(stderr,"jackknife_it> You are using a different format than I expected.\n");
		fprintf(stderr,"jackknife_it> You are missing a field.  I cannot fix this. You have asked too much of me.\n");

		return ;


	}
	for (k=0;k<n_circ;k++) {
                                fscanf(fp1,"%lf %lf %lf %lf ",&x_poly,&y_poly,&z_poly,&dot_poly);

                                }

		

	for(i=1;i<n_masks;i++) {
        	fscanf(fp1,"%s %d %s %d %s %lf %s %d %s %lf %s",
               		polygon_tag,&mask_id[i],single_tag,&n_circ,caps_tag,&weight[i],weight_tag,&pixel,pixel_tag,&area[i],str_tag);
			
			cap_counter+=n_circ;

			for (k=0;k<n_circ;k++) {
                        	fscanf(fp1,"%lf %lf %lf %lf ",&x_poly,&y_poly,&z_poly,&dot_poly);

				}	
		}


	

	double *sector_area,*sector_weight;
	int *unique_ids,*jackknife_number;
	*area_tot=0;
	int ids;


	sector_area=(double *)calloc(n_masks,sizeof(double));
	sector_weight=(double *)calloc(n_masks,sizeof(double));
	unique_ids=(int *)calloc(n_masks,sizeof(int));
	jackknife_number=(int *)calloc(n_masks,sizeof(int));


	flag=count=ids=0;

        for(i=0;i<n_masks;i++){
                unique_ids[i]=-1;
		jackknife_number[i]=-1;
	}


       
	double *xaverage,*yaverage,*zaverage;
	double *sect_center_ra,*sect_center_dec;
	sect_center_ra=(double *)calloc(n_masks,sizeof(double));
	sect_center_dec=(double *)calloc(n_masks,sizeof(double));
	xaverage=(double *)calloc(n_masks,sizeof(double));
	yaverage=(double *)calloc(n_masks,sizeof(double));
	zaverage=(double *)calloc(n_masks,sizeof(double));



	int n_unique_ids=0;
	flag=0;
	count=0;


		

	for(i=0;i<n_masks;i++){ 
	        for(k=0;k<Ngal;k++){
                        if(mask_id[i]==Galaxy_Sector_Ids[k]){
                                for(j=0;j<count;j++){
                                        if((mask_id[i]==unique_ids[j])){
                                                sector_area[j]+= area[i];
                                                
						flag=1;
                                        }
                                }

                                if((flag==0)){
                                        unique_ids[count]=mask_id[i];
                                        sector_area[count]=area[i];
                                        sector_weight[count]=weight[i];
                                        count++;
                                }

                                flag=0;
                        break;
                        }
                }
        }


	flag=0;


	
	
	fprintf(stderr,"There are %d unique_ids\n",count);
	n_unique_ids=count;

	for(i=0;i<n_unique_ids;i++){
		xaverage[i]=-1.;
		yaverage[i]=-1.;
		zaverage[i]=-1.;

	}


    	for(i=0;i<n_unique_ids;i++){
                for(j=0;j<Ngal;j++){
                        if(Galaxy_Sector_Ids[j]==unique_ids[i]){
                                
				xaverage[i]+=x[j];
                                yaverage[i]+=y[j];
                                zaverage[i]+=z[j];
                                flag++;

                        }
                }

		if(flag==0){
			fprintf(stderr,"Something terrible has happened with xaverage sector id[%d]=%d\n",i,unique_ids[i]);
		}

                xaverage[i]/=flag;
                yaverage[i]/=flag;
                zaverage[i]/=flag;

                flag=0;
         }

	for(i=0;i<n_unique_ids;i++)
	{
		if(xaverage[i]==-1. || yaverage[i]==-1. || zaverage[i]==-1)
			fprintf(stderr,"Something terrible has happened with sector_id[%d]=%d\n",i,unique_ids[i]);

	}


        for(i=0;i<n_unique_ids;i++){
                sect_center_ra[i]=atan2(yaverage[i],xaverage[i]);
		if(sect_center_ra[i] < 0)
			sect_center_ra[i] =2.0*PI + sect_center_ra[i];
		sect_center_ra[i]=180./PI*sect_center_ra[i];
                sect_center_dec[i]=90.- 180./PI * acos(zaverage[i]/SQRT(SQR(xaverage[i]) + SQR(yaverage[i]) + SQR(zaverage[i])));

        }




	if((n_masks - n_unique_ids) < 2)
        	fprintf(stderr,"Jackknife_it> WARNING: You are using a mangle file with the polygon id's instead of the sector id's. Make sure your galaxy file is consitent.\n");
	


	#define MULTIPLE_ARRAY_EXCHANGER(type,a,i,j) { SGLIB_ARRAY_ELEMENTS_EXCHANGER(double,sect_center_ra,i,j); SGLIB_ARRAY_ELEMENTS_EXCHANGER(double,sector_area,i,j); SGLIB_ARRAY_ELEMENTS_EXCHANGER(double,sector_weight,i,j);SGLIB_ARRAY_ELEMENTS_EXCHANGER(int,unique_ids,i,j);SGLIB_ARRAY_ELEMENTS_EXCHANGER(double,sect_center_dec,i,j)}


	SGLIB_ARRAY_QUICK_SORT(double,sect_center_ra, count, SGLIB_NUMERIC_COMPARATOR , MULTIPLE_ARRAY_EXCHANGER);

	for(i=0;i<count;i++){
                *area_tot+=sector_area[i];
        }
	
        double area_bin=*area_tot/N_Jackknife;

	fprintf(stderr,"Total area = %lf, Jackknife area=%lf\n",*area_tot,area_bin);

	*area_tot=0;
	double dec_max=65.,dec_min=-4.;
	int n_dec_bins=floor(SQRT(N_Jackknife));
	double dec_bins_size=(dec_max-dec_min)/n_dec_bins;
	fprintf(stderr,"n_dec_bins=%d\n",n_dec_bins);       
 
      for(j=0;j<n_dec_bins;j++){
		for(i=0;i<n_unique_ids;i++){

			if((sect_center_dec[i] >=(dec_min + j*dec_bins_size)) && (sect_center_dec[i] < (dec_min + (j+1)*dec_bins_size))){

				*area_tot+=sector_area[i];
                		jackknife_number[i]=(int)floor(*area_tot/area_bin);

				if(jackknife_number[i]==N_Jackknife){
					fprintf(stderr,"Jackknife binning slight screw up.  Hacking the Fix\n");
					jackknife_number[i]=N_Jackknife-1;
				}
			}
		}		
	}


	for(i=0;i<Ngal;i++){
                for(j=0;j<n_unique_ids;j++){
                        if(Galaxy_Sector_Ids[i]==unique_ids[j]){
                               	Galaxy_Jackknife_Ids[i]=jackknife_number[j];
								
                                break;
                        }
                }
	
        }



	fclose(fp1);
	free(x);
	free(y);
	free(z);
	free(sect_center_ra);
	free(sect_center_dec);
	free(area);
        free(weight);
        free(mask_id);
    	free(sector_area);
        free(sector_weight);
        free(unique_ids);
        free(jackknife_number);
        free(xaverage);
        free(yaverage);
        free(zaverage);


}
