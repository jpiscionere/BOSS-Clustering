#include "jackknife_it.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "sglib.h"

#define PI (3.141592)
#include <string.h>
#define SQR(x) ((x)*(x))
#define SQRT(x) (pow(x,0.5))


void jackknife_it(int N_Jackknife, char *Polygon_File, int *Galaxy_Sector_Ids, int *Galaxy_Jackknife_Ids, int Ngal)
{

	int Nmax=5E5;
	int n_masks;
	int i=0,j=0,k=0;
	char *polygons_tag,*polygon_tag,*single_tag,*caps_tag, *weight_tag, *str_tag,*pixel_tag;
	char *real,*pixelization,*ptag,*snapped,*balkanized;	
	FILE *fp1;
	int  cap_counter=0;
	int count=0,flag=0;
	int n_jackknife=10;


	polygons_tag = malloc(8*sizeof(char));
  	real=malloc(4*sizeof(char));
  	ptag=malloc(2*sizeof(char));
  	pixelization=malloc(12*sizeof(char));
  	snapped=malloc(7*sizeof(char));
  	balkanized=malloc(10*sizeof(char));


  	polygon_tag = malloc(7*sizeof(char));
  	caps_tag = malloc(5*sizeof(char));
  	weight_tag = malloc(7*sizeof(char));
  	str_tag = malloc(5*sizeof(char));
  	single_tag = malloc(sizeof(char));
  	pixel_tag = malloc(6*sizeof(char));

	fp1=fopen(Polygon_File,"r");
	
	fscanf(fp1,"%d %s\n",&n_masks,polygons_tag);
	fscanf(fp1,"%s %d\n",real,&i);
        	if(real[4]!=0){
                	fprintf(stderr,"Jackknife_it > WARNING: You might be using the wrong mangle file, please check!\n");
                	fscanf(fp1,"%s\n",snapped);
                	fscanf(fp1,"%s\n",balkanized);
        	}else{
                	fscanf(fp1,"%s %s\n",pixelization,ptag);
                	fscanf(fp1,"%s\n",snapped);
                	fscanf(fp1,"%s\n",balkanized);
        	}


	int mask_id[Nmax],pixel,n_circ;
	double weight[Nmax],area[Nmax];
	double x_poly,y_poly,z_poly,dot_poly;

	for(i=0;i<n_masks;i++) {
        	fscanf(fp1,"%s %d %s %d %s %lf %s %d %s %lf %s\n",
               		polygon_tag,&mask_id[i],single_tag,&n_circ,caps_tag,&weight[i],weight_tag,&pixel,pixel_tag,&area[i],str_tag);
			cap_counter+=n_circ;
			
			for (k=0;k<n_circ;k++) {
                        	fscanf(fp1,"%lf %lf %lf %lf\n",&x_poly,&y_poly,&z_poly,&dot_poly);
			}
		}


	double sector_area[n_masks],sector_weight[n_masks];
	int unique_ids[n_masks];
	int jackknife_number[n_masks];
	double area_tot=0;

	flag=count=0;

        for(i=0;i<n_masks;i++)
                unique_ids[i]=-1;


	for(i=0;i<n_masks;i++){
                for(k=0;k<Ngal;k++){
                        if(mask_id[i]==Galaxy_Sector_Ids[k]){
                                for(j=0;j<n_masks;j++){
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

	if((n_masks - count) < 2)
        	fprintf(stderr,"Jackknife_it> WARNING: You are using a mangle file with the polygon id's instead of the sector id's. Make sure your galaxy file is consitent.\n");


	#define MULTIPLE_ARRAY_EXCHANGER(type,a,i,j) { SGLIB_ARRAY_ELEMENTS_EXCHANGER(int,unique_ids,i,j); SGLIB_ARRAY_ELEMENTS_EXCHANGER(double,sector_area,i,j); SGLIB_ARRAY_ELEMENTS_EXCHANGER(double,sector_weight,i,j); }

	SGLIB_ARRAY_QUICK_SORT(int,unique_ids, count, SGLIB_NUMERIC_COMPARATOR , MULTIPLE_ARRAY_EXCHANGER);

	for(i=0;i<count;i++){
                area_tot+=sector_area[i];
        }

        double area_bin=area_tot/n_jackknife;

	area_tot=0;

        for(i=0;i<count;i++){
                area_tot+=sector_area[i];
                jackknife_number[i]=floor(area_tot/area_bin);
	}


	for(i=0;i<Ngal;i++){
                for(j=0;j<count;j++){
                        if(Galaxy_Sector_Ids[i]==unique_ids[j]){
                               	Galaxy_Jackknife_Ids[i]=jackknife_number[j]; 
                                break;
                        }
                }
        }

	
	fclose(fp1);


}
