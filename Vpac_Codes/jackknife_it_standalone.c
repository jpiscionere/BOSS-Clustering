/*

This program takes a mangle file and a list of sectors that is the same length as the galaxy list, unique's them and then returns the jackknife sample for each input sector so that it may be matched to the galaxy.

*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "sglib.h"

#define PI (3.141592)
#include <string.h>
#define SQR(x) ((x)*(x))
#define SQRT(x) (pow(x,0.5))


int main(int argc, char *argv[])
{



int Nmax=5E5;
int observed_polygon[Nmax],gxy_jack_id[Nmax],observed_polygon_ids[Nmax],n_observed_polygons;
//int observed_polygon_ids[Nmax],n_observed_polygons;
int ngal_sectors; 
int a;
long n_check,n_masks;
long i=0, j, k=0,l;
char *polygons_tag,*polygon_tag,*single_tag,*caps_tag, *weight_tag, *str_tag,*pixel_tag,*InputFile;
char *real,*pixelization,*ptag,*snapped,*balkanized,*gxy_sect_file;
FILE *fp1;
int  cap_counter=0;
int count=0,flag=0;
int n_jackknife=10;	



  gxy_sect_file=argv[1];  
  sscanf(argv[2],"%d",&n_jackknife);
  fp1=fopen(gxy_sect_file,"r");
/*
	while(fscanf(fp1,"%d",&observed_polygon_ids[k])!=EOF){
				k++;
	}	
  	
	n_observed_polygons=k;

	fclose(fp1);
*/	

  while(fscanf(fp1,"%d",&observed_polygon[k])!=EOF){

		for(i=0;i<count;i++){
			if(observed_polygon[k]==observed_polygon_ids[i])
				flag=1;
		}
	
		if(flag==0){
			observed_polygon_ids[count]=observed_polygon[k];
			count++;
		}
		
		flag=0;
		k++;
	}

   fclose(fp1);

  ngal_sectors=k;
  n_observed_polygons=count;
  fprintf(stderr,"Number of Observed Polygons = %d\n",n_observed_polygons);

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



fscanf(stdin,"%ld %s\n",&n_masks,polygons_tag);
fscanf(stdin,"%s %ld\n",real,&i);
	if(real[4]!=0){
		fprintf(stderr,"You might be using the wrong mangle file, please check!\n"); 
		fscanf(stdin,"%s\n",snapped);
		fscanf(stdin,"%s\n",balkanized);
	}else{
		fscanf(stdin,"%s %s\n",pixelization,ptag);
		fscanf(stdin,"%s\n",snapped);
		fscanf(stdin,"%s\n",balkanized);
	}

int mask_id[Nmax],pixel,n_circ;
double weight[Nmax],area[Nmax];
double x_poly,y_poly,z_poly,dot_poly;

  

		 	
    for(i=0;i<n_masks;i++) {
        fscanf(stdin,"%s %d %s %d %s %lf %s %d %s %lf %s\n",
               polygon_tag,&mask_id[i],single_tag,&n_circ,caps_tag,&weight[i],weight_tag,&pixel,pixel_tag,&area[i],str_tag);
		cap_counter+=n_circ;
		fprintf(stderr,"polygon                 %d ( %d caps,            %8.7lf weight,   %17.16lf str):\n",mask_id,n_circ,weight,area);			
//		fprintf(stdout,"%d %17.15 %17.15lf\n",mask_id,weight,area);		
		for (k=0;k<n_circ;k++) {
       		 	fscanf(stdin,"%lf %lf %lf %lf\n",&x_poly,&y_poly,&z_poly,&dot_poly);
//			fprintf(stdout,"%21.20lf	%21.20lf	%21.20lf	%21.20lf\n",x_poly,y_poly,z_poly,dot_poly);
				
	


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
		for(k=0;k<n_observed_polygons;k++){
			if(mask_id[i]==observed_polygon_ids[k]){
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
	fprintf(stderr,"You are using a mangle file with the polygon id's instead of the sector id's. Make sure your galaxy file is consitent.\n");	

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
//		fprintf(stdout,"%d %8.7lf %17.16lf %d\n",unique_ids[i],sector_weight[i],sector_area[i],jackknife_number[i]);

	}

	for(i=0;i<n_observed_polygons;i++){
		for(j=0;j<count;j++){
			if(observed_polygon_ids[i]==unique_ids[j]){
				fprintf(stdout,"%d %d\n",observed_polygon_ids[i],jackknife_number[j]);
				break;
			}
		}
	}

	


return(0);
}
