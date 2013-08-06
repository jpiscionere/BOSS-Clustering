/*

This is a little utility that takes standard mangle polygon files and turns them into something that can be run thru the pixmap utility.

*/
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI (3.141592)
#include <string.h>
#define SQR(x) ((x)*(x))
#define SQRT(x) (pow(x,0.5))


int main(int argc, char *argv[])
{


int Nmax=1E6;
int *data_polygons,*data_sectors;
int ijunk;
int n_data_polygons;
long n_check,n_masks;
long i=0, j, k,l;
char *polygons_tag,*polygon_tag,*single_tag,*caps_tag, *weight_tag, *str_tag,*pixel_tag,*InputFile;
FILE *fp1;
char *poly2sec;
int flag=0; 



  polygons_tag = malloc(8*sizeof(char));
  polygon_tag = malloc(7*sizeof(char));
  caps_tag = malloc(5*sizeof(char));
  weight_tag = malloc(7*sizeof(char));
  str_tag = malloc(5*sizeof(char));
  single_tag = malloc(sizeof(char));
  pixel_tag = malloc(6*sizeof(char));




   poly2sec=argv[1];
   fp1=fopen(poly2sec,"r") ;
   assert( fp1 != NULL );


	data_polygons=(int *)calloc(Nmax, sizeof(int));
	data_sectors=(int *)calloc(Nmax, sizeof(int));
	
	while(fscanf(fp1,"%d %d ",&data_sectors[i],&data_polygons[i])!=EOF){
		i++;
	}
fprintf(stderr,"HIIIIIIIIIIIIIIIII\n"); 	

n_data_polygons=i;
fscanf(stdin,"%ld %s\n",&n_masks,polygons_tag);
fprintf(stdout,"%ld polygons\n",n_masks);

int mask_id,pixel,n_circ;
double weight,area;
double x_poly,y_poly,z_poly,dot_poly;

flag=0;

    for(i=0;i<n_masks;i++) {
		fscanf(stdin,"%s %d %s %d %s %lf %s %d %s %lf %s\n",
               		polygon_tag,&mask_id,single_tag,&n_circ,caps_tag,&weight,weight_tag,&pixel,pixel_tag,&area,str_tag);

	
                fprintf(stderr,"polygon                 %d ( %d caps,            %8.7lf weight, 1517 pixel,  %17.16lf str):\n",mask_id,n_circ,weight,area);


		for(j=0;j<n_data_polygons;j++){
			if(mask_id == data_polygons[j]){	
				mask_id=data_sectors[j];
				flag=1;
				fprintf(stderr,"HIIIIIIIIIIIIII\n");
				break;
			}
		}	

		if(flag==1){
			fprintf(stdout,"polygon                 %d ( %d caps,            %8.7lf weight, 1517 pixel,  %17.16lf str):\n",data_sectors[j],n_circ,weight,area);			
			for (k=0;k<n_circ;k++) {
       		 		fscanf(stdin,"%lf %lf %lf %lf\n",&x_poly,&y_poly,&z_poly,&dot_poly);
				fprintf(stdout,"%21.20lf	%21.20lf	%21.20lf	%21.20lf\n",x_poly,y_poly,z_poly,dot_poly);
			}		
		}else{
			for (k=0;k<n_circ;k++) {
                               fscanf(stdin,"%lf %lf %lf %lf\n",&x_poly,&y_poly,&z_poly,&dot_poly);
			}			
			
		}		
		flag=0;
	}


	



/*
double sector_area[n_masks],sector_weight[n_masks];
int unique_ids[n_masks],flag=0,count=0;


	for(i=0;i<n_masks;i++)
		unique_ids[i]=-1;


	for(i=0;i<n_masks;i++){
		for(j=0;j<n_masks;j++){
			if(mask_id[i]==unique_ids[j]){
				sector_area[j]+= area[i];
				flag=1;
			}
		}

		if(flag==0){
			unique_ids[count]=mask_id[i];
			sector_area[count]=area[i];
			sector_weight[count]=weight[i];
			count++;
		}
		flag=0;

	}		
	




	
	for(i=0;i<count;i++)
		fprintf(stdout,"%d %8.7lf %17.16lf\n",unique_ids[i],sector_weight[i],sector_area[i]);

*/
return(0);
}
