/*

This is a little utility that takes standard mangle polygon files and turns them into something that can be run thru the pixmap utility.

*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI (3.141592)
#include <string.h>
#define SQR(x) ((x)*(x))
#define SQRT(x) (pow(x,0.5))


int main(int argc, char *argv[])
{


int Nmax=1E5;
int *overlap_polygon,n_overlap_polygons;
int ijunk;
long n_check,n_masks;
long i, j, k,l;
char *polygons_tag,*polygon_tag,*single_tag,*caps_tag, *weight_tag, *str_tag,*pixel_tag,*InputFile;
FILE *fp1;

  

  polygons_tag = malloc(8*sizeof(char));
  polygon_tag = malloc(7*sizeof(char));
  caps_tag = malloc(5*sizeof(char));
  weight_tag = malloc(7*sizeof(char));
  str_tag = malloc(5*sizeof(char));
  single_tag = malloc(sizeof(char));
  pixel_tag = malloc(6*sizeof(char));






fscanf(stdin,"%ld %s\n",&n_masks,polygons_tag);
fprintf(stdout,"%ld polygons\n",n_masks);

int mask_id,pixel,n_circ;
double weight,area;
double x_poly,y_poly,z_poly,dot_poly;


		 	
    for(i=0;i<n_masks;i++) {
        fscanf(stdin,"%s %d %s %d %s %lf %s %d %s %lf %s\n",
               polygon_tag,&mask_id,single_tag,&n_circ,caps_tag,&weight,weight_tag,&pixel,pixel_tag,&area,str_tag);
//		fprintf(stdout,"polygon                 %d ( %d caps,            %8.7lf weight, 1517 pixel,  %17.16lf str):\n",mask_id,n_circ,weight,area);			
	
		fprintf(stderr,"polygon                 %d ( %d caps,            %8.7lf weight, 1517 pixel,  %17.16lf str):\n",mask_id,n_circ,weight,area);			
		fprintf(stdout,"%d %17.15 %17.15lf\n",mask_id,weight,area);		
		for (k=0;k<n_circ;k++) {
       		 	fscanf(stdin,"%lf %lf %lf %lf\n",&x_poly,&y_poly,&z_poly,&dot_poly);
//			fprintf(stdout,"%21.20lf	%21.20lf	%21.20lf	%21.20lf\n",x_poly,y_poly,z_poly,dot_poly);
			
		}

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
