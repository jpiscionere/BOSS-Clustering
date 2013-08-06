#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI (3.141592)
#include <string.h>
#include <assert.h>
#define SQR(x) ((x)*(x))
#define SQRT(x) (pow(x,0.5))


int main(int argc, char *argv[])
{

int i=0,j,k,Nmax;
FILE *fp1;
char *Gxy_Spectro;
int sector[1000000];
double weight[1000000];
double ra,dec,z,weight_cp;
int sector_in;

	Gxy_Spectro=argv[1];
        fp1=fopen(Gxy_Spectro,"r") ;
        assert( fp1 != NULL );

	while(fscanf(fp1,"%d %lf ",&sector[i],&weight[i])!=EOF)
	{
		i++;
	}

	while(fscanf(stdin,"%lf %lf %lf %lf %d\n",&ra,&dec,&z,&weight_cp,&sector_in)!=EOF)
	{
		for(j=0;j<i;j++)
		{
			if(sector_in==sector[j]){
				fprintf(stdout,"%lf  %lf %lf %lf %d\n",ra,dec,z,weight,sector_in);
				break;	
			}
		
		}	

	
	}
	






return 0;
}
