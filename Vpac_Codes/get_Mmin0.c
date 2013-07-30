/* PROGRAM get_Mmin0
   
   --- get_Mmin0 ngal siglogM logM0 logM1 alpha BGCfile < Halocenters > logMmin
   --- Computes Mmin required to get a desired space density for a particular halo list and HOD.
   --- Version 0 uses the BGC halo format and the Zheng, Coil, Zehavi 2007 HOD model).

      * ngal    = desired space density of galaxies (in h^3/Mpc^3)
      * siglogM = width of cutoff at Mmin (scatter in the mass-luminosity relation)
      * logM0   = minimum halo mass that can contain a satellite galaxy (in units of Msun/h)
      * logM1   = halo mass for which <Nsat>=1 (in units of Msun/h)
      * alpha   = slope of the <Nsat> - M relation
      * BGCfile (BGC format)
      * < Halo center file (output of calc_group_stats.c)
      * > logMmin
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "/home/piscioja/NRC/include/nr.h"
#include "/home/piscioja/halobias/bgc_read_utils.c"
#define MAXLEN (500)
#define sqr(x) ((x)*(x))
#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define mabs(A) ((A) < 0.0 ? -(A) : (A))
#define cnint(x) ((x-floor(x)) < 0.5 ? floor(x) : ceil(x))
#define csign(x) (x < 0.0 ? -1 : 1)
#define PI (3.141592)
#define epsilon 1e-9


int main()
{

double get_Mmin(double , double , double , double , double , char* , char*);
double logMmin;
char bgc_file[MAXLEN],halos_file[MAXLEN];

snprintf(bgc_file,MAXLEN,"/net/bender/data0/LasDamas/Esmeralda/3011/fof_b0p2/Esmeralda_3011_z0p000_fof_b0p2.0000.bgc");
snprintf(halos_file,MAXLEN,"/net/bender/data0/LasDamas/Esmeralda/3011/Esmeralda_3011_z0p000_fof_b0p2.dpp.halos");


logMmin=get_Mmin(0.03,0.1,11.0,12.5,0.1,bgc_file,halos_file);

fprintf(stderr,"get_Mmin0> logMmin = %5.2f\n",logMmin) ;

return 0;


}


double get_Mmin(double rhogal, double siglogM, double logM0, double logM1, double alpha, char *bgc_file, char *halos_file)

{
  int i,j ;
/*---Arguments-------------------------*/

  double M0,M1 ;
  FILE *fp1,*fp2 ;
/*---Halo-files------------------------*/

  OUTPUT_HEADER hdr ;
  int junki,Nhalos,Npart ;
  double junkf,mp,Lbox,Mhalo,logMhalo ;
  int Nbin,*Nhalosbin ;
  double logMbin1,logMbin2,dlogMbin ;
/*---Bias-parameters-------------------*/
  double logMh,Mh,Ncenavg,Nsatavg,Ngalaxies,Ngaltarget ;
  double logMmin,fNgal ;
  void Printhelp(void) ;

/*---Read-Arguments------------------------------------------------------------*/
  M0=pow(10.,logM0) ;
  M1=pow(10.,logM1) ;


/*---Read-BGC-header-----------------------------------------------------------*/

  fp1=fopen(bgc_file,"r") ;
  assert(fp1 != NULL);
  bgc_read_header(fp1,&hdr);
  fclose(fp1) ;
  Nhalos=hdr.ngroups_total ;
  mp=hdr.part_mass*1e10 ;
  Lbox=hdr.BoxSize ;

  Ngaltarget = rhogal*Lbox*Lbox*Lbox ;

/*---Read-halo-center-file-and-build-mass-function-----------------------------*/

  logMbin1 = 10 ;
  logMbin2 = 16 ;
  dlogMbin = 0.01 ;
  Nbin = (int)((logMbin2-logMbin1)/dlogMbin) ;
  Nhalosbin=(int *) calloc(Nbin,sizeof(int)) ;

  fp2=fopen(halos_file,"r") ;
  assert(fp2 != NULL);



  for(i=0;i<Nhalos;i++)
    {
      fscanf(fp2,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %d",&junki,&Npart,&junkf,&junkf,&junkf,&junkf,&junkf,&junkf,&junkf,&junkf,&junki) ;

/*---Apply-Warren-et-al-correction-to-mass---*/

      Mhalo = mp*((double)Npart*(1.-pow(((double)Npart),-0.6))) ;
      logMhalo = log10(Mhalo) ;

/*---Build-mass-function---*/
      
      j = (int)((logMhalo - logMbin1)/dlogMbin) ;

      Nhalosbin[j]++ ;
    }

fclose(fp2);

    for(i=0;i<Nbin;i++){
	logMh = logMbin1 + i*dlogMbin + 0.5*dlogMbin ;
	

	fprintf(stdout,"%lf %d\n",logMh,Nhalosbin[i]);

	}

/*---Loop-over-Mmin-values-------------------------------------------------------*/

  logMmin = 11. ;
  fNgal = 100. ;

  while(fNgal>1.01)
    {
      logMmin += 0.01 ;
      Ngalaxies = 0. ;
      
      for(j=0;j<Nbin;j++)
	{
	  if(Nhalosbin[j]>0)
	    {
	      logMh = logMbin1 + j*dlogMbin + 0.5*dlogMbin ;
	      Mh = pow(10.,logMh) ;
	
/*---Navg(M)--------------------------------------------------------------------*/

	      Ncenavg = 0.5*(1 + erff((logMh-logMmin)/siglogM)) ;
	      
	      if(Mh>M0)
		{
		  Nsatavg = 0.5*(1 + erff((logMh-logMmin)/siglogM))*pow(((Mh-M0)/M1),alpha) ;
		}
	      else
		{
		  Nsatavg = 0. ;
		}
	      
	      Ngalaxies += (Ncenavg + Nsatavg)*(double)Nhalosbin[j] ;
	    }
	}
      fNgal = Ngalaxies/Ngaltarget ;
    }


  return logMmin ;
}






