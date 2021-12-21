//
//  Shared_Code_lat_input.c
//
//  Created by boris Sauterey on 05/10/2021.
//  Copyright (c) 2014 boris. All rights reserved.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>


////// Model of trait diffusion ///////
int main(int argc, char *argv[])
{

	omp_set_num_threads(1);

	// Output files
	clock_t begin, intermediate;
	begin = clock();
	double time_spent;
	time_t now;
	double time_b = omp_get_wtime();
	struct tm *today;
	char date[100];
	time(&now);
	today = localtime(&now);
	strftime(date, 100, "S_and_T_%d.%m.%Y-%Hh-%Mmin-%Ssec", today);
	char directory_name[100];
	snprintf(directory_name,100,"S_and_T_lat%d",atoi(argv[1]));
	mkdir(directory_name, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	chdir(directory_name);

	/// OUTPUTS:
	FILE *TIME;
	char filename1[100];
	snprintf(filename1, 100, "time.txt");
	TIME = fopen(filename1,"w");
	if (TIME == NULL){
		printf("Can't open file\n");
		exit(0);
	}

	FILE *ABUNDANCEN;
	char filename2[100];
	snprintf(filename2, 100, "abundancen.txt");
	ABUNDANCEN = fopen(filename2,"w");
	if (ABUNDANCEN == NULL){
		printf("Can't open file\n");
		exit(0);
	}

	FILE *ABUNDANCEP;
	char filename3[100];
	snprintf(filename3, 100, "abundancep.txt");
	ABUNDANCEP = fopen(filename3,"w");
	if (ABUNDANCEP == NULL){
		printf("Can't open file\n");
		exit(0);
	}

	FILE *ABUNDANCEZ;
	char filename4[100];
	snprintf(filename4, 100, "abundancez.txt");
	ABUNDANCEZ = fopen(filename4,"w");
	if (ABUNDANCEZ == NULL){
		printf("Can't open file\n");
		exit(0);
	}

	fclose(TIME);
	fclose(ABUNDANCEN);
 	fclose(ABUNDANCEP);
	fclose(ABUNDANCEZ);

	//Code parameters
	int 	  i;                    	     //Indexes
	int     i2;                     	   //
	int     j;                   	       //
	int     k;                 	         //

	int     Sres     = 200;       	     //Number of size-classes
	int     Spaceres = 1;         	     //Number of patches
	double  sizep[Sres];          	     //Phytoplankton size classes
	double  sizez[Sres];           	     //Zooplankton size classes

	double 	t;
	double  time   = 10000;         	   //Max time
	double  deltat = 0.01;       	       //Time step
	double  deltat_ad;            	     //Adaptive time step
	double  minNP;                 	     //Variable with the lowest X/dX
	int     ti     = 0;            	     //Keeps track of the number of time steps
	int     out    = 400;  	  	    	   //Number of time steps in outputs
	double  prompt = time/out;      	   //Tracks when to write outputs
	int     t_temp = 0;             	   //Time ratchet to signal output writing
	double  dt     = 0.01;            	 //Alternative time step (keeping both is not necessary in that version of the code)

	//State variable
	double   N[Spaceres];                       //Nutrient concentration
	double **P = malloc(Sres*out * sizeof *P);  //Spatial distribution of the phytoplankton size-classes
	for (i = 0; i < Sres*out; i++){P[i] = malloc(Spaceres * sizeof *P[i]);}
	double **Z = malloc(Sres*out * sizeof *Z);  //Spatial distribution of the zooplankton size-classes
	for (i = 0; i < Sres*out; i++){Z[i] = malloc(Spaceres * sizeof *Z[i]);}

	//Environmental and Biological parameters
	double  lat[Spaceres];               //Latitudinal transect
	double  T[Spaceres];                 //Temperature
	double  Tref = 18;           	       //Reference temperature from Marañon et al 2013
	double  N0[Spaceres];        	       //Spatial gradient of nutrient inputT
	double  muzmax[Sres];        	       //Max attack rate of zooplankton
	double **Vmax = malloc(Sres*out * sizeof *Vmax);  //Max uptake rate
	for (i = 0; i < Sres*out; i++){Vmax[i] = malloc(Spaceres * sizeof *Vmax[i]);}
	double  K[Sres];             	       //Half-saturation constant for photosynthesis
	double  Qmin[Sres];          	       //Min internal quota of nutrient
	double **mu = malloc(Sres*out * sizeof *mu);  //Max growth rate
	for (i = 0; i < Sres*out; i++){mu[i] = malloc(Spaceres * sizeof *mu[i]);}
	double **mumax = malloc(Sres*out * sizeof *mumax);  //Composite max growth rate of phytoplankton when solving the equilibrium of the internal nutrient quota
	for (i = 0; i < Sres*out; i++){mumax[i] = malloc(Spaceres * sizeof *mumax[i]);}
	double **Kn = malloc(Sres*out * sizeof *Kn);  //Composite half-saturation constant of phytoplankton when solving the equilibrium of the internal nutrient quota
	for (i = 0; i < Sres*out; i++){Kn[i] = malloc(Spaceres * sizeof *Kn[i]);}
	double  Kg;                  	       //Half-saturation constant for Zooplankton
	double  mortp;                		   //Phytoplankton mortality rate
	double  mortz;                		   //Zooplankton mortality rate

	double  Di;                           //Dilution rate of the chemostat model (I in the manuscript)
	double  stp;                  		    //Size evolutionary diffusion constant
	double  Spacetp;             		      //Spatial diffusion constant
	double **predprey = malloc(Sres*out * sizeof *predprey);  //Size ratio matrix
	for (i = 0; i < Sres*out; i++){predprey[i] = malloc(Sres * sizeof *predprey[i]);}
	double  difterm;             		      //Evolutionary diffusion term
	double  Sdifterm;             		    //Spatial diffusion term for biological populations
	double  Ndifterm;              		    //Spatial diffusion term for nutrients (not in the original code)

	double **pref = malloc(Sres*out * sizeof *pref);  //Size-based grazing kernel
	for (i = 0; i < Sres*out; i++){pref[i] = malloc(Sres * sizeof *pref[i]);}
	double  conso;                        //Total nutrient consumption of the local phytoplankton community
	double  sum_graz[Sres];        		    //Total prey availability for a specific zooplankton size-class (in a given patch)
	double  graz_eff[Sres];        		    //Total grazing effort of a specific zooplankton size-class (in a given patch)
	double **graz_mat = malloc(Sres*out * sizeof *graz_mat);  //Grazing matrix
	for (i = 0; i < Sres*out; i++){graz_mat[i] = malloc(Sres * sizeof *graz_mat[i]);}
	double  consoz;                		    //Loss term for a phytoplankton size-class due to grazing
	double  consop;                		    //Growth term for a zooplankton size-class from grazing

	double  dN[Spaceres];          		  //Derivative for nutrient
	double **dP = malloc(Sres*out * sizeof *dP);  //Derivative for phytoplankton abundances
	for (i = 0; i < Sres*out; i++){dP[i] = malloc(Spaceres * sizeof *dP[i]);}
	double **dZ = malloc(Sres*out * sizeof *dZ);  //Derivative for zooplankton abundances*/
	for (i = 0; i < Sres*out; i++){dZ[i] = malloc(Spaceres * sizeof *dZ[i]);}

//DEBUG
	unsigned int    i_min = 0;
	unsigned int    j_min = 0;

	/////* Size Vectors */////////////////////////////////////////////////////////////////////
	for (i=0; i<Sres; i++)
	{
		sizep[i] = exp(log(1e-5)+(double)i*(log(1e6)-log(1e-5))/(Sres-1));   //Volume in µm-3
		sizez[i] = exp(log(1e-3)+(double)i*(log(1e8)-log(1e-3))/(Sres-1));   //Volume in µm-3
	}

	/////* Size ratio and predation pref matrixes *///////////////////////////////////////////
	for(i=0; i<Sres; i++)
	{
		for(j=0; j<Sres; j++)
		{
			predprey[i][j]=sizez[i]/sizep[j];
			pref[i][j]=1.0/(2.0*0.5)*exp(-pow(log(predprey[i][j]/100.0),2.0)
                             /(2.0*pow(0.5,2.0)));
		}
	}

	/////* Other parameters */////////////////////////////////////////////////////////////////


	Di = 0.1;  /* Dilution rate */
	stp     = 5e-4*pow(((float)(Sres)/200),4);  /* Diffusion rates */
	Spacetp = 0;

	/////* Initial conditions *///////////////////////////////////////////////////////////////
	t = 0.0; /* Time */
	for (i=0; i<Spaceres; i++) /* N & T */
	{
		lat[i] = atoi(argv[1]);
		N0[i]  = 1e-1+10/(1+exp(-0.12*(lat[i]-50))) + 2e0*exp(-0.08*lat[i]);
		T[i]   = 40-lat[i]/80*25;
		N[i]   = N0[i];
	}

	for(i=0; i<Sres; i++) /* Physiological parameter for phytos and zoos */
	{
		Qmin[i]   = 0.14*pow(sizep[i],-0.04);                     //Marañon et al 2013 (mmol N mmol C-1)
		K[i]      = 0.15*pow(sizep[i],0.33);                      //Edward et al 2010  (mmol m-3)
		for(j=0; j<Spaceres; j++)
		{
			Vmax[i][j]   = 0.10*exp(0.063*(T[j]-Tref))*pow(sizep[i],0.09);  //Marañon et al 2013 (mmol N mmol C-1 d-1)
			mu[i][j]     = 6.40*exp(0.063*(T[j]-Tref))*pow(sizep[i],-0.27); //Ward et al 2017    (d-1)
			mumax[i][j]  = mu[i][j]*Vmax[i][j]/(Vmax[i][j]+mu[i][j]*Qmin[i]);
			Kn[i][j]     = mu[i][j]*K[i]*Qmin[i]/(Vmax[i][j]+mu[i][j]*Qmin[i]);
		}

		muzmax[i] = 21.9*pow(sizez[i],-0.16);                   // Hansen et al 97
	}
	mortp = 0.15;
	Kg    = 2;
	mortz = 0.15;

	for(i=0; i<Sres; i++) /* Phytos */
	{
		for (j=0; j<Spaceres; j++)
		{
			if(i==(int)floor((Sres-1)/2))
			{
				P[i][j] = 0.1;
			}
			else
			{
				P[i][j] = 1e-20;
			}
		}
	}

	for(i=0; i<Sres; i++) /* Zoos */
	{
		for (j=0; j<Spaceres; j++)
		{
			if(i==(int)floor((Sres-1)/2))
			{
				Z[i][j] = 0.1;
			}
			else
			{
				Z[i][j] = 1e-20;
			}
		}

	}

	////////////////////////// LOOP //////////////////////////////////////////////////////////
	while(t <= time)
	{

		////* Derivative */
		for (j=0; j<Spaceres; j++)
		{
			conso = 0.0;
			for(i=0; i<Sres; i++)
			{
				/*pop i*/ /*upt. i*/  /*sat. i*/
				conso += P[i][j] * Vmax[i][j] * N[j]/(N[j]+Kn[i][j]);
			}

			/* Spatial diffusion term for nutrient */
			if(j==0)
			{
				//Ndifterm = Spacetp * (N[j+1] - N[j]);
				Ndifterm = 0;
			}
			else if(j==Spaceres-1)
			{
				//Ndifterm = Spacetp * (N[j-1] - N[j]);
				Ndifterm = 0;
			}
			else
			{
				//Ndifterm = Spacetp * (N[j-1] + N[j+1] - 2.0*N[j]);
				Ndifterm = 0;
			}

			dN[j]=Di*(N0[j]-N[j])-conso+Ndifterm;

			///// Grazing matrix /////
			/* Sum of preys weighted by pref function and saturation function */
			for(i=0; i<Sres; i++) /* pred loop */
			{
				sum_graz[i] = 0;
				for(i2=0; i2<Sres; i2++) /* prey loop */
				{
					sum_graz[i] += pref[i][i2]*P[i2][j];
				}
				graz_eff[i] = sum_graz[i]/(sum_graz[i]+Kg)*muzmax[i]*(1.0-exp(-1.0*sum_graz[i]));
			}
			/* Grazing matrix */
			for(i=0; i<Sres; i++)
			{
				for(i2=0; i2<Sres; i2++)
				{
					if (sum_graz[i] > 0)
					{
						graz_mat[i][i2] = graz_eff[i]*pref[i][i2]/sum_graz[i];
					}
					else
					{
						graz_mat[i][i2] = 0;
					}
				}
			}

			//////////////////////////

			/* Phytos */


			for(i=0; i<Sres; i++)
			{
				/* Sum over zoos of the grazing */
				consoz = 0.0;

				for(i2=0; i2<Sres; i2++)
				{
					consoz = consoz + graz_mat[i2][i]*Z[i2][j];
				}
				/* Trait diffusion term */
				if(i==(int)(0))
				{
					difterm = stp * (P[i+1][j] - 2.0*P[i][j]);
				}
				else if(i==(int)(Sres-1))
				{
					difterm = stp * (P[i-1][j] - 2.0*P[i][j]);
				}
				else
				{
					difterm = stp * (P[i-1][j] + P[i+1][j] - 2.0*P[i][j]);
				}

				/* Spatial diffusion term */
				if(j==(int)(0))
				{
					//Sdifterm = Spacetp * (P[i][j+1] - P[i][j]);
					Sdifterm = 0;
				}
				else if(j==(int)(Spaceres-1))
				{
					//Sdifterm = Spacetp * (P[i][j-1] - P[i][j]);
					Sdifterm = 0;
				}
				else
				{
					//Sdifterm = Spacetp * (P[i][j-1] + P[i][j+1] - 2.0*P[i][j]);
					Sdifterm = 0;
				}

								   /*	   conso	   */                  /*		 mortality		  */
				dP[i][j] = (mumax[i][j]*N[j]/(N[j]+Kn[i][j]) - (1*exp(-1.0e10/Sres*1e2*P[i][j])+mortp)
					                  /*grazing*/
                    - consoz) * P[i][j] + difterm + Sdifterm;
			}

			/* Zoos */
			for(i=0; i<Sres; i++)
			{
				/* Sum over phytos of the grazing */
				consop = 0.0;
				for(i2=0; i2<Sres; i2++)
				{
					consop = consop + graz_mat[i][i2]*P[i2][j];
				}

				/* Trait diffusion term */
				if(i==0)
				{
					difterm = stp * (Z[i+1][j] - 2.0*Z[i][j]);
				}
				else if(i==Sres-1)
				{
					difterm = stp * (Z[i-1][j] - 2.0*Z[i][j]);
				}
				else
				{
					difterm = stp * (Z[i-1][j] + Z[i+1][j] - 2.0*Z[i][j]);
				}

				/* Spatial diffusion term */
				if(j==0)
				{
					//Sdifterm = Spacetp * (Z[i][j+1] - Z[i][j]);
					Sdifterm = 0;
				}
				else if(j==Spaceres-1)
				{
					//Sdifterm = Spacetp * (Z[i][j-1] - Z[i][j]);
					Sdifterm = 0;
				}
				else
				{
					//Sdifterm = Spacetp * (Z[i][j-1] + Z[i][j+1] - 2.0*Z[i][j]);
					Sdifterm = 0;
				}
				           /*grazing*/ /*		   mortality			*/
				dZ[i][j] = (consop - (1*exp(-1.0e10/Sres*1e2*Z[i][j])+mortz)) * Z[i][j] + difterm + Sdifterm;
			}

		}

		minNP = fabs(N[0]/dN[0]);
		for(i=Spaceres; i--; )
		{
			if(fabs(N[i]/dN[i]) < minNP)
			minNP = fabs(N[i]/dN[i]);
		}

		  for(i=0; i<Sres; i++)
		  {
				for(j=0; j<Spaceres; j++)
				{
					if(P[i][j] < 1e-20 & dP[i][j] < 0)
					{
						dP[i][j] = 0;
					}
					if(Z[i][j] < 1e-20 & dZ[i][j] < 0)
					{
						dZ[i][j] = 0;
					}
					if(P[i][j] != 0.0 && dP[i][j] != 0.0 && fabs(P[i][j]/dP[i][j]) < minNP)
					{
						minNP = fabs(P[i][j]/dP[i][j]);
						i_min = i;
						j_min = j;
					}
					if(Z[i][j] != 0.0 && dZ[i][j] != 0.0 && fabs(Z[i][j]/dZ[i][j]) < minNP)
					{
						minNP = fabs(Z[i][j]/dZ[i][j]);
						i_min = i;
						j_min = j;
					}
				}
		  }

		deltat_ad = minNP/10; /* Time step */
		dt = deltat_ad;

		////* Integration */

		t = t+dt; /* Time */

		for (j=0; j<Spaceres; j++)
		{
			N[j] = N[j]+dN[j]*dt; /* N */
			if (N[j] != N[j] || N[j] < 0)
			{
				printf("No more nutrients...\n");
				exit(0);
			}

			for(i=0; i<Sres; i++) /* Phytos and zoos */
			{
				P[i][j] = P[i][j] + dP[i][j]*dt;
				Z[i][j] = Z[i][j] + dZ[i][j]*dt;

				if (P[i][j] != P[i][j] || P[i][j] < 0)
				{
					printf("Crash P[%i][%i] = %e\n",i,j,P[i][j]);
					exit(0);
				}

				if (Z[i][j] != Z[i][j] || Z[i][j] < 0)
				{
					printf("Crash Z[%i][%i] = %e\n",i,j,Z[i][j]);
					exit(0);
				}
			}
		}


		/* prompt */
		if((int)ceil(t/prompt) > t_temp)
		{

			printf("t = %f, or %f %%\n",t,t/time*100);

			TIME = fopen(filename1,"a+");
			ABUNDANCEN = fopen(filename2,"a+");
			ABUNDANCEP = fopen(filename3,"a+");
			ABUNDANCEZ = fopen(filename4,"a+");

			for (j=0; j<Spaceres; j++)
			{
				fprintf(ABUNDANCEN,"%e\n",N[j]);
				for(i=0; i<Sres; i++)
				{
					if (i < (Sres-1))
					{
						fprintf(ABUNDANCEP,"%e ",P[i][j]);
						fprintf(ABUNDANCEZ,"%e ",Z[i][j]);
					}
					else
					{
						fprintf(ABUNDANCEP,"%e\n",P[i][j]);
						fprintf(ABUNDANCEZ,"%e\n",Z[i][j]);
					}
				}
			}
			fprintf(TIME,"%f\n",t);

			fclose(TIME);
			fclose(ABUNDANCEN);
			fclose(ABUNDANCEP);
			fclose(ABUNDANCEZ);

			t_temp += 1;
		}

		ti = ti + 1;
	}

	return 0;
}
