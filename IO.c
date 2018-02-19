//
//  IO.c
//  CPC
//
//  Created by Sarah Brotman on 05/01/2018.
//  Copyright © 2018 Sarah Brotman. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include "IO.h"


void initSystem(System *system) {
    // SYSTEM PARAMETERS THAT WILL EVENTUALLY BE READ IN FROM INPUT FILE
    system->To = 300;           // deg K
    system->t_final = 3E-2;     // s
    system->TMD = 0.5;          // percentage
    system->stoich = 1;         // stoichiometric ratio of reactants
    system->dt = 1E-12;
    system->v = 3000;           // speed of oxidation reaction for Al (nearly instantaneous, therefore high)
    system->flux = 0;
    system->dflux = 0;
    system->heatTime = 450E-6;
    system->heatPower = 30E3;
    system->timeToWrite = 1E-7;
    system->Cs_ratio = 0.1;
    
    // THERMITE PARAMETERS
    system->Ro_Al = 68.4E-8;         // m
    system->Ro_Al2O3 = 6.6E-9;       // m
    // system->Ro_CuO= 100;
    
    system->D_Cu2O = 0.0;
    system->D_Al2O3 = 0.0;


    // GEOMETRICAL PARAMETERS
    //system->theta = (3/4)*PI;
    system->theta = 2.35619449;
    
    system->compo = malloc(sizeof(Species)*6);
    // SPECIES PARAMETERS AND CONSTANTS
    // make array of composition of thermite, going in expected species order
    // compo[0] : CuO
    // compo[1] : Cu2O
    // compo[2] : Cu
    // compo[3] : Al2O3
    // compo[4] : Al
    Species CuO, Cu2O, Cu, Al2O3, Al, O2;
    
    CuO.name = "CuO";
    Cu2O.name = "Cu2O";
    Cu.name = "Cu";
    Al2O3.name = "Al2O3";
    Al.name = "Al";
    O2.name = "O2";
    
    CuO.a = 2.7556E-10; // m
    Al2O3.a = 3.49E-10; // m
    Al.a = 2.5504E-10;  //m
    
    CuO.M = 79.545;     // g/mol
    Al.M = 26.9815386;         // g/mol
    Al2O3.M = 101.9613;
    Cu2O.M = 143.091;
    Cu.M = 63.546;
    O2.M = 31.998;
    
    // densities
    Al.p = 2.7E6;        // g/m^3
    Al2O3.p = 3.950E6;
    CuO.p = 6.310E6;
    Cu2O.p = 6.000E6;
    Cu.p = 8.960E6;
    O2.p = 1380;
    
    // molar volumes
    Al.Vm = Al.M/Al.p;       // m^3/mol
    Al2O3.Vm = Al2O3.M/Al2O3.p;
    CuO.Vm = CuO.M/CuO.p;
    Cu2O.Vm = Cu2O.M/Cu2O.p;
    Cu.Vm = Cu.M/Cu.p;
    O2.Vm = O2.M/O2.p;
    
    
    //activation energies
    Cu2O.Ea = 6.73E4;
    Cu.Ea = Cu2O.Ea;
    Al2O3.Ea = 1.4E5;
    
    // prefactors
    Cu2O.Do = 1.16E-6;
    Cu.Do = Cu2O.Do;
    Al2O3.Do = 9E-5;
    
    // enthalpies of formations
    // J/mol
    Cu2O.Hf = -168600;
    CuO.Hf = -157300;
    Al2O3.Hf = -1675700;
    Cu.Hf = 0;
    Al.Hf = 0;
    
    // enthalpies of fusion/melting
    // J/mol
    CuO.Hfus = 17.47E3;
    Cu2O.Hfus = 17.47E3;
    Cu.Hfus = 13E3;
    Al2O3.Hfus = 111E3;
    Al.Hfus = 10.79E3;
    O2.Hfus = 0.0;

    // temperatures of melting
    CuO.Tfus = 1599;
    Cu2O.Tfus = 1517;
    Cu.Tfus = 1358;
    Al2O3.Tfus = 2315;
    Al.Tfus = 933;
    O2.Tfus = 0.0;


    CuO.molarRatio = 0.0;
    Cu2O.molarRatio = 0.0;
    Cu.molarRatio = 0.0;
    Al2O3.molarRatio = 0.0;
    Al.molarRatio = 0.0;

    CuO.Nmol = 0.0;
    Cu2O.Nmol = 0.0;
    Cu.Nmol= 0.0;
    Al2O3.Nmol= 0.0;
    Al.Nmol = 0.0;

    CuO.flagMelt = 0;
    Cu2O.flagMelt = 0;
    Cu.flagMelt= 0;
    Al2O3.flagMelt= 0;
    Al.flagMelt = 0;
    O2.flagMelt = 0;

    system->compo[0] = CuO;
    system->compo[1]= Cu2O;
    system->compo[2]= Cu;
    system->compo[3]= Al2O3;
    system->compo[4] = Al;
    system->compo[5] = O2;
}


// PROBABLY NOT NECESSARY

void initInteraction(System *system) {
    
    double R1 = system->R1;
    double R2 = system->R2;
    double d = system->d;
    double theta = system->theta;

    if (pow(d,2.0) > (pow(R1,2.0)+pow(R2,2.0))) {
        system->d1 = (R1*(R1 - R2*cos(theta)))/d;
        system->d2 = (R2*(R2 - R1*cos(theta)))/d; 
    }

    /******* DEBUGGING    ******
    printf("CHECKKKKKKK");
    printf("d1 = %e\n\n\n",system->d1);
    printf("d2 = %e\n\n\n",system->d2);
    printf("d = %e\n\n\n",system->d);
    printf("R2 = %e\n\n\n",system->R2);
    printf("R1 = %e\n\n\n",system->R1);
    printf("theta = %e\n\n\n",system->theta);
    *****************************/

}

void initResults(System * system) {
	char folder[256] = "Results";

	mode_t process_mask = umask(0);
	mkdir(folder, S_IRWXU);
	umask(process_mask);

	FILE * file = NULL;
	file = fopen("Results/Temperature.txt", "w+");

	if (file != NULL) {
		fprintf(file, "Temperature Results\nTime - t (s)\t Temperature - T (K)\n");
		fprintf(file, "%e\t%e\n", system->t, system->To);
	} else {
		printf("Cannot open the file : 'Temperature.txt'\n");
	}

	fclose(file);

	file = fopen("Results/Concentration.txt", "w+");

	if (file != NULL) {
		fprintf(file, "Concentration of Oxygen Diffusion Results\nTime - t (s)\t Concentration - Cs ()\t Concentration - C1 ()\t Concentration - C0 ()\t Concentration - C2()\n");
		fprintf(file, "\t\t Interface CuO/Cu2O \t Interface Cu2O/Al2O3 \t Interface of particles \t Interface Al2O3/Al \n");
		fprintf(file, "%e\t%e\t%e\t%e\t%e\n",system->t,system->Cs_Cu,system->C1,system->C0,system->C2);
	} else {
		printf("Cannot open the file : 'Concentration.txt'");
	}

	fclose(file);

	file = fopen("Results/Composition.txt","w+");

	if (file != NULL) {
		fprintf(file, "Composition of System Results\nOrigin is at center of CuO (left) particle\n");
		fprintf(file, "Time -t (s)\t Left Boundary (m) \t a0 (m) \t a1 (m) (m) \t a2 (m) \t Right Boundary (m)\n");
		fprintf(file, "%e\t%e\t%e\t%e\t%e\t%e\n",system->t,system->leftBound,-1*system->a0,system->a1,system->a2,system->rightBound);
	} else {
		printf("Cannot open the file: 'Composition.txt'");
	}

	fclose(file);

	file = fopen("Results/VolumePerParticle.txt","w+");

		if (file != NULL) {
			fprintf(file, "Composition of System Results\n");
			fprintf(file, "Time -t (s)\t  Vol_CuO\t dVol_CuO\t Vol_Cu2O\t Vol_Cu\t Vol_Al2O3\t Vol_Al\n");
			fprintf(file, "%e\t%e\t%e\t%e\t%e\t%e\n",system->t,system->V_CuO,system->V_Cu2O,system->V_Cu,system->V_Al2O3,system->V_Al);
		} else {
			printf("Cannot open the file: 'VolumePerParticle.txt'");
		}

	fclose(file);

	file = fopen("Results/Flux.txt","w+");
	if(file != NULL) {
		fprintf(file, "Oxygen Flux Results\nTime -t (s)\t Flux - phi ()\t Change in Flux - dPhi\n");
		fprintf(file, "%e\t%e\n",system->t,system->flux);
	} else {
		printf("Cannot open the file: 'Flux.txt'");
	}

	fclose(file);

	file = fopen("Results/Diffusion.txt","w+");
	if(file != NULL) {
		fprintf(file, "Diffusion Coefficient Results\nTime -t (s)\t D_Cu2O ()\t D_Al2O3\n");
		fprintf(file, "%e\t%e\t%e\n",system->t,system->D_Cu2O,system->D_Al2O3);
	} else {
		printf("Cannot open the file: 'Diffusion.txt'");
	}

	fclose(file);

	file = fopen("Results/Quantities Per Particle.txt","w+");
	if(file != NULL) {
		fprintf(file, "Time -t (s)\t CuO\t Cu2O\t Cu\t Al2O3\t Al\t O2\n");
		double molCuO = system->compo[0].Nmol/system->n_Al;
		double molCu2O = system->compo[1].Nmol/system->n_Al;
		double molCu = system->compo[2].Nmol/system->n_Al;
		double molAl2O3 = system->compo[3].Nmol/system->n_Al;
		double molAl = system->compo[4].Nmol/system->n_Al;
		double molO2 = system->compo[5].Nmol/system->n_Al;
		fprintf(file, "%e\t%1.15e\t%1.15e\t%1.15e\t%1.15e\t%1.15e\t%1.15e\n",system->t,molCuO,molCu2O,molCu,molAl2O3,molAl,molO2);
	} else {
		printf("Cannot open the file: 'Diffusion.txt'");
	}

	fclose(file);
}

void writeResults(System * system) {
		FILE * file = NULL;
		file = fopen("Results/Temperature.txt", "a");

		if (file != NULL) {
			fprintf(file, "%e\t%5.25lf\t%e\t%e\n", system->t, system->T,system->dT,system->C);
		} else {
			printf("Cannot open the file : 'Temperature.txt'\n");
		}

		fclose(file);

		file = fopen("Results/Concentration.txt", "a");

		if (file != NULL) {
			fprintf(file, "%e\t%e\t%5.20e\t%5.20e\t%5.20e\n",system->t,system->Cs_Cu,system->C1,system->C0,system->C2);
		} else {
			printf("Cannot open the file : 'Concentration.txt'");
		}

		fclose(file);

		file = fopen("Results/Composition.txt","a");

		if (file != NULL) {
			fprintf(file, "%e\t%e\t%e\t%e\t%e\t%e\n",system->t,system->leftBound,-1*system->a0,system->a1,system->a2,system->rightBound);
		} else {
			printf("Cannot open the file: 'Composition.txt'");
		}
		fclose(file);

		file = fopen("Results/VolumePerParticle.txt","a");

		if (file != NULL) {
			fprintf(file, "%e\t%e\t%e\t%e\t%e\t%e\n",system->t,system->V_CuO,system->V_Cu2O,system->V_Cu,system->V_Al2O3,system->V_Al);
		} else {
			printf("Cannot open the file: 'VolumePerParticle.txt'");
		}
		fclose(file);

		file = fopen("Results/Flux.txt","a");
		if(file != NULL) {
			fprintf(file, "%e\t%e\t%e\n",system->t,system->flux,system->dflux);
		} else {
			printf("Cannot open the file: 'Flux.txt'");
		}

		fclose(file);

		file = fopen("Results/Diffusion.txt","a");
		if(file != NULL) {
			fprintf(file, "%e\t%e\t%e\n",system->t,system->D_Cu2O,system->D_Al2O3);
		} else {
			printf("Cannot open the file: 'Diffusion.txt'");
		}

		fclose(file);

		file = fopen("Results/Quantities Per Particle.txt","a");
		if(file != NULL) {
			double molCuO = system->compo[0].Nmol/system->n_Al;
			double molCu2O = system->compo[1].Nmol/system->n_Al;
			double molCu = system->compo[2].Nmol/system->n_Al;
			double molAl2O3 = system->compo[3].Nmol/system->n_Al;
			double molAl = system->compo[4].Nmol/system->n_Al;
			double molO2 = system->compo[5].Nmol/system->n_Al;
			fprintf(file, "%e\t%1.15e\t%1.15e\t%1.15e\t%1.15e\t%1.15e\t%1.15e\n",system->t,molCuO,molCu2O,molCu,molAl2O3,molAl,molO2);
		} else {
			printf("Cannot open the file: 'Quantities Per Particle.txt'");
		}

		fclose(file);
}

/*************************************/
/* 	int n = sizeof(system.compo/system.compo[0]);

	if (file != NULL) {
		for (int i = 0; i < n; i++) {

		}

		*******************/

