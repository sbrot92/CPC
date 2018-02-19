//
//  System.c
//  CPC
//
//  Created by Sarah Brotman on 05/01/2018.
//  Copyright © 2018 Sarah Brotman. All rights reserved.
//

#include "System.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include "IO.h"


void init(System *system){
    
    // Read input files and assign necessary variables
    initSystem(system);
    
    // Calculate element volume
    // currently hard coded with Vincen't given input file volume
    // donc using a cylindrical chamber/pellet
    system->V_chamber = PI*pow(.002,2.0)*.0007;
    //printf("V_chamber = %e\n",system->V_chamber);

    
    // Use stoichiometric ratio to determine an acceptable nanoparticle size for CuO until
    // calculated like we talked about with n CuO particles per Al
    
    double C_Al,C_CuO,C_Al2O3;
    // Initial quantities for one mole of thermite according to stoich
    if (system->stoich > 0.0) {
        C_Al = (1.0/(1.0 + 3.0/(2.0 * system->stoich)));
        C_CuO = (3.0/(2.0 * system->stoich)) * C_Al;
    } else if (system->stoich == 0) {
        C_Al = 1.0;
        C_CuO = 0.0;
    }
    
    C_Al2O3 = pow((system->compo[4].a/system->compo[3].a),3.0) * (pow(((system->Ro_Al2O3/system->Ro_Al)+1.0),3.0) - 1) * C_Al;
    C_Al -= C_Al2O3;
    //printf("C_aluminum = %e\nC_cuo = %e\nC_alumina = %e\n",C_Al,C_CuO,C_Al2O3);

    system->compo[0].molarRatio = C_CuO;
    system->compo[3].molarRatio = C_Al2O3;
    system->compo[4].molarRatio = C_Al;
    
    // calc volume of one mole of thermite
    system->V_thermite = (4*PI/3)*(C_Al*pow(system->compo[4].a,3.0) + C_CuO*pow(system->compo[0].a,3.0) + C_Al2O3*pow(system->compo[3].a,3.0));
    //printf("V_thermite = %e\n",system->V_thermite);

    system->V_thermite *= Na;
    //printf("Vm_thermite = %e\n",system->V_thermite);
    
    // calculate number of mol necessary for desired TMD (%TMD)
    system->Nmol_thermite = (system->V_chamber/system->V_thermite)*system->TMD;
    //printf("Nmol_thermite = %e\n",system->Nmol_thermite);

    C_Al *= system->Nmol_thermite;
    C_Al2O3 *= system->Nmol_thermite;
    C_CuO *= system->Nmol_thermite;
    //printf("C_aluminum = %e\nC_cuo = %e\nC_alumina = %e\n",C_Al,C_CuO,C_Al2O3);

    system->compo[0].Nmol = C_CuO;
    system->compo[3].Nmol = C_Al2O3;
    system->compo[4].Nmol = C_Al;

    //printf("nmol 1 Alumine = %5.20e\n",system->compo[3].Nmol);
    //printf("nmol 2 Alumine = %5.20e\n",system->compo[4].Nmol);
    //printf("nmol 3 Alumine = %5.20e\n",system->compo[0].Nmol);

    
    // calc number of nanoparticles
    // This pmust be equal for n=1 sintering model
    system->n_Al = (C_Al*Na) * (3.0/(4.0*PI)) * (pow((system->compo[4].a/system->Ro_Al),3.0));
    system->n_CuO = system->n_Al;
    printf("nParticles = %e\n",system->n_Al);
    
    // Use number of nanoparticles to determine a radius of CuO nanoparticle for this stoichiometric ratio
    system->Ro_CuO = (system->Ro_Al) * (pow((C_CuO/C_Al),(1.0/3.0))/(system->compo[4].a/system->compo[0].a));
    //double hey = ((4*PI/3)*system->n_CuO)/(C_CuO*Na);
    //double Ro2 = system->compo[0].a/(pow(hey,(1.0/3.0)));
    //system->Ro_CuO = 120e-9;
    printf("Ro_CuO = %5.20e\n",system->Ro_CuO);
    //printf("Ro_CuO try 2 = %5.20e\n",Ro2);

    
    system->t = 0;
    system->dT = 0;
    
    system->R1 = system->Ro_CuO;
    system->R2 = system->Ro_Al + system->Ro_Al2O3;
    system->d = sqrt(pow(system->R1,2.0) + pow(system->R2,2.0) - (2*system->R1*system->R2*cos(system->theta)));
    
    initInteraction(system);

    // Initial widths for all species
    system->w_CuO = system->Ro_CuO + system->d1;
    system->w_Cu2O = 0;
    system->w_Cu = 0;
    system->w_Al2O3 = system->Ro_Al2O3;
    system->w_Al = system->Ro_Al + system->d2 - system->w_Al2O3;
    
    /********* DEBUGGING *****
    printf("R_CuO = %e\n",system->Ro_CuO);
    printf("d = %e\n",system->d);
    printf("w_CuO = %e\n",system->w_CuO);
    printf("w_Cu2O = %e\n",system->w_Cu2O);
    printf("w_Cu = %e\n",system->w_Cu);
    printf("w_Al2O3 = %e\n",system->w_Al2O3);
    printf("w_Al = %e\n",system->w_Al); */


    // Initial boundary limits based upon species width's
    // INITIAL SPECIES : CuO, Al2O3, Al
    // Origin is center of particle 1 (CuO)
    
    //LEFTMOST BOUNDARY
    double leftBound = system->Ro_CuO * -1 - system->d1;
    //Interface CuO/Cu2O or Cu2O/Cu
    double a0 = leftBound + system->w_CuO;
    // Interface CuO/Cu2O/Cu with Al2O3
    double a1 = a0 + system->w_Cu2O;
    // Interface Al2O3/ Al
    double a2 = a1 + system->w_Al2O3;
    //RIGHTMOST BOUNDARY
    double rightBound = a2 + system->w_Al;

    /********* DEBUGGING ****
    printf("leftBound = %e\n",leftBound);
    printf("a0 = %e\n",a0);
    printf("a1 = %e\n",a1);
    printf("a2=  %e\n",a2);
    printf("rightBound = %e\n",rightBound);
    printf("d1 = %e\n",system->d1);
    printf("d2= %e\n",system->d2);
    ***********************/

    system->leftBound = leftBound;
    system->a0 = a0;
    system->a1 = a1;
    system->a2 = a2;
    system->rightBound = rightBound;
    
    system->Cs_Cu = 0;
    system->Cs_Al2O3 = 0;
    system->C0 = 0;
    system->C1 = 0;
    system->C2 = 0;

    double d2 = system->d2;
    double R2 = system->R2;
    double R1 = system->R1;
    double d1 = system->d1;

    // LOGIC WRITTEN IN RED NOTEBOOK (SCRAP #2)  LE 01/02/18
    // per particle pair
    double V_R2 = (4*PI/3)*pow(R2,3.0) - (PI/3)*pow(R2-d2,2.0)*(3*R2 - (R2-d2));
    double V_core = (4*PI/3)*pow((R2 - a2),3.0) - (PI/3)*pow((R2-d2-a2),2.0)*(3*(R2-a2) - (R2-d2-a2));
    double a = calcInteractionRadius(system,0);
    double b = calcInteractionRadius(system,system->Ro_Al2O3);
    double V_linear = (PI/6)*a2*(3*pow(a,2.0) + 3*pow(b,2.0) + pow(a2,2.0));
    system->V_Al2O3 = V_R2 - V_core + V_linear;

    system->V_CuO = (4*PI/3)*pow(R1,3.0) - (PI/3)*pow(R1-d1,2.0)*(3*R1 - (R1-d1));
    system->V_Al = V_core - V_linear;
    system->V_Cu2O = 0.0;
    system->V_Cu = 0.0;
    /*********************
    printf("a2 = %5.20e\n",a2);
    printf("R2 = %5.20e\n",R2);
    printf("d2 = %5.20e\n",d2);
    printf("V_Al2O3 = %5.20e\n",system->V_Al2O3);
	**************************/

    double nmol2 = (system->V_Al2O3/system->compo[3].Vm)*system->n_Al;
    double nmol2Al = (system->V_Al/system->compo[4].Vm)*system->n_Al;
    double nmol2CuO = (system->V_CuO/system->compo[0].Vm)*system->n_Al;
    //printf("nmol2 = %5.20e\n",nmol2);
    //printf("nmol2Al = %5.20e\n",nmol2Al);

    //printf("nmol2CuO = %5.20e\n",nmol2CuO);

    system->VreleasedCu = 0.0;
    system->VreleasedAl = 0.0;
    system->dCu2O_Nmol = 0.0;
    system->dCu_Nmol = 0.0;
    system->dAl2O3_Nmol = nmol2/system->n_Al;
    //printf("dAl2O3_Nmol = %5.15e\n",system->dAl2O3_Nmol);


    system->currHfus = 0.0;

    system->flagNoMoreAl = 0;
    system->flagNoMoreCuO = 0;
    system->flagNoMoreCu2O = 0;
    system->heatFlag = 0;

    initResults(system);
}

double calcDiffCoeff(double D0, double Ea, double T) {
	//printf("D0 = %e\t Ea = %e\tR= %e\t T = %e\n",D0,Ea,R,T);
	double D = D0*exp(-Ea/(R*T));
	//printf("D = %e\n",D);
    return D;
}

void calcConcentration(System *system) {
    double T = system->T;
    //printf("T here = %e\n",T);
    double Cs_Cu = system->Cs_ratio * ((system->compo[1].Nmol/system->n_Al)/system->V_Cu2O);
    double Cs_Al2O3 = system->Cs_ratio * ((system->compo[3].Nmol/system->n_Al)/system->V_Al2O3);
    double v = system->v;
    double oldFlux = system->flux;
    
    /**************************
    //if ((system->t - system->lastWrite > system->timeToWrite)) {
    	printf("Temp = %e\n",T);
        printf("Nmol Cu2O per particle pair = %5.15e\n",(system->compo[1].Nmol/system->n_Al));
        printf("V_Cu2O = %5.15e\n",system->V_Cu2O);
        printf("Nmol Al2O3 per particle pair = %5.15e\n",(system->compo[3].Nmol/system->n_Al));
        printf("V_Al2O3 = %5.15e\n",system->V_Al2O3);
        printf("Cs_Al2O3 = %5.15e\n",Cs_Al2O3);
        printf("Cs_Cu = %5.15e\n",Cs_Cu);
    //}
    *************/

    double R1 = system->R1;
    double R2 = system->R2;
    double d1 = system->d1;
    double d2 = system->d2;
    double a0  = system->a0;
    double a1 = system->a1;
    double a2 = system->a2;

    // calculate diffusion coefficients D(T) for three relevant species
    double D_Cu2O = calcDiffCoeff(system->compo[1].Do,system->compo[1].Ea,T);
    // double D_Cu = calcDiffCoeff(system->compo[2].Do,system->compo[2].Ea,T);
    double D_Al2O3 = calcDiffCoeff(system->compo[3].Do,system->compo[3].Ea,T);

    system->D_Cu2O = D_Cu2O;
    system->D_Al2O3 = D_Al2O3;

    /********** DEBUGGING   ***

    if (system->t - system->lastWrite > system->timeToWrite) {
        printf("D_Al2O3 = %5.15e\n", system->D_Al2O3);
        printf("D_Cu2O = %5.15e\n", system->D_Cu2O);
    }
    *************/

    double brack = ((1/(R1*D_Cu2O))*log(((R1+d1-a1)/(R1-d1+a1))*((R1-d1+a0)/(R1+d1-a0)))) + ((1/(R1*D_Al2O3))*log(((R1+d1)/(R1-d1))*((R1-d1+a1)/(R1+d1-a1)))) + ((1/(R2*D_Al2O3))*log(((R2-d2+a2)/(R2+d2-a2))*((R2+d2)/(R2-d2))));
    
    /*********** DEBUGGING   ***

    if ((system->t > 3.800988e-04) || (system->t - system->lastWrite > system->timeToWrite)) {
        printf("R1 = %5.15e\n",R1);
        printf("d1 = %5.15e\n",d1);
        printf("a0 = %5.15e\n",a0);
        printf("a1 = %5.15e\n",a1);
        printf("brack1 = %5.15e\n", (1/(R1*D_Cu2O)));
        printf("brack1.5 = %5.15e\n", ((R1+d1-a1)/(R1-d1+a1))*((R1-d1+a0)/(R1+d1-a0)));
        printf("brack2 = %5.15e\n", log(((R1+d1-a1)/(R1-d1+a1))*((R1-d1+a0)/(R1+d1-a0))));
        printf("brack3 = %5.15e\n", (1/(R1*D_Al2O3)));
        printf("brack4 = %5.15e\n", log(((R1+d1)/(R1-d1))*((R1-d1+a1)/(R1+d1-a1))));
        printf("brack5 = %5.15e\n", (1/(R2*D_Al2O3)));
        printf("brack6 = %5.15e\n", log(((R2-d2+a2)/(R2+d2-a2))*((R2+d2)/(R2-d2))));
        printf("brack = %5.15e\n", brack);

    }
     *************/
    
    // Use known Cs to find flux(T)
    double flux;
    if(system->compo[1].Nmol == 0.0) {
    	flux = (Cs_Al2O3*2*PI)/(brack + (2/(v*(pow(R2,2.0)-pow((a2),2.0)))));
    } else {
    	flux = (Cs_Cu*2*PI)/(brack + (2/(v*(pow(R2,2.0)-pow((a2),2.0)))));
    }

    system->dflux = flux - oldFlux;
    /*********** DEBUGGING   ***
    if ((system->t > 3.800988e-04) || (system->t - system->lastWrite > system->timeToWrite)) {

    printf("D_Al2O3 = %5.15e\n", D_Al2O3);
    printf("D_Cu2O = %5.15e\n", D_Cu2O);
    printf("R1 = %5.15e\n",R1);
    printf("d1 = %5.15e\n",d1);
    printf("a1 = %5.15e\n",a1);
    printf("a0 = %5.15e\n",a0);
    printf("R2 = %5.15e\n",R2);
    printf("d2 = %5.15e\n",d2);
    printf("a2 = %5.15e\n",a2);
    printf("brack = %5.15e\n", brack);
    printf("flux  = %5.15e\n", flux);
    }
	*******/


    // Calculate new oxygen concentrations on boundaries between species
    double C0, C1, C2;
    if(system->compo[1].Nmol == 0.0) {
    	Cs_Cu = Cs_Al2O3;
    	C1 = Cs_Al2O3 - (flux/D_Cu2O)*(1/(2*PI*R1))*log(((R1+d1-a1)/(R1+d1-a0))*((R1-d1+a0))/(R1-d1+a1));
    } else {
    	C1 = Cs_Cu - (flux/D_Cu2O)*(1/(2*PI*R1))*log(((R1+d1-a1)/(R1+d1-a0))*((R1-d1+a0))/(R1-d1+a1));
    }
    if (C1 >= Cs_Al2O3) {
    	brack -= ((1/(R1*D_Cu2O))*log(((R1+d1-a1)/(R1-d1+a1))*((R1-d1+a0)/(R1+d1-a0))));
    	flux = (Cs_Al2O3*2*PI)/(brack + (2/(v*(pow(R2,2.0)-pow((a2),2.0)))));
    	C1 = Cs_Al2O3;
    }

    C0 = C1 - (flux/D_Cu2O)*(1/(2*PI*R1))*log(((R1+d1)/(R1+d1-a1))*((R1-d1+a1)/(R1-d1)));

    C2 = C0 -  (flux/D_Al2O3)*(1/(2*PI*R2))*log(((R2-d2+a2)/(R2-d2))*((R2+d2)/(R2+d2-a2)));
    double otherC2 = (flux/(v*PI*(pow(R2,2.0)-pow(a2,2.0))));
    double otherC0 = otherC2 + (flux/D_Al2O3)*(1/(2*PI*R2))*log(((R2-d2+a2)/(R2-d2))*((R2+d2)/(R2+d2-a2)));

    /********** DEBUGGING ****
    if (system->t - system->lastWrite > system->timeToWrite) {
    printf("CsAl2O3 = %5.15e\n",Cs_Al2O3);
    printf("CsCu = %5.15e\n",Cs_Cu);
    printf("flux = %5.15e\n",flux);
    printf("D_Cu2O = %5.15e\n", D_Cu2O);
    printf("R1 = %5.15e\n",R1);
    printf("d1 = %5.15e\n",d1);
    printf("a1 = %5.15e\n",a1);
    printf("a0 = %5.15e\n",a0);
    printf("expression a = %5.15e\n", ((R1-d1+a0)/(R1+d1-a0))*((R1+d1-a1)/(R1-d1+a1)));
    printf("expression b = %5.15e\n", log(((R1-d1+a0)/(R1+d1-a0))*((R1+d1-a1)/(R1-d1+a1))));
    printf("expression c = %5.15e\n", (flux/D_Cu2O));
    printf("expression d = %5.15e\n", (1/(2*PI*R1)));
    printf("C1 = %5.15e\n",C1);
    printf("C0 = %5.15e\n",C0);
    printf("expression 1 = %5.15e\n", (flux/D_Al2O3) );
    printf("expression 2 = %5.15e\n", (1/(2*PI*R2)) );
    printf("expression 3 = %5.15e\n", ((R2-d2+a2)));
    printf("expression 4 = %5.15e\n", (R2+d2-a2));
    printf("expression 5 = %5.15e\n", (R2+d2));
    printf("expression 6 = %5.15e\n", (R2-d2));
    printf("expression 7 = %5.15e\n", ((R2-d2+a2)/(R2-d2)));
    printf("expression 8 = %5.15e\n", ((R2+d2)/(R2+d2-a2)));
    printf("expression 9 = %5.15e\n", log(((R2-d2+a2)/(R2-d2))*((R2+d2)/(R2+d2-a2))));
    printf("expression 10 = %5.15e\n", (1/(2*PI*R2))*log(((R2-d2+a2)/(R2-d2))*((R2+d2)/(R2+d2-a2))));
    printf("expression 11 = %5.15e\n", (flux/D_Al2O3)*(1/(2*PI*R2))*log(((R2-d2+a2)/(R2-d2))*((R2+d2)/(R2+d2-a2))));
    printf("C2 = %5.15e\n",C2);
    printf("otherC2 = %5.15e\n",otherC2);
    printf("otherC0 = %5.15e\n",otherC0);
    printf("brack' = %5.15e\n", brack);
    printf("flux'  = %5.15e\n", flux);
    }
    *************/
    
    system->flux = flux;
    system->C0 = C0;
    system->C1 = C1;
    system->C2 = C2;
    system->Cs_Cu = Cs_Cu;
    system->Cs_Al2O3 = Cs_Al2O3;

    
}


void calcMassTrans(System *system) {
    double Vm_CuO = system->compo[0].Vm;
    double Vm_Cu2O = system->compo[1].Vm;
    double Vm_Cu = system->compo [2].Vm;
    double Vm_Al2O3 = system->compo[3].Vm;
    double Vm_Al = system->compo[4].Vm;
    double Vm_O2 = system->compo[5].Vm;
    
    double flux = system->flux;
    double R1 = system->R1;
    double R2 = system->R2;
    double d1 = system->d1;
    double d2 = system->d2;
    double a0 = system->a0;
    double a1 = system->a1;
    double a2 = system->a2;
    
    double da0dt, da2dt;
    double VreleasedCu,VreleasedAl;

    double olda0 = a0;
    double olda2 = a2;

    if (system->flagNoMoreCuO == 0) {
    	 da0dt = (flux*Vm_CuO)/(PI*(pow(R1,2.0) - pow((d1-a0),2.0)));
    	 da2dt = (flux*Vm_Al)/(PI*(pow(R2,2.0) - pow((d2-a2),2.0)));

    	 a0 += da0dt*system->dt;
    	 a2 += da2dt*system->dt;


    	 // dV on one particle pair
    	 VreleasedCu = calc_dV(system,(-1*olda0),(-1*a0));
    	 //printf("VreleasedCu = %e", VreleasedCu);
    	 VreleasedAl = calc_dV(system,olda2,a2);

    	 // dNmol
    	 double molReleasedCuO = VreleasedCu/Vm_CuO;
    	 double molCreatedCu2O = molReleasedCuO*.5;
    	 double molCreatedO2 = molReleasedCuO*.25;
    	 double molReleasedAl = VreleasedAl/Vm_Al;
    	 double molCreatedAl2O3 = molReleasedAl*.5;
    	 double molUsedO2 = molReleasedAl*.75;


    	 // by particle pair
    	 system->dCuO_Nmol = molReleasedCuO;
    	 system->dCu2O_Nmol = molCreatedCu2O;
    	 system->dAl2O3_Nmol = molCreatedAl2O3;
    	 system->dAl_Nmol = molReleasedAl;
    	 system->dO2_Nmol = (molCreatedO2 - molUsedO2);

    	 // total system
    	 system->compo[0].Nmol -= (molReleasedCuO * system->n_Al);
    	 system->compo[1].Nmol += (molCreatedCu2O * system->n_Al);
    	 system->compo[3].Nmol += (molCreatedAl2O3 * system->n_Al);
    	 system->compo[4].Nmol -= (molReleasedAl * system->n_Al);
    	 system->compo[5].Nmol += ((molCreatedO2 - molUsedO2) *  system->n_Al);

    	 system->Nmol_O2used = molUsedO2;

    	 system->dV_Cudt = VreleasedCu;
    	 system->V_CuO -= VreleasedCu;
    	 system->V_Cu2O += molCreatedCu2O * Vm_Cu2O;
    	 system->V_Al2O3 += molCreatedAl2O3 * Vm_Al2O3;
    	 system->V_Al -= molReleasedAl;
    } else {
   	 da0dt = (flux*Vm_Cu2O)/(PI*(pow(R1,2.0) - pow((d1-a0),2.0)));
   	 da2dt = (flux*Vm_Al)/(PI*(pow(R2,2.0) - pow((d2-a2),2.0)));

   	 a0 += da0dt*system->dt;
   	 a2 += da2dt*system->dt;


   	 // dV on one particle pair
   	 VreleasedCu = calc_dV(system,(-1*olda0),(-1*a0));
   	 //printf("VreleasedCu = %e", VreleasedCu);
   	 VreleasedAl = calc_dV(system,olda2,a2);

   	 // dNmol
   	 double molReleasedCu2O = VreleasedCu/Vm_Cu2O;
   	 double molCreatedCu = molReleasedCu2O*2;
   	 double molCreatedO2 = molReleasedCu2O*.5;
   	 double molReleasedAl = VreleasedAl/Vm_Al;
   	 double molCreatedAl2O3 = molReleasedAl*.5;
   	 double molUsedO2 = molReleasedAl*.75;


   	 // by particle pair
   	 system->dCu2O_Nmol = molReleasedCu2O;
   	 system->dCu_Nmol = molCreatedCu;
   	 system->dAl2O3_Nmol = molCreatedAl2O3;
   	 system->dAl_Nmol = molReleasedAl;
   	 system->dO2_Nmol = (molCreatedO2 - molUsedO2);

   	 // total system
   	 system->compo[1].Nmol -= (molReleasedCu2O * system->n_Al);
   	 system->compo[2].Nmol += (molCreatedCu * system->n_Al);
   	 system->compo[3].Nmol += (molCreatedAl2O3 * system->n_Al);
   	 system->compo[4].Nmol -= (molReleasedAl * system->n_Al);
   	 system->compo[5].Nmol += ((molCreatedO2 - molUsedO2) *  system->n_Al);

   	 system->Nmol_O2used = molUsedO2;

   	 system->dV_Cudt = VreleasedCu;
   	 system->V_Cu2O -= VreleasedCu;
   	 system->V_Cu += molCreatedCu * Vm_Cu;
   	 system->V_Al2O3 += molCreatedAl2O3 * Vm_Al2O3;
   	 system->V_Al -= molReleasedAl;
    }
    
    /******* DEBUGGING***
    if (a0 > 1.0E-7) {

    //printf("Vm_CuO = %e\n",Vm_CuO);
    //printf("Vm_Al = %e\n",Vm_Al);
    printf("flux = %e\n",flux);
    printf("a0 = %e\n",a0);
    printf("a2 = %e\n",a2);
    printf("R1 = %e\n",R1);
    printf("R2 = %e\n", R2);
    printf("d1  = %e\n", d1);
    printf("d2 = %e\n",d2);
    printf("da0dt = %e\n", da0dt);
    printf("da2dt  = %e\n", da2dt);
    printf("PI  = %e\n", PI);
    printf("dt  = %e\n", system->dt);
    }
	***************/
    
    if (system->flagNoMoreCuO == 0) {
    	if (a0 * -1 < system->leftBound) {
        	printf("Time %e : NO MORE CuO due to a0\n",system->t);
        	system->flagNoMoreCuO = 1;
    	} else if (system->compo[0].Nmol < 1e-7) {
        	printf("Time %e : NO MORE CuO due to nMol\n",system->t);
        	system->flagNoMoreCuO = 1;
    	} else if (system->V_CuO < 0) {
        	printf("Time %e : NO MORE CuO due to V_CuO \n",system->t);
        	system->flagNoMoreCuO = 1;
    	}
    } else if (system->flagNoMoreAl == 0) {
    	if (a2 > system->rightBound) {
        	printf("Time %e : NO MORE Al due to a2\n",system->t);
        	system->flagNoMoreAl = 1;
    	} else if (system->compo[4].Nmol < 1e-7) {
        	printf("Time %e : NO MORE Al due to nMol\n",system->t);
        	system->flagNoMoreAl = 1;
    	} else if (system->V_Al < 0) {
        	printf("Time %e : NO MORE Al due to V_Al \n",system->t);
        	system->flagNoMoreAl = 1;
    	}
    }




    /******* DEBUGGING***

    if (system->t - system->lastWrite > system->timeToWrite) {
        printf("a0 = %5.20e\n",a0);
        printf("olda0 = %5.20e\n",olda0);
        printf("a2 = %5.20e\n",a2);
        printf("olda2 = %5.20e\n",olda2);
        printf("R = %5.20e\n",R1);
        printf("VreleasedCu2 = %5.20e\n",VreleasedCu2);
        printf("VreleasedAl2 = %5.20e\n",VreleasedAl2);
        printf("VreleasedCu = %5.20e\n",VreleasedCu);
        printf("VreleasedAl = %5.20e\n",VreleasedAl);
    }
	***************/

    /******************************
    if (a0 > 1.0E-7) {
    printf("VreleasedCu = %5.20e\n",VreleasedCu);
    printf("VreleasedAl = %5.20e\n",VreleasedAl);
    //printf("molReleasedCuO = %e\n",molReleasedCuO);
    //printf("molCreatedCu2O = %e\n",molCreatedCu2O);
    //printf("molCreatedO2 = %e\n",molCreatedO2);
    //printf("molUsedO2 = %e\n",molUsedO2);
    printf("flux = %5.20e\n",flux);
    //printf("otherUSED %e\n",flux*system->dt);
    printf("molReleasedAl = %5.20e\n",molReleasedAl);
    printf("molCreateAl2O3 = %5.20e\n",molCreateAl2O3);
    printf("molCreatedAl2O3 = %5.20e\n",molCreatedAl2O3);
    }
    ****************************/


}
void calcHeatTrans(System *system) {
    double totalHfus = system->currHfus;
    double T = system->T;
    double flux = system->flux;
    

    /******* DEBUGGING***
    if (system->dT < 0) {

      printf("T = %e\n",T);
      printf("flux = %e\n",flux);
      printf("before totalHfus = %e\n",totalHfus);

    }
	**************/

    // check if melting has occurred for each species
    Species * compo = system->compo;
    //printf("sizeof(&compo) = %lu",sizeof(compo));
    //printf("sizeof(compo[0]) = %lu",sizeof(compo[0]));
    //int n = sizeof(&compo)/sizeof(compo[0]);
    int n = 6;
    for (int i = 0; i < n; i++) {
        if (T >= compo[i].Tfus && compo[i].flagMelt == 0) {
        	compo[i].flagMelt = 1;
        	totalHfus += (compo[i].Hfus * compo[i].Nmol);
        	/*if (T >= 933) {
        		printf("compo[%d] melted, temp %e, delH %e \n", i, compo[i].Tfus, compo[i].Hfus);
        		printf("compo[%d] melted, qfus = %e for %e mol\n", i, compo[i].Hfus * compo[i].Nmol,compo[i].Nmol);
                printf("totalHfus = %e\n",totalHfus);
        	}*/
        }
    }
    //printf("totalHfus = %e\n",totalHfus);


    // calculate heat of reaction
    double qDecomp;

    if (compo[2].Nmol == 0.0) {
    	qDecomp = (2*compo[1].Hf - 4*compo[0].Hf);
    	//q = (delH/2) * system->dCu2O_Nmol;
    } else {
    	qDecomp = (-2*compo[1].Hf);
    	//q = (delH/4) * system->dCu_Nmol;
    }


    double qOxide = (2 * compo[3].Hf);


    //qTotal is heat produced per mol reaction per particle pair due to oxidation (exo) - decomp (endo)
    // 4Al + 3O2 + 12CuO -> 2Al2O3 + 6Cu2O + 3O2 + qTotal
    double dqTotal;
    dqTotal = (-(qOxide + 3*qDecomp));
    // because 12 mol CuO per mol reaction * dnmol_CuO
    dqTotal *= (system->dCuO_Nmol/12);

    double dqTotalSys = dqTotal * system->n_Al;



    //dqTotaldt = ((-(qAlu + 3*delH))/12)*system->dCuO_Nmol;

    //printf("AT TEMP %e : \n",T);
    calcHeatCapacity(system);
    /*for (int i = 0; i < n; i++) {
        printf("calcC of %s = %e\n",compo[i].name,compo[i].Cp);          // TODO TODO TODO ADD Cv
    }*/
    /*calcHeatCapacity(system);
    for (int i = 0; i < n; i++) {
        printf("calcHeatCapacity of %s = %e\n",compo[i].name,compo[i].Cp);          // TODO TODO TODO ADD Cv
    }*/

    // C = Cp*nmol (total in system) for each species
    double C;
    if (compo[2].Nmol == 0.0) {
    	C = compo[0].Cp + compo[1].Cp + compo[3].Cp + compo[4].Cp + compo[5].Cp;
    } else {
    	C = compo[1].Cp + compo[2].Cp + compo[3].Cp + compo[4].Cp + compo[5].Cp;
    }


    /*double C;
    for (int i = 0; i < n; i++) {
        C += compo[i].Cp;          // TODO TODO TODO ADD Cv
    }*/
    //(compo[i].p*compo[i].Vm*
    //printf("C = %e\n",C);
    system->C = C;
    
    /******* DEBUGGING***
    if ((system->t > 3.800988e-04) || (system->t - system->lastWrite > system->timeToWrite)) {

       printf("after totalHfus = %e\n",totalHfus);
       printf("C = %e\n",C);
       printf("qAlu = %e\n",qAlu);

    }
	**************/
    double qTotal;
    if (system->t < system->heatTime) {
    	    qTotal = (system->heatPower*system->dt + dqTotalSys);
    } else {
    		qTotal = dqTotalSys;
    }

    double dT;
    if (totalHfus > 0.0) {
    	//printf("qTotal = %e\n",qTotal);
   	    //printf("totalHfus = %e\n",totalHfus);
     	//printf("newTotalHfus = %e\n",totalHfus - qTotal);
  	   	//printf("T = %e\n",T);
	   	/*for (int i = 0; i < n; i++) {
  	   		printf("compo[%d] is melted? %d\n",i,compo[i].flagMelt);
 	   	}*/
	   	totalHfus -= qTotal;
	    if (system->t - system->lastWrite > system->timeToWrite) {
    		printf("species melting\n %e left\n",totalHfus);
    	}
   		dT = 0.0;
    } else {
		dT = qTotal/C;
    }
    // MULITPLY BY NUMBER OF NANOPARTICLES IN SYSTEM
    //dT *= system->n_Al;

    system->currHfus = totalHfus;

    system->dT = dT;
    T += (dT);
    system->T = T;

    /******* DEBUGGING***
    if (system->t - system->lastWrite > system->timeToWrite) {

    //printf("after totalHfus = %e\n",totalHfus);
    printf("C = %e\n",C);
    if (system->t < system->heatTime) {
    printf("Heat added : %e\n",system->heatPower*system->dt);
    }
    printf("Heat of decomp reaction : %e\n",qDecomp);
    //printf("Heat of decomp reaction per mol Cu2O: %e\n",delH/2);
    printf("mol CuO decomposed: %e\n",system->dCuO_Nmol);
    //printf("q = %e\n",q);
    printf("Heat of  oxy reaction : %e\n",qOxide);
    //printf("Heat of oxy reaction per mol O2 per sec: %e\n",qAlu/3);
    //printf("qAlu per second = %e\n",qA);
    //printf("qAlu = %e\n",qA*system->dt);
    printf("dqTotaldt = %e\n",dqTotal);
    printf("qTotal = %e\n",qTotal);
    printf("dT = %e\n",dT);
    printf("T' = %e\n",T);
    printf("system check :\n");
    printf("dt  = %e\n", system->dt);
    printf("dT = %e\n",system->dT);
    printf("T' = %e\n",system->T);
    }
	***************/
}


void calcHeatCapacity(System * system) {
	double T = system->T;

	// For Al
	double A,B,C,D,E;

	double t = T/1000;

	if (T < 933) {
		A = 28.08920;
		B = -5.414849;
		C = 8.560423;
		D = 3.427370;
		E = -0.277375;
	} else {
		A = 31.75104;
		B = 3.935826E-8;
		C = -1.786515E-8;
		D = 2.694171E-9;
		E = 5.480037E-9;
	}
	system->compo[4].Cp = calcCp(A,B,C,D,E,t);

	// Al2O3
	if (T < 2327) {
		A = 102.4290;
		B = 38.74980;
		C = -15.91090;
		D = 2.628181;
		E = -3.007551;
	} else {
		A = 192.4640;
		B = 9.519856E-8;
		C = -2.858928E-8;
		D = 2.929147E-9;
		E = 5599405E-8;
	}
	system->compo[3].Cp = calcCp(A,B,C,D,E,t);

	// Cu
	if (T < 1358) {
		A = 17.72891;
		B = 28.09870;
		C = -31.25289;
		D = 13.97243;
		E = 0.068611;
	} else {
		A = 32.84450;
		B = -0.000084;
		C = 0.000032;
		D = -0.000004;
		E = -0.000028;
	}
	system->compo[2].Cp = calcCp(A,B,C,D,E,t);

	// Cu2O
	if (T < 1100) {
		A = 59.42033;
		B = 37.84767;
		C = -26.45083;
		D = 11.07609;
		E = -0.542180;
	} else if (T < 1516.7) {
		A = 28.82153;
		B = 39.81603;
		C = 1.459388;
		D = -0.256110;
		E = 11.30291;
	} else {
		A = 95.39980;
		B = 4.892490;
		C = -1.982781;
		D = 0.284974;
		E = 1.522601;
	}
	system->compo[1].Cp = calcCp(A,B,C,D,E,t);

	// CuO
	if (T < 2000) {
		A = 48.56494;
		B = 7.498607;
		C = -0.055980;
		D = 0.013851;
		E = -0.760082;
	} else {
		A = 63.26;
		B = 0;
		C = 0;
		D = 0;
		E = 0;
	}
	system->compo[0].Cp = calcCp(A,B,C,D,E,t);

	if (T < 700) {
		A = 31.32234;
		B = -20.23531;
		C = 57.86644;
		D = -36.50624;
		E = -0.007374;
	} else if (T < 2000) {
		A = 30.03235;
		B = 8.772972;
		C = -3.988133;
		D = 0.788313;
		E = -0.741599;
	} else {
		A = 20.91111;
		B = 10.72071;
		C = -2.020498;
		D = 0.146449;
		E = 9.245722;
	}
	system->compo[5].Cp = calcCp(A,B,C,D,E,t);

	int n = 6;
	for (int i = 0; i < n; i++) {
		system->compo[i].Cp *= system->compo[i].Nmol;
	}

}

double calcCp(double A, double B, double C, double D, double E, double t) {
	double Cp;

	return Cp = A + B*t + C*pow(t,2.0) + D*pow(t,3.0) + E/(pow(t,2.0));
}

double calcInteractionRadius(System * system, double x) {
	double R1 = system->R1;
	double R2 = system->R2;

	double d = system->d;
	double d1 = system->d1;
	double d2 = system->d2;

	double a;

	if (x < 0) {
		a = sqrt(pow((d-x),2.0) - pow(R2,2.0));
	} else if (x == 0) {
		a = sqrt(pow(R1,2.0) - pow(d1,2.0));
	} else if (x > 0 && x < d2) {
		a = sqrt(pow(R1,2.0) - pow(x,2.0));
	}

	//printf("radius = %5.20e\n",a);
	return a;
}

double calc_dV(System * system, double x0, double x1) {
	double a = calcInteractionRadius(system, x1);
	double b = calcInteractionRadius(system, x0);
	double h = fabs(x1 - x0);
	double V;

	V = (1.0/6.0)*PI*h*(3*pow(a,2.0) + 3*pow(b,2.0) + h);
	/**************************
	printf("PI = %5.20e\n",PI);
	printf("x0 = %5.20e\n",x0);
	printf("x1 = %5.20e\n",x1);
	printf("a = %5.20e\n",a);
	printf("b = %5.20e\n",b);
	printf("h = %5.20e\n",h);
	printf("exp1 = %5.20e\n",3*pow(a,2.0));
	printf("exp2 = %5.20e\n",3*pow(b,2.0));
	printf("exp3 = %5.20e\n",(3*pow(a,2.0) + 3*pow(b,2.0) + h));
	printf("exp4 = %5.20e\n",PI*h*(3*pow(a,2.0) + 3*pow(b,2.0) + h));
	printf("exp5 = %5.20e\n",(1.0/6.0)*PI*h*(3*pow(a,2.0) + 3*pow(b,2.0) + h));
	***************************/

	return V;
}
