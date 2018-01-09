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
    
    // Read input files and assing neccessary variables
    initSystem(system);
    
    // Calculate element volume
    // currently hard coded with Vincen't given input file volume
    // donc using a cylindrical chamber/pellet
    system->V_chamber = PI*pow(.002E-3,2.0)*.0007;
    
    // Use stoichiometric ratio to determine an acceptable nanoparticle size for CuO until
    // calculated like we talked about with n CuO particles per Al
    
    // Initial quantities for one mole of themite according to stoich
    if (system->stoich > 0.0) {
        system->C_Al = (1/0/(1/0 + 3.0/(2.0 * system->stoich)));
        system->C_CuO = (3.0/(2.0 * system->stoich)) * system->C_Al;
    } else if (system->stoich == 0) {
        system->C_Al = 1.0;
        system->C_CuO = 0.0;
    }
    
    system->C_Al2O3 = pow((system->compo[4].a/system->compo[3].a),3.0) * (pow(((system->Ro_Al2O3/system->Ro_Al)+1.0),3.0) - 1) * system->C_Al;
    
    // calc volume of one mole of thermite
    system->V_thermite = system->C_Al*pow(system->compo[4].a,3.0) + system->C_CuO*pow(system->compo[0].a,3.0) + system->C_Al2O3*pow(system->compo[3].a,3.0);
    system->V_thermite *= Na;
    
    // calculate number of mol necessary for desired TMD (%TMD)
    system->N_mol = (system->V_chamber/system->V_thermite)*system->TMD;
    
    system->C_Al *= system->N_mol;
    system->C_Al2O3 *= system->N_mol;
    system->C_CuO *= system->N_mol;
    
    // calc number of nanoparticles
    // This pmust be equal for n=1 sintering model
    system->n_Al = (system->C_Al*Na) * (3.0/(4.0*PI)) * (pow((system->compo[4].a/system->Ro_Al),3.0));
    system->n_CuO = system->n_Al;
    
    // Use number of nanoparticles to determine a radius of CuO nanoparticle for this stoichiometric ratio
    system->Ro_CuO = (system->Ro_Al) * (pow((system->C_CuO/system->C_Al),(1.0/3.0))/(system->compo[4].a/system->compo[0].a));
    
    system->t = 0;
    system->dT = 0;
    
    system->R1 = system->Ro_CuO;
    system->R2 = system->Ro_Al;
    system->d = sqrt(pow(system->R1,2.0) + pow(system->R2,2.0) - (2*system->R1*system->R2*cos(system->theta)));
    
    initInteraction(system);
    
    // Initial widths for all species
    system->w_CuO = system->Ro_CuO + system->d1;
    system->w_Cu2O = 0;
    system->w_Cu = 0;
    system->w_Al2O3 = system->Ro_Al2O3;
    system->w_Al = system->Ro_Al + system->d2;
    
    // Initial boundary limits based upon species width's
    // INITIAL SPECIES : CuO, Al2O3, Al
    // Origin is center of particle 1 (CuO)
    
    //LEFTMOST BOUNDARY
    double leftBound = system->Ro_CuO * -1;
    //Interface CuO/Cu2O or Cu2O/Cu
    double a0 = leftBound + system->w_CuO;
    // Interface CuO/Cu2O/Cu with Al2O3
    double a1 = a0 + system->w_Cu2O;
    // Interface Al2O3/ Al
    double a2 = a1 + system->w_Al2O3;
    //RIGHTMOST BOUNDARY
    double rightBound = a2 + system->w_Al;
    
    system->Cs_Cu = 0.1;
    system->Cs_Al2O3 = 0.1;
    system->C0 = 0;
    system->C1 = 0;
    system->C2 = 0;
}

double calcDiffCoeff(double D0, double Ea, double T) {
    return D0*exp(-Ea/R*T);
}

void calcConcentration(System *system) {
    double T = system->T;
    double Cs_Cu = system->Cs_Cu;
    double Cs_Al2O3 = system->Cs_Al2O3;
    double v = system->v;
    
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
    
    double brack = ((1/(R1*D_Cu2O))*log(((R1+d1-a1)/(R1-d1+a1))*((R1-d1+a0)/(R1+d1-a0)))) + ((1/R1*(D_Al2O3))*log(((R1+d1)/(R1-d1))*((R1-d1+a1)/(R1+d1-a1)))) + ((1/(R2*D_Al2O3))*log(((R2-d2+a2)/(R2+d2-a2))*((R2+d2)/(R2-d2))));
    
    // Use known Cs to find flux(T)
    double flux = (Cs_Cu*2*PI)/(brack + (2/(v*(pow(R2,2.0)-pow(a2,2.0)))));
    
    // Calculate new oxygen concentrations on boundaries between species
    double C1 = Cs_Cu - (flux/D_Cu2O)*(1/(2*PI*R1))*log(((R1+d1-a1)/(R1+d1-a0))*((R1-d1+a0))/(R1-d1+a1));
    
    double C0,C2;
    if (C1 < Cs_Al2O3) {
        
        C0 = C1 - (flux/D_Cu2O)*(1/(2*PI*R1))*log(((R1+d1)/(R1+d1-a1))*((R1-d1+a1)/(R1-d1)));
        
        C2 = C0 - (flux/D_Al2O3)*(1/(2*PI*R2))*log(((R2-d2+a2)/(R2-d2))*((R2+d2)/(R2+d2-a2)));
        
    } else {}
    
    system->flux = flux;
    system->C0 = C0;
    system->C1 = C1;
    system->C2 = C2; 
    
}


void calcMassTrans(System *system) {
    double Vm_CuO = system->compo[0].Vm;
    double Vm_Al = system->compo[4].Vm;
    
    double flux = system->flux;
    double R1 = system->R1;
    double R2 = system->R2;
    double d1 = system->d1;
    double d2 = system->d2;
    double a0 = system->a0;
    double a2 = system->a2;
    
    double da0dt = (flux*Vm_CuO)/(PI*(pow(R1,2.0) - pow((d1-a0),2.0)));
    double da2dt = (flux*Vm_Al)/(PI*(pow(R2,2.0) - pow((d2-a2),2.0)));
    
    a0 += da0dt;
    a2 += da2dt;
    
    system->da2dt = da2dt;
    system->da0dt = da0dt;
    system->a0 = a0;
    system->a2 = a2;
}
void calcHeatTrans(System *system) {
    double totalHfus = 0;
    double T = system->T;
    double flux = system->flux;
    
    
    // check if melting has occurred for each species
    Species * compo = system->compo;
    int n = sizeof(compo)/sizeof(compo[0]);
    
    for (int i = 0; i < n; i++) {
        if (T > compo[i].Tfus) {
            totalHfus += compo[i].Hfus;
        }
    }
    
    // calculate heat of reaction
    // q = 
    double q = (compo[2].Hf + compo[3].Hf) - (compo[0].Hf + compo[4].Hf);
    
    double C;
    for (int i = 0; i < n; i++) {
        C += (compo[i].p*compo[i].Vm);          // TODO TODO TODO ADD Cv
    }
    
    double dT = ((flux*q) + totalHfus)/C;
    system->dT = dT;
}
