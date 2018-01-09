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
    system->dt = 1E-20;
    system->v = 3000;           // speed of oxidation reaction for Al (nearly instantaneous, therefore high)

    
    // THERMITE PARAMETERS
    system->Ro_Al = 8E-8;         // m
    system->Ro_Al2O3 = 2E-9;       // m
    // system->Ro_CuO= 100;
    
    // GEOMETRICAL PARAMETERS
    system->theta = (3/2)*PI;
    
    
    // SPECIES PARAMETERS AND CONSTANTS
    // make array of composition of thermite, going in expected species order
    // compo[0] : CuO
    // compo[1] : Cu2O
    // compo[2] : Cu
    // compo[3] : Al2O3
    // compo[4] : Al
    system->compo = malloc(sizeof(Species)*5);
    
    Species CuO, Cu2O, Cu, Al2O3, Al;
    
    CuO.name = "CuO";
    Cu2O.name = "Cu2O";
    Cu.name = "Cu";
    Al2O3.name = "Al2O3";
    Al.name = "Al";
    
    CuO.a = 2.7556E-10; // m
    Al2O3.a = 3.49E-10; // m
    Al.a = 2.5504E-10;  //m
    
    CuO.M = 79.545;     // g/mol
    Al.M = 26.9815386;         // g/mol
    Al2O3.M = 101.9613;
    Cu2O.M = 143.091;
    Cu.M = 63.546;
    
    // densities
    Al.p = 2.7;        // g/m^3
    Al2O3.p = 3.950;
    CuO.p = 6.310;
    Cu2O.p = 6.000;
    Cu.p = 8.960;
    
    // molar volumes
    Al.Vm = Al.M/Al.p;       // m^3/mol
    Al2O3.Vm = Al2O3.M/Al2O3.p;
    CuO.Vm = CuO.M/CuO.p;
    Cu2O.Vm = Cu2O.M/Cu2O.p;
    Cu.Vm = Cu.M/Cu.p;
    
    
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
    Cu.Hf = -11860;
    Al.Hf = -105600;
    
    // enthalpies of fusion/melting
    // J/mol
    CuO.Hfus = 17.47E3;
    Cu2O.Hfus = 17.47E3;
    Cu.Hfus = 13E3;
    Al2O3.Hfus = 111E3;
    Al.Hfus = 10.79E3;

    system->compo[0] = CuO;
    system->compo[1]= Cu2O;
    system->compo[2]= Cu;
    system->compo[3]= Al2O3;
    system->compo[4] = Al;
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
    
}
