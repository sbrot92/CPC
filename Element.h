//
//  Element.h
//  CPC
//
//  Created by Sarah Brotman on 19/02/2018.
//  Copyright © 2018 Sarah Brotman. All rights reserved.
//

#ifndef Element_h
#define Element_h

#include <stdio.h>


typedef struct System {
    // system state
    double t, t_final, To, dT, dTdt, T, dt, C, timeToWrite;
    double heatTime, heatPower;
    
    // initial variables
    double Ro_CuO, Ro_Al, Ro_Al2O3, TMD, stoich;
    double V_chamber, V_thermite, Nmol_thermite;
    double n_Al, n_Al2O3, n_CuO;
    
    // diffusion
    double Cs_ratio,Cs_Cu, Cs_Al2O3, C0, C1, C2, v, flux, dflux;
    double VreleasedCu,VreleasedAl,V_Al2O3;
    double D_Cu2O,D_Al2O3;
    
    // mass transfer
    double w_Al, w_Al2O3, w_Cu, w_CuO, w_Cu2O;
    double dw_Al, dw_Al2O3, dw_Cu, dw_CuO, dw_Cu2O;
    double dCu2O_Nmol, dCu_Nmol, dO2_Nmol, dCuO_Nmol, dAl_Nmol, dAl2O3_Nmol, Nmol_O2used;
    double a0, a1, a2, leftBound, rightBound;
    double da0dt, da1dt, da2dt;
    
    // geometric considerations
    double R1, R2, r, d, d1, d2, theta, alpha,a,S;
    
    double lastWrite;
    
    double currHfus;
    
    // Species information
    Species * compo;
    
} System;

#endif /* Element_h */
