//
//  System.h
//  CPC
//
//  Created by Sarah Brotman on 05/01/2018.
//  Copyright © 2018 Sarah Brotman. All rights reserved.
//

#ifndef System_h
#define System_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <errno.h>

#define Na 6.02214129E23
#define R 8.3144621
#define PI 3.14159265358979323846

typedef struct Species {
    char *name;
    double Do;
    double Ea;
    // mole. radius, molar mass, melt temp, heat of melt
    double a, M, Tfus, Hfus, p, Vm;
    int flagMelt;
    // heat of formation
    double Hf, Cv[40],TCv[40];
    double Cp;
    double molarRatio, Nmol;
} Species;


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
    double V_CuO, V_Cu2O, V_Cu, V_Al;
    double VreleasedCu,VreleasedAl,V_Al2O3,dV_Cudt;
    double D_Cu2O,D_Al2O3;
    
    // mass transfer
    double w_Al, w_Al2O3, w_Cu, w_CuO, w_Cu2O;
    double dw_Al, dw_Al2O3, dw_Cu, dw_CuO, dw_Cu2O;
    // per pair
    double dCu2O_Nmol, dCu_Nmol, dO2_Nmol, dCuO_Nmol, dAl_Nmol, dAl2O3_Nmol, Nmol_O2used;
    double a0, a1, a2, leftBound, rightBound;
    double da0dt, da1dt, da2dt;
    
    // geometric considerations
    double R1, R2, r, d, d1, d2, theta, alpha,a,S;
    
    double lastWrite;

    double currHfus;

    int flagNoMoreAl;
    int flagNoMoreCuO;
    int flagNoMoreCu2O;
    int heatFlag;



    // Species information
    Species * compo;
    
} System;

void init(System *system);

double calcDiffCoeff(double D0, double Ea, double T);

void calcConcentration(System *system);
void calcMassTrans(System *system);
void calcHeatTrans(System *system);
void calcHeatCapacity(System *system);
double calcCp(double A, double B, double C, double D, double E, double t);
double calcInteractionRadius(System * system, double x);
double calc_dV(System * system, double x0, double x1);





#endif /* System_h */

