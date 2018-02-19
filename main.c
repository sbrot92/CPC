//
//  main.c
//  CPC
//
//  Created by Sarah Brotman on 05/01/2018.
//  Copyright © 2018 Sarah Brotman. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include "System.h"
#include "IO.h"


int main() {
    
    System iniSystem;
    System results;

    init(&iniSystem);
    init(&results);
    
    results.t = 0;
    results.T = results.To;

    results.lastWrite = 0;
    int lastPercent = 0;

    while (results.t < results.t_final) {
    	if (results.t > results.heatTime && results.heatFlag == 0) {
    		printf("HEATING COMPLETE at %e",results.t);
    		results.heatFlag = 1;
    	}
    	double percent = (results.t/results.t_final)*100;
    	//printf("percent = %lf\n", percent);
    	int intPercent = (int)floor(percent);
    	//printf("intpercent = %d\n", intPercent);
    	if (intPercent != lastPercent && intPercent % 5 == 0) {
    		printf("Simulation completion : %d%%\n",intPercent);
    		lastPercent = intPercent;
    	}

        /*********** DEBUGGING ***********
    	 if ((results.t > 3.800988e-04) || (results.t - results.lastWrite > results.timeToWrite)) {
                printf("fluxOUT HERE A = %e\n",results.flux);
                printf("Nmol CuO = %e\n",results.compo[0].Nmol);
             printf("Nmol Cu2O = %e\n",results.compo[1].Nmol);
             printf("Nmol Cu = %e\n",results.compo[2].Nmol);
             printf("Nmol Al2O3 = %e\n",results.compo[3].Nmol);
             printf("Nmol Al = %e\n",results.compo[4].Nmol);
             printf("Nmol O2 = %e\n",results.compo[5].Nmol);
             printf("VreleaseCu %e\n",results.VreleasedCu);
             printf("VreleaseAl %e\n",results.V_Al2O3);
    	    		 }
         ***********/
        
    	//if (results.T > 1300) {
    		calcConcentration(&results);
    		 //if ((results.t > 3.800988e-04) || (results.t - results.lastWrite > results.timeToWrite)) {
    		   //      printf("fluxOUT HERE B = %e\n",results.flux);
    		 //}
    		calcMassTrans((&results));
    	//}
            //if ((results.t > 3.800988e-04) || (results.t - results.lastWrite > results.timeToWrite)) {
            //	printf("fluxOUT HERE = %e\n",results.flux);
            //}
        calcHeatTrans((&results));

        if (results.t - results.lastWrite > results.timeToWrite) {
            writeResults(&results);
            results.lastWrite = results.t;
        }
        results.t += results.dt;
    }
    
    printf("Simulation completion : 100%%\n");
    //free(results.compo);
    return 0;
}
