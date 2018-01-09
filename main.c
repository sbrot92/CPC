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
#include <sys/timeb.h>
#include "System.h"
#include "IO.h"


int main() {
    
    System iniSystem;
    System results;

    init(&iniSystem);
    init(&results);
    
    while (results.t < results.t_final) {
        calcConcentration(&results);
        calcMassTrans((&results));
        calcHeatTrans((&results));
        results.t += results.dt;
    }
    
    
    
    return 0;
}
