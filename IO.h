//
//  IO.h
//  CPC
//
//  Created by Sarah Brotman on 05/01/2018.
//  Copyright © 2018 Sarah Brotman. All rights reserved.
//

#ifndef IO_h
#define IO_h

#include <stdio.h>
#include "System.h"

#define Na 6.02214129E23
#define R 8.3144621
#define PI 3.14159265358979323846

void initSystem(System *system);

void initInteraction(System *system);

void initResults(System *system);

void writeResults(System *system);

#endif /* IO_h */
