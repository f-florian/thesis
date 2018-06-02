/***************************************************************************
 * Copyright (C) 2018 Francesco Florian
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ***************************************************************************/

#include <cstdlib>
#include <cstdio>
#include <gsl/gsl_vector.h>
#include "eigen.h"
#include "differential.h"
#include "functions.h"

Differential::Type inttype=Differential::Type::ClenshawCurtis;
size_t maxsize=100;
size_t step=1;
size_t start=2;
size_t order=0;
double length=1;
double c0=1, D0=1, beta0=1, mu0=1, c1=0, D1=0;
int p=1;
int main(int argc, char**argv)
{
    fprintf(stderr,"usage: order start step maxsize type(g/cc) length D_1 c_1 mu0 beta0 D0 c0 p  \n");
    switch (argc-1){
    default:
        fprintf(stderr,"Too many arguments (%d supplied); ignoring the last ones\n",argc);
        [[fallthrough]];
    case 13:
        p=atoi(argv[13]);
        fprintf(stderr, "p: %d ", p);
        [[fallthrough]];
    case 12:
        c0=atof(argv[12]);
        fprintf(stderr, "c_0: %e ", c0);
        [[fallthrough]];
    case 11:
        D0=atof(argv[11]);
        fprintf(stderr, "D_0: %e ", D0);
        [[fallthrough]];
    case 10:
        beta0=atof(argv[10]);
        fprintf(stderr, "beta: %e ", beta0);
        [[fallthrough]];
    case 9:
        mu0=atof(argv[9]);
        fprintf(stderr, "mu: %e ", mu0);
        [[fallthrough]];
    case 8:
        c1=atof(argv[8]);
        fprintf(stderr, "c_1: %e ", c1);
        [[fallthrough]];
    case 7:
        D1=atof(argv[7]);
        fprintf(stderr, "D_1: %e ", D1);
        [[fallthrough]];
    case 6:
        length=atof(argv[6]);
        fprintf(stderr, "length: %e ", length);
        [[fallthrough]];
    case 5:
        switch(atoi(argv[5])){
        case 0:
            inttype=Differential::Type::Gauss;
            fprintf(stderr, "Integration type: Gauss; ");
            break;
        case 1:
            inttype=Differential::Type::ClenshawCurtis;
            fprintf(stderr, "Integration type: Clenshaw-Curtis; ");
            break;
        default:
            fprintf(stderr,"Invalid integration type specified (%d), using default 'Clenshaw-Curtis (1)'\n", atoi(argv[6]));
        }
        [[fallthrough]];
    case 4:
        maxsize=atoi(argv[4]);
        [[fallthrough]];
    case 3:
        step=atoi(argv[3]);
        if(step==0) {
            step=1;
            fprintf(stderr,"Invlid step 0, setting 1\n");
        }
        [[fallthrough]];
    case 2:
        start=atoi(argv[2]);
        if (start>maxsize)
            fprintf(stderr,"minimum size (%ld) > maximum size (%ld); nothing to compute", start, maxsize);
        [[fallthrough]];
    case 1:
        order=atoi(argv[1]);
        [[fallthrough]];
    case 0:
        ;
    }
    fprintf(stderr, "\n");
        
    //init

    // eigen::solveExplicit();
    // return 0;
    
    //alloc mesh
    double* mesh;
    if (order)
        mesh=(double*) malloc((maxsize+1)*sizeof(double));
    else
        mesh=NULL;


    parameters::init(length, c0, D0, beta0, mu0, c1, D1, p);
    for(size_t i=start; i<=maxsize;i+=step){
        eigen::init<2>(i+1, inttype, {0, length});
        auto a=eigen::computeSpectralRadius(0);
        printf("%4ld %.20e\n", i, a.first);
    }

    // r_0 vs \bar D
    // for(double dbar=0; dbar<1; dbar+=1e-3){
    //     parameters::init(length, c0, dbar*D0, beta0, mu0, c1, dbar*D1, p);
    //     eigen::init<2>(400, inttype, {0, length});
    //     auto a=eigen::computeSpectralRadius(0);
    //     printf("%.20e %.20e\n", dbar, a.first);
    // }
    free(mesh);
    eigen::freemem();
}
