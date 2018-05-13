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

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ***************************************************************************/

#include <cstdlib>
#include <cstdio>

#include <gsl/gsl_vector.h>
#include "eigen.h"
#include "interpolation.h"
#include "differential.h"

constexpr int sizes=21;
constexpr int sizemuv[3]={113,116,23};

constexpr double s[sizes]={5187,5310,5582,5995,6254,6512,7310,8317,9802,8744,8059,7556,8429,9839,7713,6369,5063,3196,1408,394,65};
constexpr double mu0[sizemuv[0]]={0.06764,0.00038,0.00024,0.00019,0.00013,0.00010,0.00010,0.00010,0.00009,0.00008,0.00007,0.00007,0.00008,0.00009,0.00012,0.00015,0.00019,0.00024,0.00029,0.00036,0.00042,0.00047,0.00050,0.00053,0.00054,0.00055,0.00055,0.00055,0.00054,0.00055,0.00057,0.00059,0.00061,0.00064,0.00068,0.00072,0.00075,0.00078,0.00082,0.00090,0.00100,0.00110,0.00120,0.00129,0.00140,0.00155,0.00171,0.00190,0.00211,0.00233,0.00255,0.00280,0.00308,0.00339,0.00373,0.00412,0.00454,0.00498,0.00540,0.00585,0.00639,0.00709,0.00795,0.00886,0.00973,0.01070,0.01182,0.01295,0.01415,0.01543,0.01684,0.01846,0.02025,0.02205,0.02388,0.02610,0.02889,0.03233,0.03649,0.04134,0.04680,0.05307,0.06025,0.06839,0.07766,0.08810,0.09941,0.11156,0.12500,0.14002,0.15698,0.17602,0.19751,0.22205,0.24801,0.27055,0.29434,0.31910,0.34485,0.37165,0.39954,0.42855,0.45874,0.49015,0.52284,0.55684,0.59223,0.62905,0.66736,0.70722,0.74869,0.79185,0.83675};
constexpr double mu1[sizemuv[1]]={0.99822,0.99968,0.99980,0.99988,0.99992,0.99992,0.99992,0.99992,0.99993,0.99993,0.99993,0.99993,0.99993,0.99993,0.99992,0.99990,0.99988,0.99987,0.99985,0.99984,0.99983,0.99981,0.99980,0.99978,0.99977,0.99976,0.99975,0.99973,0.99971,0.99970,0.99969,0.99968,0.99966,0.99963,0.99961,0.99959,0.99957,0.99954,0.99950,0.99943,0.99937,0.99931,0.99926,0.99920,0.99913,0.99905,0.99895,0.99884,0.99874,0.99864,0.00148,0.00162,0.00178,0.00193,0.00208,0.00221,0.00233,0.00245,0.00259,0.00279,0.00304,0.00333,0.00362,0.00391,0.00422,0.00460,0.00502,0.00547,0.00597,0.00654,0.00722,0.00800,0.00888,0.00977,0.01077,0.01202,0.01363,0.01559,0.01789,0.02057,0.02361,0.02714,0.03136,0.03630,0.04188,0.04816,0.05538,0.06373,0.07352,0.08454,0.09695,0.11100,0.12757,0.14619,0.16620,0.18796,0.20798,0.22839,0.24917,0.27030,0.29175,0.31351,0.33554,0.35783,0.38033,0.40301,0.42585,0.44879,0.47179,0.49482,0.51783,0.54076,0.56358,0.58622,0.60864,0.63079};
constexpr double mu2[sizemuv[2]]={0.06764,0.00010,0.00007,0.00015,0.00042,0.00055,0.00057,0.00072,0.00100,0.00155,0.00255,0.00412,0.00639,0.01070,0.01684,0.02610,0.04680,0.08810,0.15698,0.27055,0.39954,0.55684,0.74869};

int maxsize=2000;
int sizemu=sizemuv[2];
const double *mu=mu2;

int main(int argc, char**argv)
{
    switch (argc){
    default:
        printf("Ignoring arguments beyond the second\n");
    case 2:
        maxsize=atoi(argv[2]);
        if(maxsize==0)
            printf("You selected 0 as maximum size, no matrix radius will be computed");
    case 1:
        switch(atoi(argv[1])){
        case 0:
            mu=mu0;
            sizemu=sizemuv[0];
            break;
        case 1:
            mu=mu1;
            sizemu=sizemuv[1];
            break;
        case 2:
            break;
        default:
            printf("Ignoring invalid argument selecting mu");
        }
    case 0:
        ;
    }
    
    double xm[sizemu];
    double xs[sizes];

    for(int i=0;i<sizes;i++)
        xs[i]=5*i;
    for(int i=0;i<sizemu;i++)
        xm[i]=5*i;

    interpolation::splineInit(sizes, sizemu, xm, mu, xs, s);

    for(int i=10; i<=maxsize;){
        // size,start,delta,order
        eigen::init(i-1, 100./i, 100./i);
        printf("%4d %.20e\n", i, eigen::computeSpectralRadius());
        if(i<100)
            i+=10;
        else
            i+=100;
    }
    eigen::freemem();
}
