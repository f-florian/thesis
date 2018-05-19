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

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex_math.h>

#include "eigen.h"
#include "interpolation.h"
#include "differential.h"

#include <cstdlib>
#include <vector>

#define HAVE_INLINE
// #define GSL_RANGE_CHECK_OFF


using namespace interpolation;
using namespace std;

namespace eigen {
    namespace {
        gsl_matrix *A,*F;
        double * pi=NULL;
    }
    void freemem()
    {
        if(A!=NULL)
            gsl_matrix_free(A);
        if(F!=NULL)
            gsl_matrix_free(F);
        if(pi!=NULL)
            free(pi);
    }
    void init(const size_t size, const Differential::Type type, const double points[], const size_t npts)
    {
        Differential d(size+1, type);
        size_t size_m=(npts-1)*size;
        freemem();
        // allocation
        A=gsl_matrix_calloc(size_m,size_m);
        F=gsl_matrix_alloc(size_m,size_m);
        // pi=(double *)malloc(size_m*sizeof(double));
        // double nodenew,nodeold=0;
        // size_t maxidx;
        // for (size_t k=1; k < npts; ++k)
        //     for(size_t i=0; i < size; i++){
        //         if(k==1 && i==0)
        //             ++i;
        //         nodenew=d.nodes(i,points[k-1],points[k]);
        //         //compute integral (nodeold,nodenew)
        //         if(pi[size*(k-1)+i]/pi[1]!=0) //check 0/1?
        //             maxidx=size*(k-1)+i;
        //         nodeold=nodenew;
        //     }
        // // check meaningful idx: 0/1?
        // if(maxidx>5)
        //     maxidx-=5;
        // else
        //     maxidx=0;
        // for (size_t i=0; i<= maxidx; ++i)
        //     pi[i]/=pi[maxidx];

        for (size_t k=1; k < npts; ++k)
            for(size_t i=0; i < size; i++){
                auto nodei=d.nodes(i+1,points[k-1],points[k]);
                auto dasi=S0(nodei);
                for (size_t l=1; l < npts; ++l)
                    for(size_t j=0; j<=size; j++){
                        if(j==0) {
                            if(l==1)
                                continue;
                            else
                                gsl_matrix_set(F,size*(k-1)+i,size*(l-1)-1,dasi*interpolation::beta(nodei,d.nodes(0,points[l-1],points[l]))*d.quadratureWeights(0,points[l-1],points[l]));
                        } else {
                            gsl_matrix_set(F,size*(k-1)+i,size*(l-1)+j-1,dasi*interpolation::beta(nodei,d.nodes(j,points[l-1],points[l]))*d.quadratureWeights(j,points[l-1],points[l]));
                        }
                    }
                for(size_t j=0; j<=size; j++){
                    if(k==1 && j==0)
                        ++j;
                    gsl_matrix_set(A,size*(k-1)+i,size*(k-1)+j-1,d.differentiationWeights(j,i,points[k-1],points[k]));
                    // gsl_matrix_set(A,size*(k-1)+i-1,size*(k-1)+j-1,pi[size*(l-1)+i-1]*d.differentiationWeights(j,i,points[k-1],points[k]));
                }
                (*gsl_matrix_ptr(A,size*(k-1)+i,size*(k-1)+i))+=interpolation::gamma(nodei)+interpolation::mu(nodei);
                // (*gsl_matrix_ptr(A,size*(k-1)+i-1,size*(k-1)+i-1))+=interpolation::gamma(nodei)*pi[size*(l-1)+i-1];
            }
        // for(int i=0;i<size;i++){
        //     for (int j=0; j<size; ++j)
        //         printf("%.10e ", gsl_matrix_get(A,i,j));
        //     printf(";\n");
        // }
    }
    double computeSpectralRadius(double hint)
    {
        double a, r0=0;
        gsl_vector_complex * alpha=gsl_vector_complex_alloc(A->size1);
        // general eigenvalue problem
        gsl_eigen_gen_workspace *ws=gsl_eigen_gen_alloc(A->size1);
        gsl_vector * beta=gsl_vector_alloc(A->size1);
        gsl_eigen_gen(F,A,alpha,beta,ws);
        for(int i=0; i<A->size1; i++){
            a=gsl_complex_abs(gsl_vector_complex_get(alpha,i))/abs(gsl_vector_get(beta,i));
            if (hint) {
                if(abs(a-hint)<abs(r0-hint))
                    r0=a;
            } else if (r0<a) {
                r0=a;
            }
        }
        gsl_vector_free(beta);
        gsl_vector_complex_free(alpha);
        gsl_eigen_gen_free(ws);
        // return make_pair(r1,r2);
        return r0;
    }
}
