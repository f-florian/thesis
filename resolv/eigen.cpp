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
#include <cmath>

#define HAVE_INLINE
// #define GSL_RANGE_CHECK_OFF

using namespace interpolation;
using namespace std;

namespace eigen {
    namespace {
        gsl_matrix *A,*F;
    }
    void freemem()
    {
        if(A!=NULL)
            gsl_matrix_free(A);
        if(F!=NULL)
            gsl_matrix_free(F);
    }
    void init(const size_t size, const Differential::Type type, const double points[], const size_t npts)
    {
        Differential d(size+1, type);
        size_t size_m=(npts-1)*size;
        freemem();
        // allocation
        A=gsl_matrix_calloc(size_m,size_m);
        F=gsl_matrix_alloc(size_m,size_m);

        for (size_t k=0; k < npts-1; ++k)
            for(size_t i=0; i < size; i++) {
                auto nodei=d.nodes(i+1,points[k],points[k+1]);
                auto dasi=S0(nodei);
                for (size_t l=0; l < npts-1; ++l)
                    for(size_t j=0; j<=size; j++) {
                        if(j==0) {
                            if(l==0)
                                continue;
                            else
                                (*gsl_matrix_ptr(F,size*k+i,size*l-1))+=dasi*interpolation::beta(nodei,d.nodes(0,points[l],points[l+1]))*d.quadratureWeights(0,points[l],points[l+1]);
                        } else {
                            gsl_matrix_set(F,size*k+i,size*l+j-1,dasi*interpolation::beta(nodei,d.nodes(j,points[l],points[l+1]))*d.quadratureWeights(j,points[l],points[l+1]));
                        }
                    }
                for(size_t j=0; j<=size; j++){
                    if(k==0 && j==0)
                        ++j;
                    gsl_matrix_set(A,size*k+i,size*k+j-1,d.differentiationWeights(j,i+1,points[k],points[k+1]));
                }
                (*gsl_matrix_ptr(A,size*k+i,size*k+i))+=interpolation::gamma(nodei)+interpolation::mu(nodei);
            }

        // double acc[1000];
        // for (size_t i = 0; i < size; ++i) {
        //     acc[i]=4*pow(d.nodes(i+1,points[0],points[1]),3);
        //     for (size_t j = 0; j < size; ++j) {
        //         // printf("%.10e ", d.differentiationWeights(j+1,i+1,points[0],points[1]));
        //         acc[i]-=pow(d.nodes(j+1,points[0],points[1]),4)*gsl_matrix_get(A,i,j);
        //     }
        //     // printf(";\n");        
        // }
        // printf("\n\n\n");
        // for (size_t i = 0; i < size; ++i)
        //     fprintf(stderr, "%.5e\t", acc[i]);

        // for (size_t i = 0; i < size; ++i) {
        //     for (size_t j = 0; j < size; ++j) {
        //         printf("%.10e ", gsl_matrix_get(A,i,j));
        //     }
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
