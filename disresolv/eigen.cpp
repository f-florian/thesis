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

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
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
        Differential *d;
    }
    void freemem()
    {
        gsl_matrix_free(A);
        gsl_matrix_free(F);
        delete d;
    }
    void alloc(size_t orderp1, Differential::Type type)
    {
        A=gsl_matrix_calloc(orderp1,orderp1);
        F=gsl_matrix_calloc(orderp1,orderp1);
        d=new Differential(orderp1, type);
    }
    void init(const double agemax)
    {
        auto size=A->size1;
        for(size_t i = 1; i < size; ++i) {
            auto nodei=d->nodes(i,0,agemax);
            auto dasi=S0(nodei);
            for(size_t j = 1; j < size; ++j) {
                gsl_matrix_set(F,i,j,dasi*interpolation::beta(nodei,d->nodes(j,0,agemax))*d->quadratureWeights(j,0,agemax));
                if((j==0)||(j==size-1))
                    continue;
                gsl_matrix_set(A,i,j,d->differentiationWeights(j,i,0,agemax));
            }
            (*gsl_matrix_ptr(A,i,i))+=interpolation::gamma(nodei)+interpolation::mu(nodei);
        }
        gsl_matrix_set(A,0,0,1);
            
        // auto col=gsl_matrix_alloc(size_m,1);
        // for (size_t i = 0; i < size_m; ++i) {
        //     // gsl_matrix_set(A,size*k+i,size*k+j-1,d->differentiationWeights(j,i+1,points[k],points[k+1]));
        //     gsl_matrix_set(col,i,0,d->differentiationWeights(0,i+1,points[0],points[1]));
        // }
        // auto row=gsl_matrix_alloc(1,size_m);
        // for (size_t i = 0; i < size_m; ++i) {
        //     // gsl_matrix_set(A,size*k+i,size*k+j-1,d->differentiationWeights(j,i+1,points[k],points[k+1]));
        //     gsl_matrix_set(row,0,i,d->evalPolynomial(i+1,0,points[0],points[1])/d->evalPolynomial(1,0,points[0],points[1]));
        // }

        // gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,col,row,1,A);

        // double acc[1000];
        // for (size_t i = 0; i < size; ++i) {
        //     acc[i]=4*pow(d->nodes(i+1,points[0],points[1]),3);
        //     for (size_t j = 0; j < size; ++j) {
        //         // printf("%.10e ", d->differentiationWeights(j+1,i+1,points[0],points[1]));
        //         acc[i]-=pow(d->nodes(j+1,points[0],points[1]),4)*gsl_matrix_get(A,i,j);
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
    double computeSpectralRadius()
    {
        double a, r0=0;
        gsl_vector_complex * alpha=gsl_vector_complex_alloc(A->size1);
        // general eigenvalue problem
        gsl_eigen_gen_workspace *ws=gsl_eigen_gen_alloc(A->size1);
        gsl_vector * beta=gsl_vector_alloc(A->size1);
        gsl_eigen_gen(F,A,alpha,beta,ws);
        for(size_t i=0; i<A->size1; i++){
            a=gsl_complex_abs(gsl_vector_complex_get(alpha,i))/abs(gsl_vector_get(beta,i));
            if (r0<a)
                r0=a;
        }
        gsl_vector_free(beta);
        gsl_vector_complex_free(alpha);
        gsl_eigen_gen_free(ws);
        return r0;
    }
}
