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

#define HAVE_INLINE
// #define GSL_RANGE_CHECK_OFF

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>

#include "eigen.h"
#include "functions.h"
#include <differential.h>

#include <cstdlib>
#include <cstring>
#include <limits>
#include <vector>


using namespace std;
using namespace parameters;

namespace eigen {
    namespace {
        gsl_matrix *B,*M;
    }
    void freemem()
    {
        if(B!=NULL)
            gsl_matrix_free(B);
        if(M!=NULL)
            gsl_matrix_free(M);
    }
    void init(const size_t size, const Differential::Type type, const double points[], const size_t npts)
    {
        freemem();
        Differential d(size, type);
        size_t size_m=(npts-1)*size-2;
        gsl_matrix *H=gsl_matrix_calloc(size_m+2,size_m+2);
        gsl_matrix *tmp=gsl_matrix_alloc(size_m+2,size_m+2);
        gsl_matrix *tmp0=gsl_matrix_alloc(size_m+2,size_m+2);

        if (D(0))
            size_m=(npts-1)*size-2;
        else
            size_m=(npts-1)*size-1;

        // allocation
        M=gsl_matrix_alloc(size_m,size_m);
        B=gsl_matrix_calloc(size_m,size_m);

        // construct differentiotion matrix
        for (size_t k=1; k < npts; ++k)
            for(size_t i=0; i < size; i++)
                for(size_t j=0; j<size; j++)
                    gsl_matrix_set(H,size*(k-1)+i,size*(k-1)+j,d.differentiationWeights(j,i,points[k-1],points[k]));

        if (D(0)) {
            gsl_matrix_memcpy(tmp,H);

            for (size_t k=1; k < npts; ++k)
                for(size_t i=0; i < size; i++) {
                    auto nodei=d.nodes(i,points[k-1],points[k]);
                    auto scale=-parameters::D(nodei);
                    for(size_t j=0; j<size; j++)
                        (*gsl_matrix_ptr(tmp,size*(k-1)+i,size*(k-1)+j))*=scale;
                    (*gsl_matrix_ptr(tmp,size*(k-1)+i,size*(k-1)+i))+=parameters::c(nodei);
                }
            
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, H, tmp, 0, tmp0);
            
            for (size_t k=0; k < npts-1; ++k)
                for(size_t i=0; i < size; i++) {
                    auto nodei=d.nodes(i,points[k],points[k+1]);
                    (*gsl_matrix_ptr(tmp0,size*k+i,size*k+i))+=parameters::beta(nodei)+parameters::mu(nodei);
                    if((k!=0 || i!=0) && (k!=npts-2 || i!= size-1))
                        gsl_matrix_set(B,size*k+i-1,size*k+i-1,2*parameters::beta(nodei));
            }
            // allocation       
            gsl_matrix *eqmatrix=gsl_matrix_alloc(2,2);
            gsl_matrix *nonhom=gsl_matrix_alloc(2,size_m);
            gsl_matrix *nonhom0=gsl_matrix_alloc(2,size_m);

            //copy; TODO: memcpy rows segments
            for(size_t i=0; i<M->size1; ++i)
                for(size_t j=0; j<M->size2; ++j)
                    gsl_matrix_set(M,i,j,gsl_matrix_get(tmp0,i+1,j+1));

            gsl_matrix_set(eqmatrix,1,1,c(points[0])-D(points[0])*gsl_matrix_get(H,0,0));
            gsl_matrix_set(eqmatrix,0,1,D(points[0])*gsl_matrix_get(H,0,size_m+1));
            gsl_matrix_set(eqmatrix,1,0,D(points[npts-1])*gsl_matrix_get(H,size_m+1,0));
            gsl_matrix_set(eqmatrix,0,0,c(points[npts-1])-D(points[npts-1])*gsl_matrix_get(H,size_m+1,size_m+1));
            double scale=gsl_matrix_get(eqmatrix,0,0)*gsl_matrix_get(eqmatrix,1,1)-gsl_matrix_get(eqmatrix,1,0)*gsl_matrix_get(eqmatrix,0,1);
            gsl_matrix_scale(eqmatrix, 1./scale);

            for (size_t i = 0; i < size_m; ++i) {
                gsl_matrix_set(nonhom,0,i,D(points[0])*gsl_matrix_get(H,0,i+1));
                gsl_matrix_set(nonhom,1,i,D(points[npts-1])*gsl_matrix_get(H,size_m+1,i+1));
            }
            gsl_blas_dgemm(CblasNoTrans,  CblasNoTrans, 1, eqmatrix, nonhom, 0, nonhom0);

            for (size_t i=0; i < size_m; ++i) {
                gsl_matrix_set(nonhom,0,i,gsl_matrix_get(tmp0,i+1,0));
                gsl_matrix_set(nonhom,1,i,gsl_matrix_get(tmp0,i+1,size_m+1));
            }
            gsl_blas_dgemm(CblasTrans,  CblasNoTrans, 1, nonhom, nonhom0, 1, M);

            // for(int i=0;i<size_m;i++){
            //     for (int j=0; j<size_m; ++j)
            //         printf("%.10e ", gsl_matrix_get(M,i,j));
            //     printf(";\n");
            // }
            // printf("\n");
            gsl_matrix_free(nonhom);
            gsl_matrix_free(nonhom0);
            gsl_matrix_free(eqmatrix);
        } else { //D=0
            //copy; TODO: memcpy rows segments
            for(size_t i=0; i<M->size1; ++i) {
                auto nodei=d.nodes(i+1,points[0],points[1]);
                for(size_t j=0; j<M->size2; ++j)
                    gsl_matrix_set(M,j,i,gsl_matrix_get(H,j+1,i+1)*c(nodei));
                gsl_matrix_set(B,i,i,2*beta(nodei));
                *(gsl_matrix_ptr(M,i,i))+=beta(nodei)+mu(nodei);
            }

            // for(int i=0;i<size_m;i++){
            //     for (int j=0; j<size_m; ++j)
            //         printf("%.10e ", gsl_matrix_get(M,i,j));
            //     printf(";\n");
            // }
            // printf("\n");
        }
        gsl_matrix_free(H);
        gsl_matrix_free(tmp);
        gsl_matrix_free(tmp0);
    }
    pair<double,gsl_vector_complex *> computeSpectralRadius(const double hint)
    {
        double a, r0=0;
        gsl_vector_complex * alpha=gsl_vector_complex_alloc(M->size1);
        // general eigenvalue problem
        gsl_vector * beta=gsl_vector_alloc(M->size1);
        gsl_matrix_complex *vect=gsl_matrix_complex_alloc(M->size1,M->size1);

        // gsl_linalg_QR_decomp(M,beta);
        // gsl_blas_dtrsm(CblasRight,CblasUpper,CblasNoTrans,CblasNonUnit,1,M,B);
        // gsl_eigen_nonsymm_workspace *ws=gsl_eigen_nonsymm_alloc(M->size1);
        // gsl_eigen_nonsymm(B,alpha,ws);
        
        gsl_eigen_genv_workspace *ws=gsl_eigen_genv_alloc(M->size1);
        gsl_eigen_genv(B,M,alpha,beta,vect,ws);
        size_t j=0;
        for(size_t i=0; i<M->size1; i++){
            a=gsl_complex_abs(gsl_vector_complex_get(alpha,i))/abs(gsl_vector_get(beta,i));
            // a=gsl_complex_abs(gsl_vector_complex_get(alpha,i));
            if (hint) {
                if(abs(a-hint)<abs(r0-hint))
                    r0=a;
            } else if (r0<a) {
                r0=a;
                j=i;
            }
            // printf("%e ", a);
        }
        // printf("\n");

        // printf("eigv:\n");
        // for (size_t i = 0; i < M->size1; ++i) {
        //         for (size_t k = 0; k < M->size1; ++k)
        //         printf("%e %e\t", gsl_matrix_complex_get(vect,i, k).dat[0], gsl_matrix_complex_get(vect,i,k).dat[1]);
        //     printf("\n");
        // }

        // gsl_eigen_nonsymm_free(ws);
        gsl_eigen_genv_free(ws);
        gsl_vector_free(beta);
        gsl_matrix_complex_get_col(alpha,vect,j);
        gsl_matrix_complex_free(vect);
        return {r0,alpha};
        gsl_vector_complex_free(alpha);
    }
    double solveExplicit()
    {
        constexpr double epsabs=1e-10;
        constexpr double epsrel=1e-10;
        constexpr size_t maxint=1000;
        constexpr double lowlimit=-1e10;
        gsl_integration_cquad_workspace *wse=gsl_integration_cquad_workspace_alloc(300);
        gsl_integration_workspace * wsi=gsl_integration_workspace_alloc(maxint);
        double lambda;
        double result, aerr;
        size_t neval;

        gsl_function fe;
        fe.function=[](double x, void * lambda) {
            double a=2*betaonc(x)/(*static_cast<double*>(lambda))-(beta(x)+mu(x))/c(x);
            if(gsl_isinf(a)==1){
                fprintf(stderr,"using %e for inf\n", numeric_limits<double>::max()/20);
                return numeric_limits<double>::max()/20;
            } else if (a<lowlimit||isinf(a)==-1) {
                // fprintf(stderr,"using %e for -inf\n", lowlimit);
                return lowlimit;
            } else if(!gsl_finite(a)){
                fprintf(stderr,"singular value at %e; inf: %d nan: %d\n", x, gsl_isinf(a), gsl_isnan(a));
            }
            return a;
        };
        fe.params=&lambda;

        auto pr=make_pair(&fe, wsi);
        gsl_set_error_handler_off();

        gsl_function lower;
        lower.function=[](double x, void * param) {
            // printf("integrating internal to %e\n", x);
            double result, aerr;
            gsl_integration_qag(static_cast<pair<gsl_function*,gsl_integration_workspace*>*>(param)->first, 0, x, epsabs, epsrel, maxint, GSL_INTEG_GAUSS61, static_cast<pair<gsl_function*,gsl_integration_workspace*>*>(param)->second, &result, &aerr);
            return exp(result)/c(x);
        };
        lower.params=&pr;

        for (lambda=.01; lambda>0; lambda-=.001) {
            gsl_integration_cquad(&lower, 0, 1, epsabs, epsrel, wse, &result, &aerr, &neval);
            printf("%f %e %e %ld\n", lambda, result, aerr, neval);
        }

        gsl_integration_cquad_workspace_free(wse);
        gsl_integration_workspace_free(wsi);
        return 0;
    }
}
