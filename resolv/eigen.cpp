#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex_math.h>

#include "eigen.h"
#include "interpolation.h"
#include "differential.h"

#include <cstdlib>
#include <vector>

#undef STDEIG

#ifdef STDEIG
#include <gsl/gsl_blas.h>
#endif

#define HAVE_INLINE
#define GSL_RANGE_CHECK_OFF


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
        Differential d(size, type);
        size_t size_m=(npts-1)*size;
        freemem();
        // allocation
        A=gsl_matrix_alloc(size_m,size_m);
        F=gsl_matrix_alloc(size_m,size_m);
        for (size_t k=1; k < npts; ++k)
            for(size_t i=0; i < size; i++){
                auto nodei=d.nodes(i,points[k-1],points[k]);
                auto dasi=S0(nodei);
                for (size_t l=1; l < npts; ++l)
                    for(size_t j=0; j<size; j++)
                        gsl_matrix_set(F,size*(k-1)+i,size*(l-1)+j,dasi*interpolation::beta(nodei,d.nodes(j,points[l-1],points[l]))*d.quadratureWeights(j,points[l-1],points[l]));
                for(size_t j=0; j<size; j++)
                    gsl_matrix_set(A,size*(k-1)+i,size*(k-1)+j,d.differentiationWeights(j,i,points[k-1],points[k]));
                (*gsl_matrix_ptr(A,size*(k-1)+i,size*(k-1)+i))+=interpolation::gamma(nodei)+interpolation::mu(nodei);
            }
        // for(int i=0;i<size;i++){
        //     for (int j=0; j<size; ++j)
        //         printf("%.10e ", gsl_matrix_get(A,i,j));
        //     printf(";\n");
        // }
    }
    double computeSpectralRadius()
    {
        double a,r0=0;
        gsl_vector_complex * alpha=gsl_vector_complex_alloc(A->size1);
#ifndef STDEIG
        // general eigenvalue problem
        gsl_eigen_gen_workspace *ws=gsl_eigen_gen_alloc(A->size1);
        gsl_vector * beta=gsl_vector_alloc(A->size1);
        gsl_eigen_gen(F,A,alpha,beta,ws);
        for(int i=0; i<A->size1; i++){
            a=gsl_complex_abs(gsl_vector_complex_get(alpha,i))/abs(gsl_vector_get(beta,i));
            if(r0<a)
                r0=a;
            // if(r1<a){
            //     r2=r1;
            //     r1=a;
            // } else if (r2<a) {
            //     r2=a;
            // }
        }
        gsl_vector_free(beta);
        gsl_vector_complex_free(alpha);
        gsl_eigen_gen_free(ws);
#else
        // standard eigenvalue problem
        gsl_blas_dtrsm(CblasRight,  CblasLower,  CblasNoTrans, CblasNonUnit, 1, A, F);
        gsl_eigen_nonsymm_workspace *ws2= gsl_eigen_nonsymm_alloc(A->size1);
        gsl_eigen_nonsymm(F, alpha, ws2);
        for(int i=0; i<A->size1; i++){
            a=gsl_complex_abs(gsl_vector_complex_get(alpha,i));
            if(r0<a)
                r0=a;
        }
        gsl_eigen_nonsymm_free(ws2);
#endif
        // return make_pair(r1,r2);
        return r0;
    }
}
