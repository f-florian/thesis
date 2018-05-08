#define HAVE_INLINE
#define GSL_RANGE_CHECK_OFF

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>

#include "eigen.h"
#include "functions.h"
#include <differential.h>

#include <cstdlib>
#include <cstring>
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
        Differential d(size, type);
        size_t size_m=(npts-1)*size-2;
        freemem();
        // allocation
        gsl_matrix *H=gsl_matrix_calloc(size_m+2,size_m+2);
        M=gsl_matrix_alloc(size_m,size_m);
        B=gsl_matrix_calloc(size_m,size_m);
        gsl_matrix *tmp=gsl_matrix_alloc(size_m+2,size_m+2);
        for (size_t k=1; k < npts; ++k)
            for(size_t i=0; i < size; i++)
                for(size_t j=0; j<size; j++)
                    gsl_matrix_set(H,size*(k-1)+i,size*(k-1)+j,d.differentiationWeights(j,i,points[k-1],points[k]));
        gsl_matrix_memcpy(tmp,H);
        for (size_t k=1; k < npts; ++k)
            for(size_t i=0; i < size; i++) {
                auto nodei=d.nodes(i,points[k-1],points[k]);
                auto scale=parameters::D(nodei);
                for(size_t j=0; j<size; j++)
                    (*gsl_matrix_ptr(tmp,size*(k-1)+i,size*(k-1)+j))*=scale;
                (*gsl_matrix_ptr(tmp,size*(k-1)+i,size*(k-1)+i))+=parameters::c(nodei);
            }
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, H, tmp, 0, tmp);

        //copy; TODO: memcpy rows segments
        for(size_t i=0; i<M->size1; ++i)
            for(size_t j=0; j<M->size2; ++j)
                gsl_matrix_set(M,i,j,gsl_matrix_get(tmp,i+1,j+1));

        // TODO: inverse
        gsl_matrix *eqmatrix=gsl_matrix_alloc(2,2);
        gsl_matrix_set(eqmatrix,0,0,c(0)-D(0)*gsl_matrix_get(H,0,0));
        gsl_matrix_set(eqmatrix,0,1,-D(0)*gsl_matrix_get(H,0,size_m+1));
        gsl_matrix_set(eqmatrix,1,0,-D(points[npts-1])*gsl_matrix_get(H,size_m+1,0));
        gsl_matrix_set(eqmatrix,1,1,c(points[npts-1])-D(points[npts-1])*gsl_matrix_get(H,size_m+1,size_m+1));

        gsl_matrix *nonhom=gsl_matrix_alloc(2,size_m);
        for (size_t i = 0; i < size_m; ++i) {
            gsl_matrix_set(nonhom,0,i,D(0)*gsl_matrix_get(H,0,i+1));
            gsl_matrix_set(nonhom,size_m-1,i,D(points[npts-1])*gsl_matrix_get(H,size_m+1,i+1));
        }
        gsl_blas_dgemm(CblasNoTrans,  CblasNoTrans, 1,
                       eqmatrix, nonhom, 0, nonhom);
        gsl_matrix *vcol = gsl_matrix_alloc(size_m,1);
        gsl_matrix *vrow = gsl_matrix_alloc(1,size_m);

        for (size_t i=0; i < size_m; ++i) {
            gsl_matrix_set(vcol,i,0,gsl_matrix_get(tmp,i+1,0));
            gsl_matrix_set(vrow,i,0,gsl_matrix_get(nonhom,0,i));
        }
        gsl_blas_dgemm(CblasNoTrans,  CblasNoTrans, 1, vcol, vrow, 1, M);
        for (size_t i=0; i < size_m; ++i) {
            gsl_matrix_set(vcol,i,0,gsl_matrix_get(tmp,i+1,size_m+1));
            gsl_matrix_set(vrow,i,0,gsl_matrix_get(nonhom,1,i));
        }
        gsl_blas_dgemm(CblasNoTrans,  CblasNoTrans, 1, vcol, vrow, 1, M);
        
        for (size_t k=1; k < npts; ++k)
            for(size_t i=0; i < size; i++)
                if(k==1 && i==0) {
                    continue;
                } else if (k==npts-1 && i==size-1) {
                    break;
                } else {
                    auto nodei=d.nodes(i,points[k-1],points[k]);
                    (*gsl_matrix_ptr(M,size*(k-1)+i,size*(k-1)+i))+=parameters::beta(nodei)+parameters::mu(nodei);
                    gsl_matrix_set(B,size*(k-1)+i,size*(k-1)+i,2*parameters::beta(nodei));
                }               

        // for(int i=0;i<size;i++){
        //     for (int j=0; j<size; ++j)
        //         printf("%.10e ", gsl_matrix_get(M,i,j));
        //     printf(";\n");
        // }
        // gsl_matrix_free(H);
        // gsl_matrix_free(tmp);
        // gsl_matrix_free(vcol);
        // gsl_matrix_free(vrow);
        // gsl_matrix_free(nonhom);
        // gsl_matrix_free(eqmatrix);
    }
    double computeSpectralRadius(const double hint)
    {
        double a, r0=0;
        gsl_vector_complex * alpha=gsl_vector_complex_alloc(M->size1);
        // general eigenvalue problem
        gsl_eigen_gen_workspace *ws=gsl_eigen_gen_alloc(M->size1);
        gsl_vector * beta=gsl_vector_alloc(M->size1);
        gsl_eigen_gen(B,M,alpha,beta,ws);
        for(size_t i=0; i<M->size1; i++){
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
        // gsl_eigen_gen_free(ws);
        return r0;
    }
}
