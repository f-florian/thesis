#define HAVE_INLINE
// #define GSL_RANGE_CHECK_OFF

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
        M=gsl_matrix_alloc(size_m,size_m);
        B=gsl_matrix_calloc(size_m,size_m);
        gsl_matrix *H=gsl_matrix_calloc(size_m+2,size_m+2);
        gsl_matrix *tmp=gsl_matrix_alloc(size_m+2,size_m+2);
        gsl_matrix *tmp0=gsl_matrix_alloc(size_m+2,size_m+2);
        gsl_matrix *eqmatrix=gsl_matrix_alloc(2,2);
        gsl_matrix *nonhom=gsl_matrix_alloc(2,size_m);
        gsl_matrix *nonhom0=gsl_matrix_alloc(2,size_m);
        gsl_matrix *vcol = gsl_matrix_alloc(size_m,1);
        gsl_matrix *vrow = gsl_matrix_alloc(1,size_m);
        for (size_t k=1; k < npts; ++k)
            for(size_t i=0; i < size; i++)
                for(size_t j=0; j<size; j++)
                    gsl_matrix_set(H,size*(k-1)+i,size*(k-1)+j,d.differentiationWeights(j,i,points[k-1],points[k]));

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

        for (size_t k=1; k < npts; ++k)
            for(size_t i=0; i < size; i++) {
                auto nodei=d.nodes(i,points[k-1],points[k]);
                (*gsl_matrix_ptr(tmp0,size*(k-1)+i,size*(k-1)+i))+=parameters::beta(nodei)+parameters::mu(nodei);
                if((k!=1 || i!=0) && (k!=npts-1 || i!= size-1))
                    gsl_matrix_set(B,size*(k-1)+i-1,size*(k-1)+i-1,2*parameters::beta(nodei));
            }


        //copy; TODO: memcpy rows segments
        for(size_t i=0; i<M->size1; ++i)
            for(size_t j=0; j<M->size2; ++j)
                gsl_matrix_set(M,i,j,gsl_matrix_get(tmp0,i+1,j+1));

        gsl_matrix_set(eqmatrix,1,1,c(0)-D(0)*gsl_matrix_get(H,0,0));
        gsl_matrix_set(eqmatrix,0,1,D(0)*gsl_matrix_get(H,0,size_m+1));
        gsl_matrix_set(eqmatrix,1,0,D(points[npts-1])*gsl_matrix_get(H,size_m+1,0));
        gsl_matrix_set(eqmatrix,0,0,c(points[npts-1])-D(points[npts-1])*gsl_matrix_get(H,size_m+1,size_m+1));
        double scale=gsl_matrix_get(eqmatrix,0,0)*gsl_matrix_get(eqmatrix,1,1)-gsl_matrix_get(eqmatrix,1,0)*gsl_matrix_get(eqmatrix,0,1);
        gsl_matrix_scale(eqmatrix, 1./scale);

        for (size_t i = 0; i < size_m; ++i) {
            gsl_matrix_set(nonhom,0,i,D(0)*gsl_matrix_get(H,0,i+1));
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
        gsl_matrix_free(H);
        gsl_matrix_free(tmp);
        gsl_matrix_free(tmp0);
        gsl_matrix_free(vcol);
        gsl_matrix_free(vrow);
        gsl_matrix_free(nonhom);
        gsl_matrix_free(nonhom0);
        gsl_matrix_free(eqmatrix);
    }
    pair<double,gsl_vector_complex *> computeSpectralRadius(const double hint)
    {
        double a, r0=0;
        gsl_vector_complex * alpha=gsl_vector_complex_alloc(M->size1);
        // general eigenvalue problem
        gsl_vector * beta=gsl_vector_alloc(M->size1);
        gsl_matrix_complex *vect=gsl_matrix_complex_alloc(M->size1,M->size1);
        gsl_eigen_genv_workspace *ws=gsl_eigen_genv_alloc(M->size1);
        gsl_eigen_genv(B,M,alpha,beta,vect,ws);
        size_t j=0;
        for(size_t i=0; i<M->size1; i++){
            a=gsl_complex_abs(gsl_vector_complex_get(alpha,i))/abs(gsl_vector_get(beta,i));
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

        gsl_eigen_genv_free(ws);
        gsl_vector_free(beta);
        gsl_matrix_complex_get_row(alpha,vect,j);
        gsl_vector_complex_free(alpha);
        gsl_matrix_complex_free(vect);
        return {r0,nullptr};
    }
}
