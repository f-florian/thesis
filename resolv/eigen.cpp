#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex_math.h>
// #include <gsl/gsl_blas.h>

#include "eigen.h"
#include "interpolation.h"
#include "differential.h"

#include <cstdlib>

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
    void init(const unsigned int size, double start, double end, Differential::Type type)
    {
        Differential d(size, type);
        freemem();
        // allocation
        A=gsl_matrix_alloc(size,size);
        F=gsl_matrix_alloc(size,size);
        for(int i=0;i<size; i++){
            auto curnode=d.nodes(i,start,end);
            auto dasi=S0(curnode);
            for(int j=0; j<size; j++){
                gsl_matrix_set(A,i,j,d.differentiationWeights(j,i,start,end));
                gsl_matrix_set(F,i,j,dasi*interpolation::beta(curnode,d.nodes(j,start,end))*d.quadratureWeights(j,start,end));
            }
            (*gsl_matrix_ptr(A,i,i))+=gamma(curnode)+mu(curnode);
        }
    }
    pair<double,double> computeSpectralRadius()
    {
        double a,r1=0,r2=0;
        // general eigenvalue problem
        gsl_eigen_gen_workspace *ws=gsl_eigen_gen_alloc(A->size1);
        gsl_vector_complex * alpha=gsl_vector_complex_alloc(A->size1);
        gsl_vector * beta=gsl_vector_alloc(A->size1);
        gsl_eigen_gen(F,A,alpha,beta,ws);
        for(int i=0; i<A->size1; i++){
            a=gsl_complex_abs(gsl_vector_complex_get(alpha,i))/abs(gsl_vector_get(beta,i));
            if(r1<a){
                r2=r1;
                r1=a;
            } else if (r2<a) {
                r2=a;
            }
        }
        gsl_vector_free(beta);
        gsl_vector_complex_free(alpha);
        gsl_eigen_gen_free(ws);

        // // standard eigenvalue problem
        // gsl_blas_dtrsm(CblasRight,  CblasLower,  CblasNoTrans, CblasNonUnit, 1, A, F);
        // gsl_eigen_nonsymm_workspace *ws2= gsl_eigen_nonsymm_alloc(A->size1);
        // gsl_eigen_nonsymm(F, alpha, ws2);
        // for(int i=0; i<A->size1; i++){
        //     a=gsl_complex_abs(gsl_vector_complex_get(alpha,i));
        //     if(radius<a)
        //         radius=a;
        // }
        // gsl_eigen_nonsymm_free(ws2);
        return make_pair(r1,r2);
    }
}
