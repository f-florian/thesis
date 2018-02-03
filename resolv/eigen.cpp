#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex_math.h>
// #include <gsl/gsl_blas.h>

#include "eigen.h"
#include "interpolation.h"

using namespace interpolation;

namespace eigen {
    namespace {
        gsl_matrix *A,*F;
        gsl_vector_complex * alpha;
        gsl_vector * beta;
    }
    void freemem()
    {
        if(A!=NULL)
            gsl_matrix_free(A);
        if(F!=NULL)
            gsl_matrix_free(F);
        if(beta!=NULL)
            gsl_vector_free(beta);
        if(alpha!=NULL)
            gsl_vector_complex_free(alpha);
    }
	void init(const unsigned int size, double start, double delta)
	{
		freemem();
		auto invd=1/delta;
		//todo: A sparse?
		//todo2: A already triangular!
		A=gsl_matrix_calloc(size,size);
		F=gsl_matrix_calloc(size,size);
		
        for(int i=0;i<size-1; i++) {
            auto curnode=start+delta*i;
            gsl_matrix_set(A,i,i,-invd-gamma(curnode)-mu(curnode));
            gsl_matrix_set(A,i+1,i,invd);
            auto dasi=delta*S0(curnode);
            for(int j=0; j<size; j++){
                gsl_matrix_set(F,i,j,dasi*interpolation::beta(curnode,start+delta*j));
            }
        }
        auto curnode=start+delta*(size-1);
        gsl_matrix_set(A,size-1,size-1,-invd-gamma(curnode)-mu(curnode));
        auto dasi=delta*S0(curnode);
        for(int j=0; j<size; j++){
            gsl_matrix_set(F,size-1,j,dasi*interpolation::beta(curnode,start+delta*j));
        }
        beta=gsl_vector_alloc(size);
        alpha=gsl_vector_complex_alloc(size);
    }
    double computeSpectralRadius()
    {
	    double a,radius=0;
	    // general eigenvalue problem
        gsl_eigen_gen_workspace *ws=gsl_eigen_gen_alloc(A->size1);
        gsl_eigen_gen(F,A,alpha,beta,ws);
        for(int i=0; i<A->size1; i++){
            a=gsl_complex_abs(gsl_vector_complex_get(alpha,i))/abs(gsl_vector_get(beta,i));
            if(radius<a)
                radius=a;
            // printf("eigenvalue %d=%lf\n", i,a);
            // fflush(stdout);
        }
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
        return radius;
    }
}
