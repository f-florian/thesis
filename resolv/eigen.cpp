#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex_math.h>

#include "eigen.h"
#include "interpolation.h"

using namespace interpolation;

namespace eigen {
    namespace {
        gsl_matrix *A,*F;
        gsl_eigen_gen_workspace *ws;
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
        gsl_eigen_gen_free(ws);

    }
    void init(const unsigned int size, double start, double delta)
    {
        freemem();
        auto invd=1/delta;
        //todo: A sparse?
        //todo2: A already triangular!
        A=gsl_matrix_calloc(size,size);
        F=gsl_matrix_calloc(size,size);
        ws=gsl_eigen_gen_alloc(size);
        
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
    }
    double computeSpectralRadius()
    {
        // printf("init\n");
        // fflush(stdout);
        beta=gsl_vector_alloc(A->size1);
        alpha=gsl_vector_complex_alloc(A->size1);
        gsl_eigen_gen(F,A,alpha,beta,ws);
        // printf("calc\n");
        // fflush(stdout);
        double radius=0;
        for(int i=0; i<A->size1; i++){
            double a=gsl_complex_abs(gsl_vector_complex_get(alpha,i))/abs(gsl_vector_get(beta,i));
            if(radius<a)
                radius=a;
            // printf("eigenvalue %d=%lf\n", i,a);
            // fflush(stdout);
        }
        return radius;
    }
}
