#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex_math.h>
// #include <gsl/gsl_blas.h>

#include "eigen.h"
#include "interpolation.h"
#include "../differential/differential.h"

using namespace interpolation;

namespace eigen {
    namespace {
        gsl_matrix *A,*F;
        gsl_vector_complex * alpha;
	    gsl_vector * beta;
	    Differential* get;
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
	void setOrder(unsigned short ni)
	{
		get=new Differential(ni);
	}
	void init(const unsigned int size, double start, double delta, unsigned short order)
	{
		freemem();
		auto invd=1/delta;
		auto dim=size*order;
		//todo: A sparse?
		//todo2: A already triangular!
		A=gsl_matrix_calloc(dim,dim);
		F=gsl_matrix_calloc(dim,dim);
		
        for(int i=0;i<size-1; i++)
	        for(int l=0;l<order;l++){
		        // node : coordinate of m-th node in static submesh   
		        auto curnode=start+delta*(i+get->nodes(l));
		        // change invd with weigths for FD folmulas
		        gsl_matrix_set(A,i,i,-invd-gamma(curnode)-mu(curnode));
		        gsl_matrix_set(A,i+1,i,invd);
		        auto dasi=delta*S0(curnode);
		        for(int j=0; j<size; j++)
			        for(int m=0;m<order;m++){
				        //weightfq quadratura; remember to scale (ie: is dasi ok? do it in Differential?)
				        gsl_matrix_set(F,i,j*order+m,dasi*interpolation::beta(curnode,start+delta*j)*get->quadratureWeights(m));
			        }
	        }
        beta=gsl_vector_alloc(dim);
        alpha=gsl_vector_complex_alloc(dim);
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
