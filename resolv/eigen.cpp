#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex_math.h>
// #include <gsl/gsl_blas.h>

#include "eigen.h"
#include "interpolation.h"
#include "../differential/differential.h"

#define HAVE_INLINE
#define GSL_RANGE_CHECK_OFF

using namespace interpolation;

namespace eigen {
	Differential* get;
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
	void setOrder(unsigned short ni)
	{
		get=new Differential(ni);
	}
	void init(const unsigned int size, double start, double delta, unsigned short order)
	{
		freemem();
		auto invd=1/delta;
		auto dim=size*order;
		// allocation
		A=gsl_matrix_calloc(dim,dim);
		F=gsl_matrix_alloc(dim,dim);
		beta=gsl_vector_alloc(dim);
		alpha=gsl_vector_complex_alloc(dim);
		
		//todo: A sparse?
		//todo2: A already triangular!
		for(int i=0;i<size; i++)
			for(int l=0;l<order;l++){
				auto curnode=start+delta*(i+get->nodes(l));
				for(int m=0;m<order;m++)
					gsl_matrix_set(A,i*order+l,i*order+m,get->differentiationWeights(m,l)*invd);
				(*gsl_matrix_ptr(A,i*order+l,i*order+l))-=gamma(curnode)+mu(curnode);
				auto dasi=delta*S0(curnode);
				for(int j=0; j<size; j++)
					for(int m=0;m<order;m++){
						// printf("F idx:%d, node:%lf, weight:%lf\n", j*order+m, get->nodes(m), get->quadratureWeights(m));
						gsl_matrix_set(F,i*order+l,j*order+m,dasi*interpolation::beta(curnode,start+delta*(j+get->nodes(m)))*get->quadratureWeights(m));
					}
			}

		// for(int i=0;i<dim;i++){
		// 	for(int j=0;j<dim;j++)
		// 		printf("%9.2le ",gsl_matrix_get(F,i,j));
		// 	printf("\n");
		// }
			
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
