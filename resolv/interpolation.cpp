#include "interpolation.h"
#include <gsl/gsl_spline.h>
#include <cstdio>


namespace interpolation
{
    namespace {
        gsl_interp_accel *accels;
        gsl_interp_accel *accelm;
        gsl_spline *splinemu;
        gsl_spline *splineS0;
    }
    const double k=3e-5;
    double beta(double x1, double x2){
	    return k*(6e-7*(10000-(x1-x2)*(x1-x2))+0.001);
    }
    double gamma(double x){
        return 52;
    }
    double S0(double x)
    {
        if((100.-x)<0.00001){
	        printf("anoeinaos\n");
	        return 65;
        }
        printf("S interpolation in %le", x);
        fflush(stdout);
        auto tmp=gsl_spline_eval(splineS0,x,accels);
        printf(" (%le)\n", tmp);
        fflush(stdout);
        return tmp;
    }
    double mu(double x)
    {
        printf("mu interpolation in %le", x);
        fflush(stdout);
        auto tmp=gsl_spline_eval(splinemu,x,accelm);
        printf(" (%le)\n", tmp);
        fflush(stdout);
        return tmp;
    }
    void splineInit(int sizes, int sizem, const double mx[], const double my[], const double sx[], const double sy[])
    {
	if(accels!=NULL)
	    gsl_interp_accel_free(accels);
	if(accelm!=NULL)
	    gsl_interp_accel_free(accelm);
	accels=gsl_interp_accel_alloc();
	accelm=gsl_interp_accel_alloc();
	splinemu=gsl_spline_alloc(gsl_interp_cspline, sizem);
	splineS0=gsl_spline_alloc(gsl_interp_cspline, sizes);
	gsl_spline_init(splinemu, mx, my, sizem);
	gsl_spline_init(splineS0, sx, sy, sizes);
    }
}
