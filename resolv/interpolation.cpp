#include "interpolation.h"
#include <gsl/gsl_spline.h>
#include <cstdio>


namespace interpolation
{
    namespace {
        bool analytic;
        gsl_interp_accel *accels;
        gsl_interp_accel *accelm;
        gsl_spline *splinemu;
        gsl_spline *splineS0;
    }
    const double k=3e-5;
    double beta(const double x1, const double x2){
        return k*(6e-7*(10000-(x1-x2)*(x1-x2))+0.001);
    }
    double gamma(const double x){
        return 52;
    }
    double S0(const double x)
    {
        if (analytic)
            return 5187+ x*(226.438-2.776583333333333-x*2.776583333333333);
        if(x>100){
            fprintf(stderr,"Requested s interpolation in %.2e, which is %.2e beyond the maximum value\n",x,x-100);
            return gsl_spline_eval(splineS0,100,accels);
        }
        if(x<0){
            fprintf(stderr,"Requested s interpolation in %.2e, which is %.2e below the minimum value\n",x,x);
            return gsl_spline_eval(splineS0,100,accels);
        }
        return gsl_spline_eval(splineS0,x,accels);
    }
    double mu(const double x)
    {
        if (analytic)
            return 8.3675/(110-x);
        if(x>100){
            fprintf(stderr,"Requested interpolation in %.2e, which is %.2e beyond the maximum value\n",x,x-100);
            return gsl_spline_eval(splinemu,100,accelm);
        }
        if(x<0){
            fprintf(stderr,"Requested mu interpolation in %.2e, which is %.2e below the minimum value\n",x,x);
            return gsl_spline_eval(splineS0,100,accels);
        }
        return gsl_spline_eval(splinemu,x,accelm);
    }
    void splineInit(const int sizes, const int sizem, const double mx[], const double my[], const double sx[], const double sy[], bool analytic_p)
    {
        analytic=analytic_p;
        free();
        accels=gsl_interp_accel_alloc();
        accelm=gsl_interp_accel_alloc();
        splinemu=gsl_spline_alloc(gsl_interp_cspline, sizem);
        splineS0=gsl_spline_alloc(gsl_interp_cspline, sizes);
        gsl_spline_init(splinemu, mx, my, sizem);
        gsl_spline_init(splineS0, sx, sy, sizes);
    }
    void free()
    {
        if(accels!=NULL)
            gsl_interp_accel_free(accels);
        if(accelm!=NULL)
            gsl_interp_accel_free(accelm);
    }
}
