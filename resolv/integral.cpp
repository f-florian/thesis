/***************************************************************************
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * Copyright (C) 2018 Francesco Florian
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ***************************************************************************/

constexpr double epsabs=1e-20;
constexpr double epsrel=1e-13;

#include "integral.h"
#include <cstdlib>
#include <cmath>

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include <stdexcept>
#include <queue>
#include <forward_list>

using namespace std;

static inline bool
subinterval_too_small (double a1, double a2, double b2)
{
  const double e = GSL_DBL_EPSILON;
  const double u = GSL_DBL_MIN;

  double tmp =  100 * e + 1000 * u;
  return (fabs (a2-a1) <= tmp && fabs (b2-a2) <= tmp);

}

namespace integrator {
    double qag (const gsl_function * f,
                const double a, const double b,
                const double epsabs, const double epsrel,
                const size_t limit,
                gsl_integration_workspace * workspace,
                double *result, double *abserr,
                gsl_integration_rule * q)
    {
        forward_list<array<double, 3>> points; //start, val, err
        struct arrcmp {
            static bool operator()(decltype(points.begin()) a, decltype(points.begin()) b){
                return(*a)[2]<(*b)[2];
            }
        };
        priority_queue<decltype(points.begin()), decltype(points.begin()), arrCmp> errcng;

        double area, errsum;
        double result, abserr, resabs, resasc;
        double tolerance;
        size_t iteration = 0;
        int roundoff_type1 = 0, roundoff_type2 = 0, error_type = 0;

        double round_off;     

        *result = 0;
        *abserr = 0;

        if (epsabs <= 0 && (epsrel < 50 * GSL_DBL_EPSILON || epsrel < 0.5e-28)) {
            throw(runtime_error("tolerance cannot be achieved with given epsabs and epsrel"));
        }
        
        /* perform the first integration */
        q(f, a, b, &result0, &abserr0, &resabs0, &resasc0);
        points.push_front({a, result, abserr});
        errcng.push(points.begin());
        
        /* Test on accuracy */
        tolerance = GSL_MAX_DBL (epsabs, epsrel * fabs (result));

        /* need IEEE rounding here to match original quadpack behavior */

        round_off = GSL_COERCE_DBL (50 * GSL_DBL_EPSILON * resabs0);

        if (abserr <= round_off && abserr > tolerance)
            throw(runtime_error("cannot reach tolerance because of roundoff error on first attempt"));
        else if ((abserr <= tolerance && abserr != resasc) || abserr == 0.0) {
                *result = result;
                *abserr = abserr;
                return GSL_SUCCESS;
        }
        
        area = result;
        errsum = abserr;

        for(iteration = 1; iteration < limit; ++iteration) {
            double start, end, middle;
            bool roundoff;
            double a_i, b_i, r_i, e_i;
            double area1 = 0, area2 = 0, area12 = 0;
            double error1 = 0, error2 = 0, error12 = 0;
            double resasc1, resasc2;
            double resabs1, resabs2;

            /* Bisect thise subinterval with the largest error estimate */
            auto cur=errcng.top();
            //pop?
            auto next=cur;
            ++next;
            
            start=(*cur)[0];
            if(next==errcng.end())
                end=b;
            else
                end=(*next)[0];
            middle=(start+end)/2;
            
            q (f, start, middle, &area, &error, &resabs, &resasc);
            roundoff=(resasc != error);
            (*cur)[1]=area;
            (*cur)[2]=error;


            q (f, middle, end, &area, &error, &resabs, &resasc);
            if (raundoff && resasc != error) {
                double delta = r_i - area12;
                
                if (fabs (delta) <= 1.0e-5 * fabs (area12) && error12 >= 0.99 * e_i)
                    {
                        roundoff_type1++;
                    }
                if (iteration >= 10 && error12 > e_i)
                    {
                        roundoff_type2++;
                    }
            }
            
            tolerance = GSL_MAX_DBL (epsabs, epsrel * fabs (area));
            
            if (errsum > tolerance)
                {
                    if (roundoff_type1 >= 6 || roundoff_type2 >= 20)
                        {
                            error_type = 2;   /* round off error */
                        }
                    
                    /* set error flag in the case of bad integrand behaviour at
                       a point of the integration range */
                    
                    if (subinterval_too_small (a1, a2, b2))
                        {
                            error_type = 3;
                        }
                }
            
            update (workspace, a1, b1, area1, error1, a2, b2, area2, error2);
            
            retrieve (workspace, &a_i, &b_i, &r_i, &e_i);
            
            iteration++;
            
        }
        while (iteration < limit && !error_type && errsum > tolerance);

        *result = sum_results (workspace);
        *abserr = errsum;

        if (errsum <= tolerance)
            {
                return GSL_SUCCESS;
            }
        else if (error_type == 2)
            {
                GSL_ERROR ("roundoff error prevents tolerance from being achieved",
                           GSL_EROUND);
            }
        else if (error_type == 3)
            {
                GSL_ERROR ("bad integrand behavior found in the integration interval",
                           GSL_ESING);
            }
        else if (iteration == limit)
            {
                GSL_ERROR ("maximum number of subdivisions reached", GSL_EMAXITER);
            }
        else
            {
                GSL_ERROR ("could not integrate function", GSL_EFAILED);
            }
    }


}


    type_m(f, start, end, &val, &err, &resabs, &resasc);
    if (val > maxabs || val>maxrel*(end-start)) {
        val=0;
        err=0;
    } else {
        val=exp(-val);
        err=expm1(err)*val;
    }
    points.push_front({err, val, start});
    
    tolerance = GSL_MAX_DBL (0,maxerror * fabs (val));
    
    /* need IEEE rounding here to match original quadpack behavior */
    // round_off = GSL_COERCE_DBL (50 * GSL_DBL_EPSILON * resabs);
    round_off = 0;
    
    if (points.front()[0] <= round_off && points.front()[0] > tolerance)
        throw(runtime_error("cannot reach tolerance because of roundoff error on first attempt"));

    decltype(points.begin()) cur=points.begin(),next;

    int i, depth=0;
    for (i = 0; i < limit; ++i) {
        if (cur==points.end())
            break;
        if ((*cur)[0]<maxerror) {
            depth=0;
            cur++;
            continue;
        }
        depth++;
        olde=(*cur)[0];
        oldv=(*cur)[1];
        next=cur;
        next++;
        double endint;
        if (next==points.end())
            endint=end_m;
        else
            endint=(*next)[2];
        double middle=((*cur)[2]+endint)/2;
        type_m(f, (*cur)[2], middle, &val, &err, &resabs, &resasc);
        bool differ=(resasc != err);
        if (val > maxabs || val>maxrel*(end-start)) {
            val=0;
            err=0;
        } else {
            val=exp(-val);
            err=expm1(err)*val;
        }
        (*cur)[0]=err;
        (*cur)[1]=val;
        type_m(f, middle, endint, &val, &err, &resabs, &resasc);
        differ=differ && (resasc != err);
        if (val > maxabs || val>maxrel*(end-start)) {
            val=0;
            err=0;
        } else {
            val=exp(-val);
            err=expm1(err)*val;
        }
        next=points.insert(next, {err, val, middle});
        
        newv=((*cur)[1]*(*next)[1]);
        newe=((*cur)[0]+(*next)[0]);

        if (differ && (olde>1e-11*oldv)) {
            if (fabs (oldv-newv) <= 1.0e-2 * olde) {
                fprintf(stderr, "rt1 %e %e\n", oldv-newv, olde);
                roundoff_type1++;
            }
            if (depth >= 4 &&  newe > olde) {
                roundoff_type2++;
            }
        }
        if (roundoff_type1 >= 6 || roundoff_type2 >= 20)
            throw runtime_error("roundoff error in integration");
            
        /*     set error flag in the case of bad integrand behaviour at
               a point of the integration range */
        
        if (subinterval_too_small ((*cur)[1], middle, endint)){
            // printf("fail en int %e %e %e\n", (*cur)[1],middle,endint);
            // break;
            throw runtime_error("Adaptive quadrature is trying to use intervals too small");
        }
    }
    double *x, *y;
    x=new double[i+3];
    y=new double[i+3];
    y[0]=1;
    x[i+2]=end_m;
    int j=0;
    for (auto it=points.begin(); it!=points.end(); ++it) {
        x[j]=(*it)[2];
        y[j+1]=(*it)[1];
        fprintf(stderr,"%2d %e %e\n",j,x[j],y[j+1]);
        ++j;
    }
    return {x,y,i+3};
}
