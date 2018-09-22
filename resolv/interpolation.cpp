/***************************************************************************
 * Copyright (C) 2018 Francesco Florian
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ***************************************************************************/

#include "interpolation.h"
#include <gsl/gsl_spline.h>
#include <cstdio>


namespace interpolation
{
    namespace {
        double gamma0;
        double spreadFactor;            //3e-5
    }
    double beta(const double x1, const double x2){
        return spreadFactor*(6e-7*(10000-(x1-x2)*(x1-x2))+0.001);
    }
    double gamma(const double x){
        return gamma0;
    }
    double S0(const double x)
    {
        return 5187 + x * (226.438 - x * 2.777);
    }
    double mu(const double x)
    {
        return 8.3675/(110-x);
    }
    void splineReInit(const double param1, const double param2)
    {
        spreadFactor=param1;
        gamma0=param2;
    }
}
