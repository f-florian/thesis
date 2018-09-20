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

#include <cstdlib>
#include <cstdio>

#include <gsl/gsl_vector.h>
#include "eigen.h"
#include "interpolation.h"
#include "differential.h"

constexpr double agemax=1;
constexpr Differential::Type inttype=Differential::Type::ClenshawCurtis;

constexpr size_t orderp1=100;

double* mainf(size_t size, double* var1, double var2)
{
    double *out=static_cast<double*>(malloc(size*sizeof(double)));
    eigen::alloc(orderp1, ittype);
    
    for(size_t i = 0; i < size; ++i) {
        eigen::init();
        *(out+i)=eigen::computeSpectralRadius();
    }
    eigen::freemem();
    interpolation::free();
    return out;
}
