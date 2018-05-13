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

#ifndef EIGEN_H
#define EIGEN_H

#include <gsl/gsl_matrix.h>
#include "differential.h"
#include <utility>
#include <array>

namespace eigen
{
    double computeSpectralRadius(double hint);
    void init(const size_t size, const Differential::Type type, const double points[], const size_t npts);
    template<size_t n> void init(const size_t size, const Differential::Type type, const std::array<double,n> points) {return init(size,type,points.data(),n);}
    void freemem();
}

#endif // EIGEN_H
