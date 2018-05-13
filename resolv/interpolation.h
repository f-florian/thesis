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

#ifndef INTERPOLATION_H
#define INTERPOLATION_H

namespace interpolation
{
    void splineInit(const int size, const int sizem, const double mx[], const double my[], const double sx[], const double sy[], bool analytic);
    double mu(const double x);
    double gamma(const double x);
    double beta(const double x1, const double x2);
    double S0(const double x);
    void free();
}

#endif // INTERPOLATION_H
