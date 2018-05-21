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

struct gsl_integration_rule;
struct gsl_function;

namespace integrator {
    double qag (const gsl_function * f,
                const double a, const double b,
                const double epsabs, const double epsrel,
                const size_t limit,
                gsl_integration_workspace * workspace,
                double *result, double *abserr,
                gsl_integration_rule * q);
}
