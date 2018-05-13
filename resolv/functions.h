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

namespace parameters{
    void init(const double length_p, const double c0_p, const double D0_p, const double beta0_p, const double mu0_p, const double c1_p, const double D1_p, const int p_p);
    double c(const double a);
    double D(const double a);
    double beta(const double a);
    double betaonc(const double a);
    double mu(const double a);
    double length();
}
