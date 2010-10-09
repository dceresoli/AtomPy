# AtomPy - atomic calcalations in python
# Copyright (C) 2010  Davide Ceresoli <dceresoli@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
Logarithmic grid.
"""

from math import log
from numpy import linspace, exp, sum

class RadialGrid:
    def __init__(self, zeta=1.0, rmin=0.0001, rmax=100.0, npoints=1000):
        """
        Initialize a radial logarithmic grid, for nuclear charge zeta.
        """
        # store parameters
        self.zeta = zeta
        self.npoints = npoints
        # construct logaritmic mesh
        self.rmin = rmin
        self.rmax = rmax
        xmin, xmax = log(zeta*rmin), log(zeta*rmax)
        self.x = linspace(xmin, xmax, npoints)
        self.dx = self.x[1] - self.x[0]
        self.r = exp(self.x)/zeta
        self.dr = self.r
        self.r2 = self.r * self.r

    def integrate_one(self, f):
        """Integrate a function of the radial grid (trapezium rule)."""
        res = f[0]*self.dr[0] + f[-1]*self.dr[-1]
        res += 2.0*sum(f[1:-1]*self.dr[1:-1])
        return res * self.dx/2.0

    def integrate_two(self, f1, f2):
        """Integrate the product of two functions (trapezium rule)."""
        res = f1[0]*f2[0]*self.dr[0] + f1[-1]*f2[-1]*self.dr[-1]
        res += 2.0*sum(f1[1:-1]*f2[1:-1]*self.dr[1:-1])
        return res * self.dx/2.0

