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
Hydrogenoic wavefunctions.
"""

import math
import numpy

# helper functions
def binomial(n, m):
    return math.factorial(n)//(math.factorial(n-m)*math.factorial(m))

def gen_laguerre(n, alpha, x):
    res = numpy.zeros_like(x)
    for i in xrange(n+1):
       res += (-1)**i * binomial(n+alpha,n-i) * (x**i) / math.factorial(i)
    return res


# hydrogenoic wave function
def hydrogen_wfc(grid, n, l):
    # energy
    ene = -0.5*grid.zeta**2/n**2

    # wave function
    rho = 2.0*grid.r*grid.zeta/n
    psi = numpy.exp(-rho/2) * (rho**(l+1)) * gen_laguerre(n-l-1, 2*l+1, rho)
    norm = grid.integrate_two(psi, psi)
    psi /= math.sqrt(norm)

    return ene, psi

