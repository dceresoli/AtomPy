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
Shooting in python.
"""
import math, numpy

def shoot(N, u, dx, c2, c1, c0):
    """
    Integrate differential equation by shooting method

        u''(x)*c2(x) + u'(x)*c1(x) + u(x)*c0(x) = 0

    on an equispaced grid (spacing dx).

    u[0:2] should be set in input according to boundary conditions.

    Return: u(x), number of nodes, turning point and discontinuity of
    derivative at the turning point
    """

    fp = c2 + 0.5*dx*c1
    fm = c2 - 0.5*dx*c1
    f0 = (dx*dx)*c0 - 2.0*c2

    # inward integration up to turning point (or one point beyond to get
    # derivative). If no turning point, integrate half-way
    u[-1] = 1.0
    u[-2] = -u[-1]*f0[-1]/fm[-1]
    all_positive = numpy.all(c0 > 0.0)
    for i in xrange(N-2, 0, -1):
        u[i-1] = (-fp[i]*u[i+1] - f0[i]*u[i]) / fm[i]
        if abs(u[i-1]) > 1e10:
            u[i-1:] *= 1e-10 # for numerical stability
        if c0[i] < 0.0:
            turn = i
            break
        if all_positive and i == N/2:
            turn = N/2
            break

    # right derivative
    uturn = u[turn]
    uturn1 = u[turn+1]
    dright = (u[turn+1] - u[turn-1]) / (2.0*dx)

    # outward integration up to turning point (one point beyond to get
    # derivative)
    for i in xrange(1, turn+1):
        u[i+1] = (-f0[i]*u[i] - fm[i]*u[i-1]) / fp[i]

    # left derivative
    dleft = (u[turn+1] - u[turn-1]) / (2.0*dx)

    # rescale u in order to make it continuos
    scale = uturn / u[turn]
    u[:turn+1] *= scale
    u[turn+1] = uturn1 # above overrode
    dleft *= scale
   
    # set the sign
    u = u*numpy.sign(u[1])

    # count nodes
    nodes = numpy.sum( (u[0:turn-1]*u[1:turn]) < 0.0 )
    derdisc = (dright - dleft) * uturn
    return u, nodes, turn, derdisc

