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
Solve the radial Schroedinger equation by the shooting method.
"""

import math, numpy

try:
    from clib.shoot import shoot
except:
    from shoot import shoot
    import sys
    sys.stderr.write("using shoot.py instead of shoot.pyx\n")


def solve_atom(grid, n = 1, l = 0, vpot = None, maxiter = 100, tol = 1e-6):
    """
    Solve the radial Schrodinger equation by shooting and bisection.

    Parameters:
    grid         (in)    radial_grid objects
    n, l         (in)    principal and azimutal quantum numbers (n > l)
    vpot         (in)    external potential (if None, assume Coulomb potential)
    maxiter      (in)    maximum number of iterations
    tol          (in)    eigenvalue convergence

    Return:
    ene                  eigenvalue
    psi                  radial wave function
    """

    # initialization
    zeta, dx = grid.zeta, grid.dx
    r, x, r2 = grid.r, grid.x, grid.r2
    x = grid.x
    dx = grid.dx
    assert 0 <= l < n
    expected_nodes = n - l - 1

    # setup potential
    if vpot == None:
        vpot = -zeta/r

    # add angular momentum
    vpot += l*(l+1)/(2.0*r2)

    # initial guess
    emin, emax = -zeta*zeta/(n*n), 0.0
    ene = (emin + emax)/2.0

    # start iteration
    iteration = 0
    while True:
        iteration += 1

        # setup coefficients
        c0 = 2.0*r2*(vpot - ene)
        c1 = numpy.ones_like(x)
        c2 = -numpy.ones_like(x)

        # setup boundary conditions
        u = numpy.zeros_like(x)
        u[0:2] = r[0:2]**(l+1)

        # integrate and normalize
        u, nodes, turn, derdisc = shoot(grid.npoints, u, dx, c2, c1, c0)
        norm = grid.integrate_two(u, u)
        u /= math.sqrt(norm)
        ##print ">>>Iteration %3i: (%i,%i)   Energy = %12.6f   Nodes = %i" % (iteration, n, l, ene, nodes)

        # check max number of iterations, if not converged, return None as an energy
        if iteration > maxiter:
            ##print "wavefunction (%i,%i) not converged" % (n, l)
            return None, u

        # bisection on number of nodes
        if nodes > expected_nodes:
            emax = ene
            ene = (emin + emax)/2.0
            continue
        elif nodes < expected_nodes:
            emin = ene
            ene = (emin + emax)/2.0
            continue
 
        # refine energy
        shift = -0.5 * derdisc / (r[turn]*norm)

        # if correction is too large, go back to refine number of nodes
        if abs(shift) > zeta*zeta/(n*n) or ene + shift > 0.0:
            if shift > 0.0:
                emin = ene
            else:
                emax = ene
            ene = (emin + emax)/2.0
            continue

        ene += shift
        ene = min(ene, emax)
        ene = max(ene, emin)

        # check for convergence
        if (abs(shift) < tol):
            return ene, u



