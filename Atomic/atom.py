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
SCF atomic calculation.
"""

import sys
import math
import numpy
import os

from periodic_table import Element, parse_configuration, tuple_to_configuration
from radial_grid import RadialGrid
from solve import solve_atom
from xc import LDA


class Orbital:
    def __init__(self, n, l, occ, npoints):
        self.n, self.l, self.occ = n, l, occ
        self.ene = 0.0
        self.psi = numpy.zeros(npoints)

class DFTAtom:
    def __init__(self, symbol, configuration = None, xc = LDA, 
                 maxiter = 200, mixing = 0.4, tol = 1e-6, output = None, 
                 write = None):
        """
        DFT atomic calculation.

        Parameters:
        symbol          chemical symbol
        configuration   '[He] 2s2 2p2', overrides default from table
        xc              XC functional (LDA default)
        """
        # store parameters
        self.element = Element(symbol)
        if configuration == None:
            self.configuration = self.element.get_configuration()
        else:
            self.configuration = configuration
        self.xc = xc
        self.maxiter = maxiter
        self.mixing = mixing
        self.tol = tol
        self.write = write
        self.set_output(output)

        # initialize grid, orbitals, energies
        self.grid = RadialGrid(self.element.get_atomic_number())
        self.setup_orbitals()
        
        # print header
        self.print_header()


    def setup_orbitals(self):
        """Initialize atomic orbitals"""
        self.config_tuple = parse_configuration(self.configuration)
        self.orbitals = []
        self.nelec = 0.0
        for orb in self.config_tuple:
            self.orbitals.append(Orbital(orb[0], orb[1], orb[2], self.grid.npoints))
            self.nelec += orb[2]

         
    def set_output(self, output):
        """Open output file"""
        if output == None:
            self.out = open(os.devnull, 'w')
        elif output == '-':
            self.out = sys.stdout
        else:
            self.out = open(output, 'a')



    def print_header(self):
        """Print header information"""
        print >> self.out, '='*78
        print >> self.out, "DFT atomic calculation: Z = %i, name = %s" % (self.element.get_atomic_number(), self.element.get_name())
        print >> self.out, '='*78
        print >> self.out, "Radial grid          : rmin = %f, rmax = %f, npoints = %i" % (self.grid.rmin, self.grid.rmax, self.grid.npoints)
        print >> self.out, "Atomic configuration :", tuple_to_configuration(self.config_tuple)
        print >> self.out, "Number of electrons  :", self.nelec
        print >> self.out, "Exchange-correlation :", self.xc.full_name
        print >> self.out, "SCF                  : maxiter = %i, mixing = %f, tol = %f" % (self.maxiter, self.mixing, self.tol)
        print >> self.out, ""



    def print_eigenvalues(self):
        """Print eigenvalues and occupations"""
        print >> self.out, "Eigenvalues:"
        print >> self.out, "    n l occ    energy (Ha)  energy (eV)"
        for orb in self.orbitals:
             print >> self.out, "    %i %s %-6g %12.6f %12.6f" % (orb.n, orb.l, orb.occ, orb.ene, orb.ene*27.211383)
        print >> self.out, ""


    
    def print_energies(self):
        """Print final energies"""
        print >> self.out, "Total energy    : %12.4f Ha" % (self.e_tot)
        print >> self.out, "Kinetic energy  : %12.4f Ha" % (self.e_kin)
        print >> self.out, "Ionic energy    : %12.4f Ha" % (self.e_ion)
        print >> self.out, "Hartree energy  : %12.4f Ha" % (self.e_h)
        print >> self.out, "XC energy       : %12.4f Ha" % (self.e_xc)



    def run_scf(self):
        """Run the SCF calculation"""
        # initialize potential (pure coulomb)
        npoints = self.grid.npoints
        self.v_ion = -self.grid.zeta/self.grid.r
        self.v_pot = self.v_ion.copy()

        # start SCF iteration
        rho_old = numpy.zeros(npoints)
        e_old = 0.0
        for it in xrange(1, self.maxiter):
            # solve hydrogenoic problem for every orbital
            for orb in self.orbitals:
                 ene, psi = solve_atom(self.grid, orb.n, orb.l, vpot = self.v_pot, maxiter = 50, tol = self.tol)
                 if ene != None:
                      orb.ene = ene
                      orb.psi = psi.copy()

            # rho mixing
            self.calculate_rho()
            if it > 1:
                self.rho = self.mixing*self.rho + (1.0 - self.mixing)*self.rho_old
            self.rho_old = self.rho.copy()

            # calculate potential
            self.calculate_v_of_rho()
            self.v_pot = self.v_ion + self.v_h + self.v_xc

            # calculate energy
            self.e_band = 0.0
            for orb in self.orbitals:
                self.e_band += orb.occ * orb.ene
            self.e_tot = self.e_band - self.e_h + self.e_xc - self.e_vxc
            self.e_kin = self.e_tot - self.e_ion - self.e_h - self.e_xc
            de_tot = self.e_tot - e_old
            e_old = self.e_tot

            if it == 1:
                print >> self.out, "Iteration     Total energy  Delta Energy"
            print >> self.out, "%6i        %12.6f  %12.4e" % (it, self.e_tot, de_tot)

            # check convergence and exit
            if abs(de_tot) < self.tol:
                break

        # end of SCF, print energy and eigenvalues
        print >> self.out, ""
        self.print_energies()
        self.print_eigenvalues()



    def calculate_rho(self):
        """Calculate the charge density"""
        self.rho = numpy.zeros(self.grid.npoints)
        for orb in self.orbitals:
             self.rho += orb.occ * (orb.psi*orb.psi) / self.grid.r2
        nelec = self.grid.integrate_two(self.rho, self.grid.r2)
        if abs(nelec - self.nelec) > 1e-4:
             print >> self.out, "integrated charge (%f) different from number of electrons (%f)" % (nelec, self.nelec)
	self.rho /= (4.0*math.pi)


    
    def calculate_hartree(self, rho):
        """Calculate the Hartree potential from a given density"""
        npoints = self.grid.npoints
        r, r2, dr, dx = self.grid.r, self.grid.r2, self.grid.dr, self.grid.dx

        # running charge
        q = numpy.zeros(npoints)
        for i in xrange(1,npoints):
	    q[i] = q[i-1] + 4.0*math.pi*r2[i] * dx*dr[i] * rho[i]

        # hartree
        v_h = numpy.zeros(npoints)
        v_h[-1] = q[-1]/r[-1]
        for i in xrange(npoints-1,0,-1):
	    v_h[i-1] = v_h[i] + dx * q[i]*self.grid.dr[i]/self.grid.r2[i]
        e_h = 0.5*self.grid.integrate_two(v_h, rho * r2) * 4.0*math.pi
        return e_h, v_h



    def calculate_xc(self, rho):
        """Calculate the XC potential from a given density"""
        npoints = self.grid.npoints
        r, r2, dr, dx = self.grid.r, self.grid.r2, self.grid.dr, self.grid.dx

        e_x = numpy.zeros(npoints)
        v_x = numpy.zeros(npoints)
        e_c = numpy.zeros(npoints)
        v_c = numpy.zeros(npoints)

        for i in xrange(len(rho)):
            (e_x[i], v_x[i]), (e_c[i], v_c[i]) = self.xc.xc(rho[i])

        v_xc = v_x + v_c            
        e_xc = self.grid.integrate_two(e_x + e_c, rho * r2) * 4.0*math.pi
        e_vxc = self.grid.integrate_two(v_xc, rho * r2) * 4.0*math.pi
        return e_xc, v_xc, e_vxc



    def calculate_v_of_rho(self):
        """Calculate Hartree and XC, potential and energy"""
        self.e_h, self.v_h = self.calculate_hartree(self.rho)
        self.e_ion = self.grid.integrate_two(self.v_ion, self.rho*self.grid.r2) * 4.0*math.pi
        self.e_xc, self.v_xc, self.e_vxc = self.calculate_xc(self.rho)


    def __getstate__(self):
	"""In order to implement pickle"""
	odict = self.__dict__.copy()
	del odict['out']
	return odict

        

if __name__ == '__main__':
    atom = DFTAtom("O", output = "-")
    atom.run_scf()

    quit()

    import pylab

    pylab.figure(1)
    pylab.plot(atom.grid.r, atom.rho*atom.grid.r2)
    pylab.xlim(0.0, 3.0)

    pylab.figure(2)
    pylab.plot(atom.grid.r, atom.v_pot)
    pylab.plot(atom.grid.r, atom.v_ion, '-')
    pylab.xlim(0.0, 3.0)
    pylab.ylim(-20.0, 0.0)

    pylab.figure(3)
    for orb in atom.orbitals:
        pylab.plot(atom.grid.r, orb.psi, label="n=%i, l=%i, ene=%f" % (orb.n, orb.l, orb.ene))
    pylab.xlim(0.0, 3.0)
    pylab.legend()

    pylab.show()
    quit()








