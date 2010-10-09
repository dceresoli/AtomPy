#!/usr/bin/env python
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

from Atomic import *

# parse command line
from optparse import OptionParser

parser = OptionParser("usage: %prog atom.data")
(options, args) = parser.parse_args()

if len(args) == 0:
    parser.print_help()
    exit(1)

# read data
import cPickle as pickle
atom = pickle.load(file(args[0]))

# plot
import pylab

pylab.figure(1)
pylab.plot(atom.grid.r, atom.rho*atom.grid.r2, label='charge density')
pylab.xlim(0.0, 3.0)
pylab.legend()

pylab.figure(2)
pylab.plot(atom.grid.r, atom.v_pot, label='SCF potential')
pylab.plot(atom.grid.r, atom.v_ion, label='ionic potential')
pylab.xlim(0.0, 3.0)
pylab.ylim(-20.0, 0.0)
pylab.legend()

pylab.figure(3)
for orb in atom.orbitals:
    pylab.plot(atom.grid.r, orb.psi, label="n=%i, l=%i, ene=%f" % (orb.n, orb.l, orb.ene))
    pylab.xlim(0.0, 3.0)
    pylab.legend()

pylab.show()




