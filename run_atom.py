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

parser = OptionParser("usage: %prog [options] atomic_symbol (e.g. C)")
parser.add_option("-o", "--output", dest="output", help="save data to file", metavar="FILE", default="")
parser.add_option("-c", "--config", dest="config", help="electronic configuration (e.g. \"[He] 2s 2p\")", default="")

(options, args) = parser.parse_args()

if len(args) == 0:
    parser.print_help()
    exit(1)

output = options.output
config = options.config

# run SCF calculation
el = Element(args[0])
if config == '':
    config = el.get_configuration()

atom = DFTAtom(args[0], output='-', configuration=config)
atom.run_scf()

if output != '':
    import cPickle as pickle
    f = file(output, "wb")
    pickle.dump(atom, f)
    print "data for plotting saved to file:", output


