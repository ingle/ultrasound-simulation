'''A program to be run from the command line.  It reads in the phantom file, and displays a backscatter coefficient curve
along with some information about the phantom.  Then it computes Faran backscatter coefficients and tries to 
match your own backscatter curve to one of the faran curves. '''

import sys
import numpy as np
from matplotlib import pyplot
from faranScattering import faranBsc
f = open(sys.argv[1])

phanSizeX =float( np.fromfile(f, np.double,1) )
phanSizeY = float( np.fromfile(f, np.double,1) )
phanSizeZ = float( np.fromfile(f, np.double, 1) )
totalScatters = int(np.fromfile(f, np.int32,1))
c0 = float(np.fromfile(f, np.double, 1))
a0 = float(np.fromfile(f, np.double,1))
a1 = float(np.fromfile(f, np.double,1))
a2 = float(np.fromfile(f, np.double,1))


print 'Sound speed of: ' + str(c0)
print 'Attenuation of: ' + str(a0) + ' + ' + str(a1) +'^' + str(a2)
print 'Phantom size of: ' + str(phanSizeX) + ', ' + str(phanSizeY) + ', ' + str(phanSizeZ)
print "Total number of scatterers is: " + str(totalScatters)

scatterers = np.fromfile(f, np.double, 3*totalScatters)
numBsc = int(np.fromfile(f, np.int32, 1))
freqStep = float(np.fromfile(f, np.double, 1))
bsc = np.fromfile(f, np.double, numBsc)
f.close()
freq = np.arange(0, numBsc)*freqStep

pyplot.plot(freq, bsc)
pyplot.show()
